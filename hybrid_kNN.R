# --------------------------------------------------------------------
# Author:
# *Enter your group number here, as well as names and student numbers*
# --------------------------------------------------------------------

## Use this code skeleton to implement the hybrid kNN algorithm.  Please submit
## your implementation together with your report and your simulation code.

## IMPORTANT: Please do not change the function name and make sure that you
##            use the correct input and output as specified below. This will
##            simplify grading because each group's code will be similarly
##            structured.  However, feel free to add other arguments to the
##            function definition as needed.


# Input:
# X .............. data matrix or data.frame with continuous (numeric) variables;
#                  missing values coded as NA
# k .............. number of neighbors (positive integer)
# weighted ....... logical; if TRUE (default), distances use random-forest
#                  permutation-importance weights (hybrid kNN); if FALSE, all
#                  variables get equal weights (standard kNN). Both variants
#                  share the same code path.
# seed ........... optional integer seed passed to ranger() for reproducible
#                  random forests (only relevant if weighted = TRUE)
# num.trees ...... number of trees in each random forest
# block_size ..... number of recipients processed at once in the vectorized
#                  distance computation (limits memory use for large n)
# ... ............ further arguments passed on to ranger::ranger()
#                  (e.g., mtry, min.node.size, num.threads)

# Output
# A list with the following components:
# X_hat ....... the imputed data set (same type as the input X)
# weights ..... a list of variable importance weights from the random forest
#               models. That is, each list element contains the weights for
#               the corresponding variable to be imputed. Use the variable name
#               as the name of the list element.
# Additional components:
# importance .. list of raw permutation importances (before truncation and
#               normalization; NULL if no forest was fitted)
# diagnostics . data.frame with, per imputed variable: number of missing
#               values, number of recipients imputed with fewer than k donors,
#               and number of recipients imputed by the observed-mean fallback
# k, weighted . the settings that were used

hybrid_kNN <- function(X, k, weighted = TRUE, seed = NULL, num.trees = 500,
                       block_size = 500, ...) {

  ## 0. Input checks and preparation
  prep <- check_and_prepare(X, k, weighted)
  Xm <- prep$Xm
  cn <- prep$cn
  k <- prep$k
  n <- nrow(Xm)
  p <- ncol(Xm)

  miss <- is.na(Xm)                    # NaN counts as NA
  n_miss <- colSums(miss)
  X_hat <- Xm                          # imputations are written into this copy
  weights <- list()
  importance <- list()
  diagnostics <- data.frame(variable = character(0), n_missing = integer(0),
                            n_fewer_than_k = integer(0), n_mean_fallback = integer(0),
                            stringsAsFactors = FALSE)

  # nothing to impute
  if (n == 0L || p == 0L || !any(miss)) {
    return(finalize_result(X_hat, weights, importance, diagnostics, prep))
  }

  ## 1. Standardization (for distance computation only)
  zdata <- standardize(Xm, miss)

  ## 2. Impute each variable x_j that has missing values
  for (j in which(n_miss > 0)) {

    name_j <- cn[j]
    rec <- which(miss[, j])            # recipients: rows where x_j is missing
    don <- which(!miss[, j])           # donors: rows where x_j is observed
    others <- setdiff(seq_len(p), j)   # variables used for the distances (x_-j)

    # fully missing column: nothing to learn from; warn and leave as NA
    if (length(don) == 0L) {
      warning("variable '", name_j, "' is completely missing and is left as NA")
      next
    }

    # p = 1: there are no other variables to compute distances on
    if (length(others) == 0L) {
      warning("only one variable available; '", name_j,
              "' is imputed by its observed mean")
      X_hat[rec, j] <- mean(Xm[don, j])
      weights[[name_j]] <- numeric(0)
      importance[name_j] <- list(NULL)
      diagnostics <- rbind(diagnostics,
                           diagnostics_row(name_j, length(rec), 0L, length(rec)))
      next
    }

    ## 2a. Variable weights
    vw <- variable_weights(j, Xm, don, others, cn, zdata$med, prep$weighted,
                           num.trees, seed, ...)
    weights[[name_j]] <- vw$w
    importance[name_j] <- list(vw$imp_raw)

    ## 2b. Weighted distances and aggregation over the k nearest donors
    imp <- impute_variable(j, rec, don, others, Xm, zdata, vw$w, k, block_size)
    X_hat[rec, j] <- imp$values
    diagnostics <- rbind(diagnostics,
                         diagnostics_row(name_j, length(rec), imp$n_fewer,
                                         imp$n_fallback))
  }

  rownames(diagnostics) <- NULL
  finalize_result(X_hat, weights, importance, diagnostics, prep)
}


## ====================================================================
## Helpers
## ====================================================================

# Input checks and column names. Returns the data as a double matrix with
# unique, non-empty column names, plus the validated k and weighted and what
# is needed to restore the input type at the end.
check_and_prepare <- function(X, k, weighted) {
  is_df <- is.data.frame(X)
  if (!is_df && !is.matrix(X)) {
    stop("'X' must be a matrix or a data.frame")
  }
  row_names <- NULL
  if (is_df) {
    is_num <- vapply(X, is.numeric, logical(1))
    if (!all(is_num)) {
      stop("all columns of 'X' must be numeric; non-numeric column(s): ",
           paste(which(!is_num), collapse = ", "))
    }
    row_names <- attr(X, "row.names")
    Xm <- as.matrix(X)
  } else {
    if (!is.numeric(X)) stop("'X' must be numeric")
    Xm <- X
  }
  storage.mode(Xm) <- "double"   # integer columns may receive non-integer imputations
  p <- ncol(Xm)

  # the weights list is indexed by name, so names must be present and unique
  cn <- colnames(Xm)
  if (is.null(cn)) cn <- rep("", p)
  empty <- is.na(cn) | cn == ""
  cn[empty] <- paste0("V", seq_len(p))[empty]
  cn <- make.unique(cn)
  colnames(Xm) <- cn

  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1) {
    stop("'k' must be a single positive integer")
  }
  k <- as.integer(floor(k))

  weighted <- isTRUE(weighted)
  if (weighted && !requireNamespace("ranger", quietly = TRUE)) {
    stop("package 'ranger' is required for weighted = TRUE")
  }

  list(Xm = Xm, cn = cn, k = k, weighted = weighted, is_df = is_df,
       row_names = row_names)
}


# Standardizes each column using its observed values (distances only).
# Z0 is the standardized data with NAs set to 0, Obs the 0/1 observed indicator.
standardize <- function(Xm, miss) {
  mu <- colMeans(Xm, na.rm = TRUE)                      # NaN for fully missing
  s  <- apply(Xm, 2, stats::sd, na.rm = TRUE)           # NA if < 2 observations
  s[is.na(s) | s == 0] <- 1                             # constant column: only center
  Z <- sweep(sweep(Xm, 2, mu, "-"), 2, s, "/")

  Z0 <- Z
  Z0[miss] <- 0
  Obs <- (!miss) * 1

  med <- apply(Xm, 2, stats::median, na.rm = TRUE)      # NA for fully missing

  list(mu = mu, Z0 = Z0, Obs = Obs, med = med)
}


# Weights of the variables x_-j for imputing x_j: normalized random-forest
# permutation importances, or equal weights for standard kNN and whenever no
# forest can be used. Also returns the raw importances (NULL without forest).
variable_weights <- function(j, Xm, don, others, cn, med, weighted, num.trees,
                             seed, ...) {
  equal_w <- setNames(rep(1 / length(others), length(others)), cn[others])
  if (!weighted) return(list(w = equal_w, imp_raw = NULL))

  # fully missing predictors cannot be used by the forest and get weight 0
  usable <- others[!is.na(med[others])]
  if (length(usable) == 0L || length(don) < 2L) {
    return(list(w = equal_w, imp_raw = NULL))
  }

  # train on the donors; NAs in predictors get a temporary median fill
  Xtr <- Xm[don, usable, drop = FALSE]
  for (l in seq_along(usable)) {
    na_l <- is.na(Xtr[, l])
    if (any(na_l)) Xtr[na_l, l] <- med[usable[l]]
  }
  Xtr <- as.data.frame(Xtr, stringsAsFactors = FALSE)
  names(Xtr) <- cn[usable]

  imp_raw <- tryCatch({
    rf <- ranger::ranger(x = Xtr, y = Xm[don, j], num.trees = num.trees,
                         importance = "permutation", seed = seed, ...)
    rf$variable.importance
  }, error = function(e) {
    warning("random forest for '", cn[j], "' failed (", conditionMessage(e),
            "); using equal weights")
    NULL
  })
  if (is.null(imp_raw)) return(list(w = equal_w, imp_raw = NULL))

  imp <- pmax(imp_raw, 0)                 # negative importance -> 0
  imp[!is.finite(imp)] <- 0
  if (sum(imp) == 0) return(list(w = equal_w, imp_raw = imp_raw))

  w <- setNames(rep(0, length(others)), cn[others])
  w[cn[usable]] <- imp / sum(imp)
  list(w = w, imp_raw = imp_raw)
}


# Imputed values of x_j for all recipients, processed in blocks of recipients
# to limit memory use. Also counts recipients imputed with fewer than k donors
# and recipients imputed by the observed mean.
impute_variable <- function(j, rec, don, others, Xm, zdata, w, k, block_size) {
  # donor side of the distance computation, shared by all recipient blocks
  Z0_d  <- zdata$Z0[don, others, drop = FALSE]
  Obs_d <- zdata$Obs[don, others, drop = FALSE]
  Z0_d2 <- Z0_d^2
  xj_don <- Xm[don, j]                            # donor values, original scale

  values <- numeric(length(rec))
  n_fewer <- 0L
  n_fallback <- 0L
  blocks <- split(seq_along(rec), ceiling(seq_along(rec) / block_size))

  for (pos in blocks) {
    b <- rec[pos]
    D <- weighted_distances(zdata$Z0[b, others, drop = FALSE],
                            zdata$Obs[b, others, drop = FALSE],
                            Z0_d, Obs_d, Z0_d2, w)

    for (i in seq_along(b)) {
      imp <- aggregate_donors(D[i, ], xj_don, k, zdata$mu[j])
      values[pos[i]] <- imp$value
      n_fewer <- n_fewer + imp$fewer_than_k
      n_fallback <- n_fallback + imp$mean_fallback
    }
  }

  list(values = values, n_fewer = n_fewer, n_fallback = n_fallback)
}


# Distances d(x_i, x_l) between a block of recipients and all donors, using the
# variables observed in both rows (S) with weights w renormalized over S:
#   d^2 = sum_S w_m (z_im - z_lm)^2 / sum_S w_m
# Z0 is 0 where missing, so the sums over S become matrix products with Obs.
# Inf where no shared variable has positive weight.
weighted_distances <- function(Z0_r, Obs_r, Z0_d, Obs_d, Z0_d2, w) {
  # multiply column m by its weight w_m
  Z0_r_w  <- sweep(Z0_r, 2, w, "*")
  Obs_r_w <- sweep(Obs_r, 2, w, "*")

  # numerator:
  sq_recipient <- tcrossprod(Z0_r_w * Z0_r, Obs_d)   # sum_m w_m z_im^2 Obs_lm
  cross        <- tcrossprod(Z0_r_w, Z0_d)           # sum_m w_m z_im z_lm
  sq_donor     <- tcrossprod(Obs_r_w, Z0_d2)         # sum_m w_m Obs_im z_lm^2
  num <- sq_recipient - 2 * cross + sq_donor

  # denominator:
  den <- tcrossprod(Obs_r_w, Obs_d)

  D <- sqrt(pmax(num, 0) / den)              # pmax guards rounding errors
  D[den <= 0] <- Inf
  D
}


# Imputed value for one recipient: mean of the donor values xj_don over the k
# nearest donors (d_i = distances to all donors). Uses fewer donors if fewer
# than k have a finite distance, and the observed mean mu_j if none has.
aggregate_donors <- function(d_i, xj_don, k, mu_j) {
  valid <- which(is.finite(d_i))
  if (length(valid) == 0L) {
    return(list(value = mu_j, fewer_than_k = FALSE, mean_fallback = TRUE))
  }
  kk <- min(k, length(valid))
  nn <- valid[order(d_i[valid])[seq_len(kk)]]   # ties: lower row index
  list(value = mean(xj_don[nn]), fewer_than_k = kk < k, mean_fallback = FALSE)
}


# One row of the diagnostics table.
diagnostics_row <- function(variable, n_missing, n_fewer_than_k, n_mean_fallback) {
  data.frame(variable = variable, n_missing = n_missing,
             n_fewer_than_k = n_fewer_than_k, n_mean_fallback = n_mean_fallback,
             stringsAsFactors = FALSE)
}


# Returns the result in the same type as the input.
finalize_result <- function(X_imp, weights, importance, diagnostics, prep) {
  if (prep$is_df) {
    X_imp <- as.data.frame(X_imp, stringsAsFactors = FALSE)
    names(X_imp) <- prep$cn
    attr(X_imp, "row.names") <- prep$row_names
  }
  list(X_hat = X_imp, weights = weights, importance = importance,
       diagnostics = diagnostics, k = prep$k, weighted = prep$weighted)
}
