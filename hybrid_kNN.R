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
  
  ## ------------------------------------------------------------------
  ## 0. Input checks and preparation
  ## ------------------------------------------------------------------
  
  is_df <- is.data.frame(X)
  if (!is_df && !is.matrix(X)) {
    stop("'X' must be a matrix or a data.frame")
  }
  if (is_df) {
    # only continuous variables are supported
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
  n <- nrow(Xm)
  p <- ncol(Xm)
  
  # assign variable names V1, V2, ... if missing (also fill in empty names and
  # make duplicated names unique, since the weights list is indexed by name)
  cn <- colnames(Xm)
  if (is.null(cn)) cn <- rep("", p)
  empty <- is.na(cn) | cn == ""
  cn[empty] <- paste0("V", seq_len(p))[empty]
  cn <- make.unique(cn)
  colnames(Xm) <- cn
  
  # check k
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1) {
    stop("'k' must be a single positive integer")
  }
  k <- as.integer(floor(k))
  
  weighted <- isTRUE(weighted)
  if (weighted && !requireNamespace("ranger", quietly = TRUE)) {
    stop("package 'ranger' is required for weighted = TRUE")
  }
  
  # helper to return the result in the same type as the input
  finalize <- function(X_imp, weights, importance, diagnostics) {
    if (is_df) {
      X_imp <- as.data.frame(X_imp, stringsAsFactors = FALSE)
      names(X_imp) <- cn
      attr(X_imp, "row.names") <- row_names
    }
    list(X_hat = X_imp, weights = weights, importance = importance,
         diagnostics = diagnostics, k = k, weighted = weighted)
  }
  
  miss <- is.na(Xm)                    # missingness indicator (NaN counts as NA)
  n_miss <- colSums(miss)
  X_hat <- Xm                          # imputations are written into this copy
  weights <- list()
  importance <- list()
  diagnostics <- data.frame(variable = character(0), n_missing = integer(0),
                            n_fewer_than_k = integer(0), n_mean_fallback = integer(0),
                            stringsAsFactors = FALSE)
  
  # nothing to impute
  if (n == 0L || p == 0L || !any(miss)) {
    return(finalize(X_hat, weights, importance, diagnostics))
  }
  
  ## ------------------------------------------------------------------
  ## 1. Standardization (for distance computation only)
  ## ------------------------------------------------------------------
  # Each variable is centered and scaled with the mean and standard deviation
  # of its observed values, so that variables measured on different scales
  # contribute comparably to the distance and the weights alone determine
  # their relative influence. Constant columns (sd = 0) and columns with a
  # single observed value (sd = NA) are only centered; they then contribute
  # zero to every distance. Fully missing columns stay NA.
  mu <- colMeans(Xm, na.rm = TRUE)                      # NaN for fully missing
  s  <- apply(Xm, 2, stats::sd, na.rm = TRUE)           # NA if < 2 observations
  s[is.na(s) | s == 0] <- 1
  Z <- sweep(sweep(Xm, 2, mu, "-"), 2, s, "/")
  
  # Z with NAs replaced by 0, and observation indicator as numeric matrix;
  # this allows the pairwise-complete distances to be computed with matrix
  # products (see step 3)
  Z0 <- Z
  Z0[miss] <- 0
  Obs <- (!miss) * 1
  
  # medians of observed values, used as temporary fill for RF predictors only
  med <- apply(Xm, 2, stats::median, na.rm = TRUE)      # NA for fully missing
  
  ## ------------------------------------------------------------------
  ## 2.-4. Loop over variables with missing values
  ## ------------------------------------------------------------------
  for (j in which(n_miss > 0)) {
    
    var_j <- cn[j]
    rec <- which(miss[, j])            # recipients: rows where j is missing
    don <- which(!miss[, j])           # donors: rows where j is observed
    others <- setdiff(seq_len(p), j)   # candidate distance variables
    
    # fully missing column: nothing to learn from; warn and leave as NA
    if (length(don) == 0L) {
      warning("variable '", var_j, "' is completely missing and is left as NA")
      next
    }
    
    # p = 1: there are no other variables to compute distances on
    if (length(others) == 0L) {
      warning("only one variable available; '", var_j,
              "' is imputed by its observed mean")
      X_hat[rec, j] <- mean(Xm[don, j])
      weights[[var_j]] <- numeric(0)
      importance[var_j] <- list(NULL)
      diagnostics <- rbind(diagnostics, data.frame(
        variable = var_j, n_missing = length(rec), n_fewer_than_k = 0L,
        n_mean_fallback = length(rec), stringsAsFactors = FALSE))
      next
    }
    
    ## --- 2. Variable weights ------------------------------------------
    equal_w <- setNames(rep(1 / length(others), length(others)), cn[others])
    imp_raw <- NULL
    
    if (!weighted) {
      w <- equal_w                     # standard kNN
    } else {
      # predictors: all other variables; fully missing predictors cannot be
      # used by the forest and receive weight 0
      usable <- others[!is.na(med[others])]
      w <- setNames(rep(0, length(others)), cn[others])
      
      if (length(usable) == 0L || length(don) < 2L) {
        # no forest can be fitted; fall back to equal weights
        w <- equal_w
      } else {
        # training data: rows where j is observed; NAs in predictors get a
        # temporary median fill (used for fitting the forest only)
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
          warning("random forest for '", var_j, "' failed (", conditionMessage(e),
                  "); using equal weights")
          NULL
        })
        
        if (is.null(imp_raw)) {
          w <- equal_w
        } else {
          imp <- pmax(imp_raw, 0)                 # negative importance -> 0
          imp[!is.finite(imp)] <- 0
          if (sum(imp) > 0) {
            w[cn[usable]] <- imp / sum(imp)       # normalize to sum to 1
          } else {
            w <- equal_w                          # all zero -> equal weights
          }
        }
      }
    }
    weights[[var_j]] <- w
    importance[var_j] <- list(imp_raw)
    
    ## --- 3. Weighted distances (vectorized, blockwise over recipients) --
    # For recipient r and donor d, with shared observed set S(r,d) and weights
    # w_l renormalized over S(r,d), the squared weighted Euclidean distance is
    #   D2(r,d) = sum_{l in S} w_l (z_rl - z_dl)^2 / sum_{l in S} w_l.
    # Because Z0 is 0 wherever a value is missing, expanding the square gives
    #   numerator   = (Z0_r^2 w) Obs_d' - 2 (Z0_r w) Z0_d' + (Obs_r w) (Z0_d^2)'
    #   denominator = (Obs_r w) Obs_d'
    # where every product only involves variables observed in both rows.
    # Only the original (incomplete) data are used, never earlier imputations,
    # so the order in which variables are imputed does not matter.
    Z0_d  <- Z0[don, others, drop = FALSE]
    Obs_d <- Obs[don, others, drop = FALSE]
    Z0_d2 <- Z0_d^2
    y_don <- Xm[don, j]                            # donor values, original scale
    k_eff <- min(k, length(don))
    
    n_fewer <- 0L
    n_fallback <- 0L
    blocks <- split(rec, ceiling(seq_along(rec) / block_size))
    
    for (b in blocks) {
      Z0_r  <- Z0[b, others, drop = FALSE]
      Obs_r <- Obs[b, others, drop = FALSE]
      Z0_rw  <- sweep(Z0_r, 2, w, "*")
      Obs_rw <- sweep(Obs_r, 2, w, "*")
      
      num <- tcrossprod(Z0_rw * Z0_r, Obs_d) - 2 * tcrossprod(Z0_rw, Z0_d) +
        tcrossprod(Obs_rw, Z0_d2)
      den <- tcrossprod(Obs_rw, Obs_d)
      
      # donors without shared observed variables (or whose shared variables all
      # have weight 0, so the weights cannot be renormalized) are excluded
      D <- sqrt(pmax(num, 0) / den)              # pmax guards rounding errors
      D[den <= 0] <- Inf
      D <- matrix(D, nrow = length(b))           # keep matrix shape if 1 row
      
      ## --- 4. Aggregation: mean of j over the k nearest donors ----------
      for (i in seq_along(b)) {
        d_i <- D[i, ]
        valid <- which(is.finite(d_i))
        if (length(valid) == 0L) {
          # no valid donor at all: fall back to the observed mean of j
          X_hat[b[i], j] <- mu[j]
          n_fallback <- n_fallback + 1L
        } else {
          kk <- min(k_eff, length(valid))        # fewer than k -> use available
          if (kk < k) n_fewer <- n_fewer + 1L
          nn <- valid[order(d_i[valid])[seq_len(kk)]]   # ties: lower row index
          X_hat[b[i], j] <- mean(y_don[nn])
        }
      }
    }
    
    diagnostics <- rbind(diagnostics, data.frame(
      variable = var_j, n_missing = length(rec), n_fewer_than_k = n_fewer,
      n_mean_fallback = n_fallback, stringsAsFactors = FALSE))
  }
  
  rownames(diagnostics) <- NULL
  finalize(X_hat, weights, importance, diagnostics)
}

