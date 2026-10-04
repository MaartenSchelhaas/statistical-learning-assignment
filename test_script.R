# --------------------------------------------------------------------
# Tests for hybrid_kNN()
# Run from the folder that contains hybrid_kNN.R:  Rscript test_hybrid_kNN.R
# Uses base R only (plus ranger, which hybrid_kNN() needs).
# --------------------------------------------------------------------

source("hybrid_kNN.R")

## Small test helper: prints PASS/FAIL and counts failures -------------
n_fail <- 0L
check <- function(desc, cond) {
  ok <- isTRUE(cond)
  cat(if (ok) "PASS" else "FAIL", "-", desc, "\n")
  if (!ok) n_fail <<- n_fail + 1L
  invisible(ok)
}
tol <- 1e-10

## Independent reference: weighted distance between rows r and q for target
## column j, computed pair by pair (no matrix tricks), on data standardized by
## the observed mean/sd of each column. Returns the renormalized weights of the
## shared variables too, so their sum can be inspected.
pair_distance <- function(X, j, r, q, w, renormalize = TRUE) {
  X <- as.matrix(X)
  mu <- colMeans(X, na.rm = TRUE)
  s  <- apply(X, 2, sd, na.rm = TRUE)
  s[is.na(s) | s == 0] <- 1
  others <- setdiff(seq_len(ncol(X)), j)
  shared <- others[!is.na(X[r, others]) & !is.na(X[q, others])]
  if (length(shared) == 0L) return(list(d = Inf, w_shared = numeric(0)))
  w_sh <- w[colnames(X)[shared]]
  if (sum(w_sh) == 0) return(list(d = Inf, w_shared = w_sh))
  if (renormalize) w_sh <- w_sh / sum(w_sh)
  z_r <- (X[r, shared] - mu[shared]) / s[shared]
  z_q <- (X[q, shared] - mu[shared]) / s[shared]
  list(d = sqrt(sum(w_sh * (z_r - z_q)^2)), w_shared = w_sh)
}

## Reference imputation of one cell (row r, column j): mean of j over the k
## nearest donors according to pair_distance(); observed mean if no donor.
reference_impute <- function(X, j, r, k, w, renormalize = TRUE) {
  X <- as.matrix(X)
  don <- which(!is.na(X[, j]))
  d <- vapply(don, function(q) pair_distance(X, j, r, q, w, renormalize)$d, 0)
  valid <- which(is.finite(d))
  if (length(valid) == 0L) return(mean(X[don, j]))
  nn <- valid[order(d[valid])][seq_len(min(k, length(valid)))]
  mean(X[don[nn], j])
}


cat("\n== (1) Hand-made example ==\n")
# Target a is missing only in row 1. Distance variables b, c, d.
# Row 3 shares only b with row 1, all other donors share b, c and d.
X1 <- data.frame(
  a = c(NA, 10, 20, 30, 40, 50),
  b = c( 0, .6, .9,  2, -1,  1),
  c = c( 0, .6, NA,  2, -1, -2),
  d = c( 0, .6, NA,  2,  1, -2))
# Observed means / sds:  b: 0.5833 / 1.0128,  c: -0.08 / 1.5271,  d: 0.32 / 1.4873
# Standardized values (3 decimals):
#   row 1: b -0.576, c  0.052, d -0.215     (recipient)
#   row 2: b  0.016, c  0.445, d  0.188
#   row 3: b  0.313, c     NA, d     NA     (shares only b)
#   row 4: b  1.399, c  1.362, d  1.130
#   row 5: b -1.563, c -0.602, d  0.457
#   row 6: b  0.411, c -1.257, d -1.560
# Equal weights (weighted = FALSE): 1/3 each.
# Squared distances to row 1, WITH renormalization over shared variables:
#   row 2: 0.223 | row 3: 0.790 (weight of b renormalized to 1) | row 4: 2.474
#   row 5: 0.619 | row 6: 1.500
#   -> ranking 2, 5, 3, 6, 4 ;  k = 1 -> a = 10 ;  k = 2 -> (10 + 40)/2 = 25
# WITHOUT renormalization row 3 would get 0.790 / 3 = 0.263, i.e. it would look
# closer than row 5 only because it shares fewer variables:
#   -> ranking 2, 3, 5, ... ;  k = 2 -> (10 + 20)/2 = 15  (wrong)
r1_k1 <- hybrid_kNN(X1, k = 1, weighted = FALSE)
r1_k2 <- hybrid_kNN(X1, k = 2, weighted = FALSE)
r1_k3 <- hybrid_kNN(X1, k = 3, weighted = FALSE)
check("k = 1 imputes value of nearest donor (row 2): 10",
      abs(r1_k1$X_hat$a[1] - 10) < tol)
check("k = 2 imputes mean of rows 2 and 5: 25 (renormalized)",
      abs(r1_k2$X_hat$a[1] - 25) < tol)
check("k = 2 does NOT give 15 (what unnormalized weights would give)",
      abs(r1_k2$X_hat$a[1] - 15) > tol)
check("k = 3 imputes mean of rows 2, 5, 3: (10 + 40 + 20)/3",
      abs(r1_k3$X_hat$a[1] - 70/3) < tol)

w_eq <- r1_k2$weights$a
d_hand <- c(0.2227046, 0.7897303, 2.474496, 0.6186239, 1.49952)
d_ref  <- sapply(2:6, function(q) pair_distance(X1, 1, 1, q, w_eq)$d^2)
check("reference squared distances match the hand values above",
      max(abs(d_ref - d_hand)) < 1e-6)
check("unnormalized reference reproduces the wrong answer 15",
      abs(reference_impute(X1, 1, 1, 2, w_eq, renormalize = FALSE) - 15) < tol)
obs1 <- !is.na(X1)
check("observed values are unchanged",
      identical(as.matrix(r1_k2$X_hat)[obs1], as.matrix(X1)[obs1]))
check("all missing cells (a[1], c[3], d[3]) are imputed", !anyNA(r1_k2$X_hat))


cat("\n== (2) data.frame without column names ==\n")
set.seed(1)
X2 <- data.frame(matrix(rnorm(60), 20, 3))
X2[c(2, 5, 9), 1] <- NA
X2[c(3, 7), 2] <- NA
names(X2) <- NULL
r2 <- tryCatch(hybrid_kNN(X2, k = 3, seed = 1), error = function(e) e)
check("runs without error", !inherits(r2, "error"))
if (!inherits(r2, "error")) {
  check("returns a data.frame", is.data.frame(r2$X_hat))
  check("columns named V1, V2, V3", identical(names(r2$X_hat), c("V1", "V2", "V3")))
  check("weights list named by imputed variables V1, V2",
        identical(names(r2$weights), c("V1", "V2")))
  check("no NAs left", !anyNA(r2$X_hat))
}
X2b <- X2
names(X2b) <- c("x", "", "z")          # one empty name only
r2b <- hybrid_kNN(X2b, k = 3, seed = 1)
check("single empty name filled in as V2", identical(names(r2b$X_hat), c("x", "V2", "z")))
r2c <- hybrid_kNN(unname(as.matrix(X2)), k = 3, seed = 1)
check("unnamed matrix input returns a matrix with names V1..V3",
      is.matrix(r2c$X_hat) && identical(colnames(r2c$X_hat), c("V1", "V2", "V3")))


cat("\n== (3) Column without NAs ==\n")
r3 <- hybrid_kNN(X2b, k = 3, seed = 1)
check("complete column z is returned unchanged", identical(r3$X_hat$z, X2b$z))
check("complete column z has no element in weights", !("z" %in% names(r3$weights)))
check("complete column z has no row in diagnostics", !("z" %in% r3$diagnostics$variable))
check("complete column z is used as distance variable (has a weight)",
      "z" %in% names(r3$weights$x))
X3 <- as.matrix(X1[-c(1, 3), ])        # rows 1 and 3 hold all NAs of X1
stopifnot(!anyNA(X3))
check("data without any NA returned unchanged",
      identical(unname(hybrid_kNN(X3, k = 2)$X_hat), unname(X3)))


cat("\n== (4) Row with only one observed variable ==\n")
set.seed(2)
X4 <- data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30), d = rnorm(30))
X4$a <- X4$a + 2 * X4$b                # make b informative for a
X4[1, c("a", "c", "d")] <- NA          # row 1: only b observed
X4[c(4, 8, 15), "a"] <- NA
X4[c(6, 20), "c"] <- NA
for (wtd in c(FALSE, TRUE)) {
  r4 <- hybrid_kNN(X4, k = 3, weighted = wtd, seed = 1)
  lab <- paste0("(weighted = ", wtd, ") ")
  check(paste0(lab, "row 1 fully imputed"), !anyNA(r4$X_hat[1, ]))
  # Row 1 shares only b with any donor, so b's renormalized weight is 1
  # (if its weight is positive) and neighbours depend on b alone.
  for (v in c("a", "c", "d")) {
    j <- match(v, names(X4))
    expected <- reference_impute(X4, j, 1, 3, r4$weights[[v]])
    check(paste0(lab, v, "[1] equals mean of 3 nearest donors in b alone"),
          abs(r4$X_hat[1, v] - expected) < tol)
  }
}


cat("\n== (5) k larger than the number of donors ==\n")
X5 <- X1                                # 5 donors for a
r5 <- tryCatch(hybrid_kNN(X5, k = 100, weighted = FALSE), error = function(e) e)
check("runs without error", !inherits(r5, "error"))
if (!inherits(r5, "error")) {
  check("imputes the mean of all 5 donors: 30", abs(r5$X_hat$a[1] - 30) < tol)
  check("diagnostics report 1 recipient with fewer than k donors",
        r5$diagnostics$n_fewer_than_k[r5$diagnostics$variable == "a"] == 1)
}
# Donor without shared variables is excluded rather than counted:
X5b <- data.frame(a = c(NA, 1, 5, 9), b = c(1, NA, 2, 3), c = c(NA, 7, NA, NA))
r5b <- hybrid_kNN(X5b, k = 100, weighted = FALSE)
check("donor with no shared variable excluded: mean of rows 3, 4 = 7",
      abs(r5b$X_hat$a[1] - 7) < tol)


cat("\n== (6) weighted = FALSE gives equal weights ==\n")
set.seed(3)
X6 <- as.data.frame(matrix(rnorm(200), 40, 5))
X6[sample(40, 8), 1] <- NA
X6[sample(40, 8), 3] <- NA
r6 <- hybrid_kNN(X6, k = 5, weighted = FALSE)
check("weights for each imputed variable are all 1/(p-1) = 0.25",
      all(vapply(r6$weights, function(w) all(abs(w - 0.25) < tol), logical(1))))
check("weights are named by the 4 other variables",
      identical(names(r6$weights$V1), c("V2", "V3", "V4", "V5")))
check("no random forest fitted (importance is NULL)",
      all(vapply(r6$importance, is.null, logical(1))))
r6w <- hybrid_kNN(X6, k = 5, weighted = TRUE, seed = 1)
check("weighted = TRUE weights sum to 1 per variable",
      all(abs(vapply(r6w$weights, sum, 0) - 1) < tol))
check("weighted = TRUE weights are non-negative",
      all(unlist(r6w$weights) >= 0))
check("weighted = TRUE is reproducible with the same seed",
      identical(r6w$X_hat, hybrid_kNN(X6, k = 5, weighted = TRUE, seed = 1)$X_hat))


cat("\n== Shared-variable weights sum to 1 for recipient-donor pairs ==\n")
# hybrid_kNN() does not return per-pair weights, so this is checked in two
# steps: (a) for some pairs with different shared sets, the renormalized RF
# weights of the shared variables sum to 1; (b) hybrid_kNN() reproduces the
# reference imputations that use these renormalized weights, for every
# recipient, but not the ones based on unnormalized weights.
set.seed(4)
n <- 60
X7 <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n), x4 = rnorm(n))
X7$x1 <- X7$x1 + 1.5 * X7$x2 + X7$x3
for (v in names(X7)) X7[sample(n, 15), v] <- NA
r7 <- hybrid_kNN(X7, k = 4, weighted = TRUE, seed = 1)
w7 <- r7$weights$x1
cat("RF weights for x1:", round(w7, 3), "\n")

rec <- which(is.na(X7$x1))
don <- which(!is.na(X7$x1))
pairs <- expand.grid(r = rec[1:4], q = don[1:6])
sums <- t(apply(pairs, 1, function(pq) {
  pd <- pair_distance(X7, 1, pq["r"], pq["q"], w7)
  c(n_shared = length(pd$w_shared), sum_w = sum(pd$w_shared))
}))
print(cbind(pairs, round(sums, 12))[1:10, ])
has_info <- sums[, "n_shared"] > 0 & is.finite(sums[, "sum_w"])
check("renormalized shared weights sum to 1 for all checked pairs",
      all(abs(sums[has_info, "sum_w"] - 1) < tol))
check("checked pairs include different numbers of shared variables",
      length(unique(sums[has_info, "n_shared"])) > 1)

ref_renorm <- sapply(rec, function(r) reference_impute(X7, 1, r, 4, w7, TRUE))
ref_raw    <- sapply(rec, function(r) reference_impute(X7, 1, r, 4, w7, FALSE))
check("hybrid_kNN() equals renormalized reference for all recipients of x1",
      max(abs(r7$X_hat$x1[rec] - ref_renorm)) < tol)
check("unnormalized reference differs for at least one recipient (renormalization matters here)",
      any(abs(ref_raw - ref_renorm) > tol))
all_vars_ok <- all(vapply(names(r7$weights), function(v) {
  j <- match(v, names(X7)); rr <- which(is.na(X7[[v]]))
  ref <- sapply(rr, function(r) reference_impute(X7, j, r, 4, r7$weights[[v]]))
  max(abs(r7$X_hat[[v]][rr] - ref)) < tol
}, logical(1)))
check("same holds for all four imputed variables", all_vars_ok)


cat("\n==", if (n_fail == 0) "ALL TESTS PASSED" else paste(n_fail, "TEST(S) FAILED"), "==\n")


