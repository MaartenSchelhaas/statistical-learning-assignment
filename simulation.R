
# --------------------------------------------------------------------
# Author:
# *Enter your group number here, as well as names and student numbers*
# --------------------------------------------------------------------
# Simulation study: standard kNN (weighted = FALSE) vs. hybrid kNN with
# random-forest importance weights (weighted = TRUE).
#
# Run from the folder containing hybrid_kNN.R:  Rscript simulation.R
# Packages: ranger (used by hybrid_kNN), ggplot2 (plots), parallel (optional).
# --------------------------------------------------------------------

source("hybrid_kNN.R")
library(ggplot2)

## ====================================================================
## 0. Simulation parameters
## ====================================================================
n          <- 500                 # observations
p          <- 10                  # variables
R          <- 200                 # replications
k          <- 5                   # neighbours in the main study
miss_vars  <- c("X1", "X2", "X3") # variables that receive missing values
miss_rates <- c(0.10, 0.30)       # target missing proportion per variable
mechanisms <- c("MCAR", "MAR", "MNAR")
scenarios  <- c("A", "B", "C")    # drop "C" to skip the nonlinear scenario

rho_A      <- 0.5                 # Scenario A: equicorrelation
rho_B      <- 0.8                 # Scenario B: correlation among X1-X3
sd_eps_C   <- 0.5                 # Scenario C: noise sd in X1 = sin(X2) + X3^2 + e
mar_driver <- "X4"                # fully observed variable driving MAR
slope      <- 2                   # logistic slope on the standardized driver

num_trees  <- 500                 # trees per random forest
seed_base  <- 20261009            # base for all seeds
n_cores    <- 1                   # > 1 uses parallel::mclapply (not on Windows)

run_k_check <- TRUE               # optional check with several k
k_grid      <- c(1, 5, 10)
R_k         <- 50                 # replications for the k check
rate_k      <- 0.30               # missing rate for the k check

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)


## ====================================================================
## 1. Data generating processes
## ====================================================================
# Multivariate normal draws via the Cholesky factor of the covariance matrix.
rmvn <- function(n, Sigma) {
  matrix(rnorm(n * ncol(Sigma)), n) %*% chol(Sigma)
}

generate_data <- function(scenario, n, p) {
  if (scenario == "A") {
    # Every variable equally correlated with every other: all are equally
    # informative, so equal weights should be (close to) optimal.
    Sigma <- matrix(rho_A, p, p)
    diag(Sigma) <- 1
    X <- rmvn(n, Sigma)
  } else if (scenario == "B") {
    # X1-X3 strongly correlated block; X4-X10 independent noise. Equal weights
    # let the 7 noise variables dominate the distance; RF weights should not.
    Sigma <- diag(p)
    Sigma[1:3, 1:3] <- rho_B
    diag(Sigma) <- 1
    X <- rmvn(n, Sigma)
  } else if (scenario == "C") {
    # Nonlinear: X1 = sin(X2) + X3^2 + e, all other variables independent N(0,1).
    # X1 is (almost) uncorrelated with X2 and X3, but a random forest can
    # still detect the dependence.
    X <- matrix(rnorm(n * p), n, p)
    X[, 1] <- sin(X[, 2]) + X[, 3]^2 + rnorm(n, sd = sd_eps_C)
  } else {
    stop("unknown scenario")
  }
  colnames(X) <- paste0("X", seq_len(p))
  X
}


## ====================================================================
## 2. Missing data mechanisms (implemented by hand)
## ====================================================================
# Logistic missingness model P(M = 1) = plogis(a + slope * z), where z is the
# standardized driver. The intercept a is calibrated with uniroot() so that the
# average missingness probability over the sample equals the target rate, i.e.
# the expected missing proportion equals the target.
calibrate_intercept <- function(z, rate, slope) {
  uniroot(function(a) mean(plogis(a + slope * z)) - rate,
          interval = c(-10, 10), extendInt = "yes", tol = 1e-10)$root
}

# Returns a logical matrix (TRUE = missing). Missingness is generated
# independently for each variable in 'vars'.
#   MCAR: independent Bernoulli(rate)
#   MAR : depends on the fully observed variable 'driver'
#   MNAR: depends on the variable's own (unobserved) value
make_missing <- function(X, vars, rate, mechanism, driver = mar_driver,
                         slope = 2) {
  M <- matrix(FALSE, nrow(X), ncol(X), dimnames = dimnames(X))
  for (v in vars) {
    if (mechanism == "MCAR") {
      prob <- rep(rate, nrow(X))
    } else {
      z <- if (mechanism == "MAR") X[, driver] else X[, v]
      z <- (z - mean(z)) / sd(z)
      a <- calibrate_intercept(z, rate, slope)
      prob <- plogis(a + slope * z)
    }
    M[, v] <- runif(nrow(X)) < prob
  }
  M
}


## ====================================================================
## 3. Evaluation: RMSE on the imputed cells only
## ====================================================================
# kNN with Euclidean distances and mean aggregation targets the conditional
# mean of the missing value, i.e. it minimizes squared error, so RMSE between
# imputed and true values on the imputed cells is the matching criterion.
rmse_imputed <- function(X_true, X_hat, M, vars) {
  sapply(vars, function(v) {
    m <- M[, v]
    if (!any(m)) return(NA_real_)
    sqrt(mean((X_hat[m, v] - X_true[m, v])^2))
  })
}


## ====================================================================
## 4. One replication of one scenario
## ====================================================================
# The same complete data set is used for all mechanisms and missing rates in
# a replication, and the same incomplete data set is given to both methods
# (paired comparison). All randomness is driven by the seed for (scenario, r),
# so results do not depend on n_cores.
run_one_rep <- function(r, scenario, k_values = k, rates = miss_rates,
                        mechs = mechanisms, store_weights = FALSE) {
  sc_id <- match(scenario, LETTERS)
  set.seed(seed_base + 10000 * sc_id + r)
  X <- generate_data(scenario, n, p)
  
  res <- list(); prop <- list(); wts <- list()
  for (mech in mechs) for (rate in rates) {
    M <- make_missing(X, miss_vars, rate, mech, slope = slope)
    X_mis <- X
    X_mis[M] <- NA
    
    # realized missing proportions
    prop[[length(prop) + 1]] <- data.frame(
      scenario = scenario, rep = r, mechanism = mech, rate = rate,
      variable = miss_vars, realized = colMeans(M[, miss_vars]))
    
    for (kk in k_values) for (wtd in c(FALSE, TRUE)) {
      fit <- hybrid_kNN(X_mis, k = kk, weighted = wtd, seed = r,
                        num.trees = num_trees)
      rmse <- rmse_imputed(X, fit$X_hat, M, miss_vars)
      res[[length(res) + 1]] <- data.frame(
        scenario = scenario, rep = r, mechanism = mech, rate = rate, k = kk,
        method = if (wtd) "hybrid" else "standard",
        variable = c(miss_vars, "average"),
        rmse = c(rmse, mean(rmse)),
        n_mean_fallback = c(fit$diagnostics$n_mean_fallback[
          match(miss_vars, fit$diagnostics$variable)], NA))
      
      if (store_weights && wtd) {
        for (v in names(fit$weights)) {
          wts[[length(wts) + 1]] <- data.frame(
            scenario = scenario, rep = r, mechanism = mech, rate = rate,
            target = v, predictor = names(fit$weights[[v]]),
            weight = unname(fit$weights[[v]]))
        }
      }
    }
  }
  list(rmse = do.call(rbind, res), prop = do.call(rbind, prop),
       weights = if (length(wts)) do.call(rbind, wts) else NULL)
}

# apply over replications, optionally in parallel
run_reps <- function(reps, FUN, ...) {
  if (n_cores > 1 && .Platform$OS.type == "unix") {
    parallel::mclapply(reps, FUN, ..., mc.cores = n_cores)
  } else {
    lapply(reps, FUN, ...)
  }
}
bind <- function(lst, el) do.call(rbind, lapply(lst, `[[`, el))


## ====================================================================
## 5. Main study
## ====================================================================
all_runs <- list()
for (sc in scenarios) {
  t0 <- Sys.time()
  cat("Scenario", sc, "...")
  all_runs[[sc]] <- run_reps(seq_len(R), run_one_rep, scenario = sc,
                             store_weights = (sc == "B"))
  cat(" done in", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
}
runs <- unlist(all_runs, recursive = FALSE)
results  <- bind(runs, "rmse")
realized <- bind(runs, "prop")
weights_B <- bind(runs, "weights")

saveRDS(results,   file.path(out_dir, "rmse_results.rds"))
saveRDS(realized,  file.path(out_dir, "realized_missing.rds"))
saveRDS(weights_B, file.path(out_dir, "weights_scenarioB.rds"))


## ====================================================================
## 6. Summaries
## ====================================================================
mc_se <- function(x) sd(x) / sqrt(length(x))
key <- c("scenario", "mechanism", "rate", "variable")

# 6a. Realized missing proportions (target vs. realized)
prop_summary <- aggregate(realized ~ scenario + mechanism + rate + variable,
                          data = realized,
                          FUN = function(x) c(mean = mean(x), sd = sd(x),
                                              min = min(x), max = max(x)))
prop_summary <- do.call(data.frame, prop_summary)
write.csv(prop_summary, file.path(out_dir, "realized_missing_summary.csv"),
          row.names = FALSE)
cat("\nRealized missing proportions (mean over replications):\n")
print(xtabs(realized.mean ~ mechanism + rate + scenario,
            data = aggregate(realized.mean ~ scenario + mechanism + rate,
                             data = prop_summary, FUN = mean)), digits = 3)

# 6b. Mean RMSE with Monte Carlo SE per method
rmse_summary <- aggregate(rmse ~ scenario + mechanism + rate + variable + method,
                          data = results,
                          FUN = function(x) c(mean = mean(x), se = mc_se(x)))
rmse_summary <- do.call(data.frame, rmse_summary)

# 6c. Paired differences hybrid - standard (negative = hybrid better)
wide <- reshape(results[, c(key, "rep", "method", "rmse")],
                idvar = c(key, "rep"), timevar = "method", direction = "wide")
wide$diff <- wide$rmse.hybrid - wide$rmse.standard
diff_summary <- aggregate(diff ~ scenario + mechanism + rate + variable,
                          data = wide,
                          FUN = function(x) c(mean = mean(x), se = mc_se(x),
                                              hybrid_better = mean(x < 0)))
diff_summary <- do.call(data.frame, diff_summary)

# one table for the report: average over X1-X3
tab <- merge(
  reshape(rmse_summary[rmse_summary$variable == "average",
                       c("scenario", "mechanism", "rate", "method",
                         "rmse.mean", "rmse.se")],
          idvar = c("scenario", "mechanism", "rate"), timevar = "method",
          direction = "wide"),
  diff_summary[diff_summary$variable == "average", -4])
tab <- tab[order(tab$scenario, match(tab$mechanism, mechanisms), tab$rate), ]
write.csv(rmse_summary, file.path(out_dir, "rmse_summary_by_variable.csv"),
          row.names = FALSE)
write.csv(diff_summary, file.path(out_dir, "rmse_difference_summary.csv"),
          row.names = FALSE)
write.csv(tab, file.path(out_dir, "rmse_table_average.csv"), row.names = FALSE)
cat("\nAverage RMSE over X1-X3 (mean and MC SE), paired difference hybrid - standard:\n")
print(format(tab, digits = 3), row.names = FALSE)

# 6d. How often did the mean fallback occur? (see note in hybrid_kNN.R)
fb <- aggregate(n_mean_fallback ~ scenario + mechanism + rate + method,
                data = results[results$variable != "average", ], FUN = sum)
fb <- fb[fb$n_mean_fallback > 0, ]
cat("\nRecipients imputed by observed-mean fallback (total over reps and X1-X3):\n")
if (nrow(fb)) print(fb, row.names = FALSE) else cat("none\n")
write.csv(fb, file.path(out_dir, "mean_fallback_counts.csv"), row.names = FALSE)

# 6e. RF weights in Scenario B: share of weight on the informative block
if (!is.null(weights_B)) {
  weights_B$informative <- weights_B$predictor %in% miss_vars
  w_share <- aggregate(weight ~ mechanism + rate + target + rep, FUN = sum,
                       data = weights_B[weights_B$informative, ])
  w_share <- aggregate(weight ~ mechanism + rate + target, data = w_share,
                       FUN = function(x) c(mean = mean(x), se = mc_se(x)))
  w_share <- do.call(data.frame, w_share)
  names(w_share)[4:5] <- c("share_X1toX3.mean", "share_X1toX3.se")
  write.csv(w_share, file.path(out_dir, "weights_B_informative_share.csv"),
            row.names = FALSE)
  cat("\nScenario B: total RF weight on the informative variables X1-X3",
      "(equal weights would give 2/9 = 0.222):\n")
  print(format(w_share, digits = 3), row.names = FALSE)
}


## ====================================================================
## 7. Plots
## ====================================================================
mech_lab <- function(d) {
  d$mechanism <- factor(d$mechanism, levels = mechanisms)
  d$rate <- factor(paste0(100 * d$rate, "% missing"),
                   levels = paste0(100 * miss_rates, "% missing"))
  d
}

# 7a. Paired RMSE differences (average over X1-X3)
d_avg <- mech_lab(wide[wide$variable == "average", ])
p_diff <- ggplot(d_avg, aes(x = mechanism, y = diff, fill = rate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.size = 0.6) +
  facet_wrap(~ scenario, scales = "free_y",
             labeller = labeller(scenario = function(s) paste("Scenario", s))) +
  labs(x = NULL, fill = NULL,
       y = "RMSE hybrid - RMSE standard\n(negative = hybrid better)") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "rmse_difference_boxplots.pdf"), p_diff,
       width = 8, height = 4)
ggsave(file.path(out_dir, "rmse_difference_boxplots.png"), p_diff,
       width = 8, height = 4, dpi = 200)

# 7b. Same, per variable (for an appendix)
d_var <- mech_lab(wide[wide$variable != "average", ])
p_diff_var <- p_diff %+% d_var +
  facet_grid(variable ~ scenario, scales = "free_y")
ggsave(file.path(out_dir, "rmse_difference_boxplots_by_variable.pdf"),
       p_diff_var, width = 8, height = 7)

# 7c. RF weights in Scenario B (mean weight per predictor and target)
if (!is.null(weights_B)) {
  w_mean <- aggregate(weight ~ mechanism + rate + target + predictor,
                      data = weights_B, FUN = mean)
  w_mean <- mech_lab(w_mean)
  w_mean$predictor <- factor(w_mean$predictor, levels = paste0("X", 1:p))
  p_w <- ggplot(w_mean, aes(x = predictor, y = weight, fill = mechanism)) +
    geom_col(position = "dodge") +
    geom_hline(yintercept = 1 / (p - 1), linetype = "dashed") +
    facet_grid(rate ~ target,
               labeller = labeller(target = function(s) paste("Imputing", s))) +
    labs(x = "Predictor", y = "Mean RF weight",
         caption = "Dashed line: equal weight 1/(p-1)") +
    theme_bw() + theme(legend.position = "bottom",
                       axis.text.x = element_text(angle = 90, vjust = 0.5))
  ggsave(file.path(out_dir, "weights_scenarioB.pdf"), p_w, width = 9, height = 5)
  ggsave(file.path(out_dir, "weights_scenarioB.png"), p_w, width = 9, height = 5,
         dpi = 200)
}


## ====================================================================
## 8. Optional: sensitivity to k
## ====================================================================
if (run_k_check) {
  cat("\nk check (k =", paste(k_grid, collapse = ", "), "; R =", R_k, ") ...\n")
  k_runs <- list()
  for (sc in scenarios) {
    k_runs[[sc]] <- run_reps(seq_len(R_k), run_one_rep, scenario = sc,
                             k_values = k_grid, rates = rate_k)
  }
  k_res <- bind(unlist(k_runs, recursive = FALSE), "rmse")
  k_res <- k_res[k_res$variable == "average", ]
  saveRDS(k_res, file.path(out_dir, "k_check_results.rds"))
  
  k_sum <- aggregate(rmse ~ scenario + mechanism + k + method, data = k_res,
                     FUN = function(x) c(mean = mean(x), se = mc_se(x)))
  k_sum <- do.call(data.frame, k_sum)
  write.csv(k_sum, file.path(out_dir, "k_check_summary.csv"), row.names = FALSE)
  print(format(k_sum, digits = 3), row.names = FALSE)
  
  k_sum$mechanism <- factor(k_sum$mechanism, levels = mechanisms)
  p_k <- ggplot(k_sum, aes(x = k, y = rmse.mean, colour = method)) +
    geom_line() + geom_point() +
    geom_errorbar(aes(ymin = rmse.mean - 2 * rmse.se,
                      ymax = rmse.mean + 2 * rmse.se), width = 0.3) +
    facet_grid(scenario ~ mechanism, scales = "free_y") +
    scale_x_continuous(breaks = k_grid) +
    labs(y = "Mean RMSE (average over X1-X3)",
         caption = paste0(100 * rate_k, "% missing; error bars: +/- 2 MC SE")) +
    theme_bw() + theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "k_check.pdf"), p_k, width = 8, height = 6)
}

saveRDS(sessionInfo(), file.path(out_dir, "sessionInfo.rds"))
cat("\nAll results saved in '", out_dir, "/'.\n", sep = "")

