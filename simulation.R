# --------------------------------------------------------------------
# Author:
# *Enter your group number here, as well as names and student numbers*
# --------------------------------------------------------------------
# Simulation study: standard kNN (weighted = FALSE) vs. hybrid kNN with
# random-forest importance weights (weighted = TRUE).
#
# Run from the folder containing hybrid_kNN.R:  Rscript simulation.R
# Packages: ranger (used by hybrid_kNN), ggplot2 (plots).
# --------------------------------------------------------------------

source("hybrid_kNN.R")
library(ggplot2)

## ====================================================================
## 0. Simulation parameters
## ====================================================================
n          <- 500                 # observations
p          <- 10                  # variables
R          <- 200                 # replications
k          <- 5                   # number of neighbours
miss_vars  <- c("X1", "X2", "X3") # variables that receive missing values
miss_rate  <- 0.30                # target missing proportion per variable
mechanisms <- c("MCAR", "MAR", "MNAR")
scenarios  <- c("A", "B")

rho_A      <- 0.5                 # Scenario A: equicorrelation
rho_B      <- 0.8                 # Scenario B: correlation among X1-X3
mar_driver <- "X4"                # fully observed variable driving MAR
slope      <- 2                   # logistic slope on the standardized driver

num_trees  <- 500                 # trees per random forest
seed_base  <- 20261009            # base for all seeds

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)


## ====================================================================
## 1. Seeds
## ====================================================================
# Every random step gets its own seed, computed only from the scenario label,
# the replication number, the mechanism and the missing rate -- never from
# the position in a loop. Adding or removing scenarios, mechanisms or rates
# therefore does not change the results of the remaining settings.
#   data:        seed_base + 10000 * scenario_id + r
#   missingness: data seed + 1e6 * mechanism_id + 1e7 * (100 * rate)
#   forests:     seed = r (passed to ranger via hybrid_kNN)
scenario_id  <- function(scenario) match(scenario, LETTERS)       # A = 1, B = 2
data_seed    <- function(scenario, r) seed_base + 10000 * scenario_id(scenario) + r
missing_seed <- function(scenario, r, mechanism, rate) {
  data_seed(scenario, r) + 1e6 * match(mechanism, c("MCAR", "MAR", "MNAR")) +
    1e7 * round(100 * rate)
}


## ====================================================================
## 2. Data generating processes
## ====================================================================
# Multivariate normal draws via the Cholesky factor of the covariance matrix.
rmvn <- function(n, Sigma) {
  matrix(rnorm(n * ncol(Sigma)), n) %*% chol(Sigma)
}

generate_data <- function(scenario, n, p) {
  Sigma <- diag(p)
  if (scenario == "A") {
    # Every variable equally correlated with every other: all are equally
    # informative, so equal weights should be (close to) optimal.
    Sigma[] <- rho_A
  } else if (scenario == "B") {
    # X1-X3 strongly correlated block; X4-X10 independent noise. Equal weights
    # let the 7 noise variables dominate the distance; RF weights should not.
    Sigma[1:3, 1:3] <- rho_B
  } else {
    stop("unknown scenario")
  }
  diag(Sigma) <- 1
  X <- rmvn(n, Sigma)
  colnames(X) <- paste0("X", seq_len(p))
  X
}


## ====================================================================
## 3. Missing data mechanisms (implemented by hand)
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
## 4. Evaluation: RMSE on the imputed cells only
## ====================================================================
# kNN with Euclidean distances and mean aggregation targets the conditional
# mean of the missing value, i.e. it minimizes squared error, so RMSE between
# imputed and true values on the imputed cells is the matching criterion.
# Per variable: RMSE over that variable's imputed cells.
# Combined ("average"): unweighted mean of the three per-variable RMSEs.
rmse_imputed <- function(X_true, X_hat, M, vars) {
  sapply(vars, function(v) {
    m <- M[, v]
    if (!any(m)) return(NA_real_)
    sqrt(mean((X_hat[m, v] - X_true[m, v])^2))
  })
}


## ====================================================================
## 5. One replication of one scenario
## ====================================================================
# The same complete data set is used for all mechanisms in a replication, and
# the same incomplete data set is given to both methods (paired comparison).
run_one_rep <- function(r, scenario, store_weights = FALSE) {
  set.seed(data_seed(scenario, r))
  X <- generate_data(scenario, n, p)
  
  res <- list(); prop <- list(); wts <- list()
  for (mech in mechanisms) {
    set.seed(missing_seed(scenario, r, mech, miss_rate))
    M <- make_missing(X, miss_vars, miss_rate, mech, slope = slope)
    X_mis <- X
    X_mis[M] <- NA
    
    # realized missing proportions
    prop[[mech]] <- data.frame(
      scenario = scenario, rep = r, mechanism = mech, rate = miss_rate,
      variable = miss_vars, realized = colMeans(M[, miss_vars]))
    
    for (wtd in c(FALSE, TRUE)) {
      fit <- hybrid_kNN(X_mis, k = k, weighted = wtd, seed = r,
                        num.trees = num_trees)
      rmse <- rmse_imputed(X, fit$X_hat, M, miss_vars)
      res[[length(res) + 1]] <- data.frame(
        scenario = scenario, rep = r, mechanism = mech, rate = miss_rate, k = k,
        method = if (wtd) "hybrid" else "standard",
        variable = c(miss_vars, "average"),
        rmse = c(rmse, mean(rmse)))
      
      if (store_weights && wtd) {
        for (v in names(fit$weights)) {
          wts[[length(wts) + 1]] <- data.frame(
            scenario = scenario, rep = r, mechanism = mech, rate = miss_rate,
            target = v, predictor = names(fit$weights[[v]]),
            weight = unname(fit$weights[[v]]))
        }
      }
    }
  }
  list(rmse = do.call(rbind, res), prop = do.call(rbind, prop),
       weights = if (length(wts)) do.call(rbind, wts) else NULL)
}


## ====================================================================
## 6. Main study
## ====================================================================
runs <- list()
for (sc in scenarios) {
  t0 <- Sys.time()
  cat("Scenario", sc, "...")
  runs <- c(runs, lapply(seq_len(R), run_one_rep, scenario = sc,
                         store_weights = (sc == "B")))
  cat(" done in", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
}
bind <- function(el) do.call(rbind, lapply(runs, `[[`, el))
results   <- bind("rmse")
realized  <- bind("prop")
weights_B <- bind("weights")

saveRDS(results,   file.path(out_dir, "rmse_results.rds"))
saveRDS(realized,  file.path(out_dir, "realized_missing.rds"))
saveRDS(weights_B, file.path(out_dir, "weights_scenarioB.rds"))


## ====================================================================
## 7. Summaries
## ====================================================================
mc_se <- function(x) sd(x) / sqrt(length(x))
key <- c("scenario", "mechanism", "rate", "variable")

# 7a. Realized missing proportions (target vs. realized)
prop_summary <- aggregate(realized ~ scenario + mechanism + rate + variable,
                          data = realized,
                          FUN = function(x) c(mean = mean(x), sd = sd(x),
                                              min = min(x), max = max(x)))
prop_summary <- do.call(data.frame, prop_summary)
write.csv(prop_summary, file.path(out_dir, "realized_missing_summary.csv"),
          row.names = FALSE)
cat("\nRealized missing proportions (mean over replications; target ",
    miss_rate, "):\n", sep = "")
print(format(prop_summary[order(prop_summary$scenario, prop_summary$variable), ],
             digits = 3), row.names = FALSE)

# 7b. Mean RMSE with Monte Carlo SE per method
rmse_summary <- aggregate(rmse ~ scenario + mechanism + rate + variable + method,
                          data = results,
                          FUN = function(x) c(mean = mean(x), se = mc_se(x)))
rmse_summary <- do.call(data.frame, rmse_summary)

# 7c. Paired differences hybrid - standard (negative = hybrid better)
wide <- reshape(results[, c(key, "rep", "method", "rmse")],
                idvar = c(key, "rep"), timevar = "method", direction = "wide")
wide$diff <- wide$rmse.hybrid - wide$rmse.standard
diff_summary <- aggregate(diff ~ scenario + mechanism + rate + variable,
                          data = wide,
                          FUN = function(x) c(mean = mean(x), se = mc_se(x),
                                              hybrid_better = mean(x < 0)))
diff_summary <- do.call(data.frame, diff_summary)

# one table for the report: combined RMSE over X1-X3
tab <- merge(
  reshape(rmse_summary[rmse_summary$variable == "average",
                       c("scenario", "mechanism", "rate", "method",
                         "rmse.mean", "rmse.se")],
          idvar = c("scenario", "mechanism", "rate"), timevar = "method",
          direction = "wide"),
  diff_summary[diff_summary$variable == "average", -4])
tab <- tab[order(tab$scenario, match(tab$mechanism, mechanisms)), ]
write.csv(rmse_summary, file.path(out_dir, "rmse_summary_by_variable.csv"),
          row.names = FALSE)
write.csv(diff_summary, file.path(out_dir, "rmse_difference_summary.csv"),
          row.names = FALSE)
write.csv(tab, file.path(out_dir, "rmse_table_average.csv"), row.names = FALSE)
cat("\nCombined RMSE over X1-X3 (mean and MC SE), paired difference hybrid - standard:\n")
print(format(tab, digits = 3), row.names = FALSE)

# 7d. RF weights in Scenario B: share of weight on the informative block
weights_B$informative <- weights_B$predictor %in% miss_vars
w_share <- aggregate(weight ~ mechanism + target + rep, FUN = sum,
                     data = weights_B[weights_B$informative, ])
w_share <- aggregate(weight ~ mechanism + target, data = w_share,
                     FUN = function(x) c(mean = mean(x), se = mc_se(x)))
w_share <- do.call(data.frame, w_share)
names(w_share)[3:4] <- c("share_X1toX3.mean", "share_X1toX3.se")
write.csv(w_share, file.path(out_dir, "weights_B_informative_share.csv"),
          row.names = FALSE)
cat("\nScenario B: total RF weight on the informative variables X1-X3",
    "(equal weights would give 2/9 = 0.222):\n")
print(format(w_share, digits = 3), row.names = FALSE)


## ====================================================================
## 8. Plots
## ====================================================================
scen_lab <- labeller(scenario = function(s) paste("Scenario", s))
wide$mechanism <- factor(wide$mechanism, levels = mechanisms)

# 8a. Paired RMSE differences, combined over X1-X3
p_diff <- ggplot(wide[wide$variable == "average", ],
                 aes(x = mechanism, y = diff)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "grey85", outlier.size = 0.6) +
  facet_wrap(~ scenario, scales = "free_y", labeller = scen_lab) +
  labs(x = NULL, y = "RMSE hybrid - RMSE standard\n(negative = hybrid better)",
       caption = paste0(100 * miss_rate, "% missing in X1-X3, k = ", k,
                        ", R = ", R)) +
  theme_bw()
ggsave(file.path(out_dir, "rmse_difference_boxplots.pdf"), p_diff,
       width = 7, height = 3.5)
ggsave(file.path(out_dir, "rmse_difference_boxplots.png"), p_diff,
       width = 7, height = 3.5, dpi = 200)

# 8b. Same, per variable
p_diff_var <- p_diff %+% wide[wide$variable != "average", ] +
  facet_grid(variable ~ scenario, scales = "free_y", labeller = scen_lab)
ggsave(file.path(out_dir, "rmse_difference_boxplots_by_variable.pdf"),
       p_diff_var, width = 7, height = 7)

# 8c. RF weights in Scenario B (mean weight per predictor and target)
w_mean <- aggregate(weight ~ mechanism + target + predictor,
                    data = weights_B, FUN = mean)
w_mean$mechanism <- factor(w_mean$mechanism, levels = mechanisms)
w_mean$predictor <- factor(w_mean$predictor, levels = paste0("X", 1:p))
p_w <- ggplot(w_mean, aes(x = predictor, y = weight, fill = mechanism)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 1 / (p - 1), linetype = "dashed") +
  facet_wrap(~ target, labeller = labeller(target = function(s) paste("Imputing", s))) +
  labs(x = "Predictor", y = "Mean RF weight",
       caption = "Dashed line: equal weight 1/(p-1)") +
  theme_bw() + theme(legend.position = "bottom",
                     axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(file.path(out_dir, "weights_scenarioB.pdf"), p_w, width = 9, height = 3.5)
ggsave(file.path(out_dir, "weights_scenarioB.png"), p_w, width = 9, height = 3.5,
       dpi = 200)

saveRDS(sessionInfo(), file.path(out_dir, "sessionInfo.rds"))
cat("\nAll results saved in '", out_dir, "/'.\n", sep = "")