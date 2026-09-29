# --------------------------------------------------------------------
# Author:
# *Enter your group number here, as well as names and student numbers*
# --------------------------------------------------------------------

## Simulation code comparing standard kNN (unweighted) and hybrid kNN
## (random forest importance weights). Submit this file together with
## hybrid_kNN.R and the report.

library(ggplot2)
library(MASS)
library(ranger)

source("hybrid_kNN.R")

set.seed(1)


# --------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------

n_runs <- 20                        # replications per configuration
n <- 200                             # number of observations
p <- 10                              # number of variables
k <- 5                               # number of neighbors
miss_prop <- c(0.1, 0.2)             # missing proportion per incomplete variable
scenarios <- c("unweighted_favoured", "weighted_favoured")
mechanisms <- c("MCAR", "MAR", "MNAR")
incomplete_cols <- c(1, 2)           # indices of variables receiving missing values


# --------------------------------------------------------------------
# Helper: Standard (Unweighted) kNN
# --------------------------------------------------------------------

# Implements standard unweighted kNN with uniform weights (w_m = 1)
# to serve as the baseline comparison.
standard_kNN <- function(X, k = 5, dist_metric = "euclidean", aggregate_method = "mean") {
  X_mat <- as.matrix(X)
  n_obs <- nrow(X_mat)
  p_vars <- ncol(X_mat)
  var_names <- colnames(X_mat)
  
  means <- colMeans(X_mat, na.rm = TRUE)
  sds <- apply(X_mat, 2, sd, na.rm = TRUE)
  sds[sds == 0 | is.na(sds)] <- 1
  X_scaled <- scale(X_mat, center = means, scale = sds)
  
  vars_with_na <- which(colSums(is.na(X_scaled)) > 0)
  
  for (target_idx in vars_with_na) {
    missing_rows <- which(is.na(X_scaled[, target_idx]))
    observed_rows <- which(!is.na(X_scaled[, target_idx]))
    pred_indices <- setdiff(seq_len(p_vars), target_idx)
    
    for (i in missing_rows) {
      row_i <- X_scaled[i, pred_indices, drop = FALSE]
      dists <- numeric(length(observed_rows))
      
      for (d in seq_along(observed_rows)) {
        j <- observed_rows[d]
        row_j <- X_scaled[j, pred_indices, drop = FALSE]
        mut_obs <- which(!is.na(row_i) & !is.na(row_j))
        
        if (length(mut_obs) == 0) {
          dists[d] <- Inf
        } else {
          diffs <- as.numeric(row_i[mut_obs] - row_j[mut_obs])
          # Uniform weights summing to 1 over mutually observed variables
          dists[d] <- sqrt(mean(diffs^2))
        }
      }
      
      valid_donors <- which(is.finite(dists))
      if (length(valid_donors) == 0) {
        val <- if (aggregate_method == "median") median(X_mat[, target_idx], na.rm = TRUE) else means[target_idx]
        X_scaled[i, target_idx] <- (val - means[target_idx]) / sds[target_idx]
        next
      }
      
      k_eff <- min(k, length(valid_donors))
      nearest_idx <- valid_donors[order(dists[valid_donors])[seq_len(k_eff)]]
      best_donors <- observed_rows[nearest_idx]
      donor_vals <- X_scaled[best_donors, target_idx]
      
      if (aggregate_method == "mean") {
        X_scaled[i, target_idx] <- mean(donor_vals)
      } else {
        X_scaled[i, target_idx] <- median(donor_vals)
      }
    }
  }
  
  X_hat <- t(apply(X_scaled, 1, function(r) r * sds + means))
  colnames(X_hat) <- var_names
  return(X_hat)
}


# --------------------------------------------------------------------
# Data generating processes
# --------------------------------------------------------------------

# Returns an n x p numeric matrix with column names.
# scenario "unweighted_favoured": all variables equally informative
# scenario "weighted_favoured": only some variables are informative
generate_data <- function(n, p, scenario) {
  if (scenario == "unweighted_favoured") {
    # All variables equally correlated (compound symmetry structure)
    rho <- 0.6
    sigma <- matrix(rho, nrow = p, ncol = p)
    diag(sigma) <- 1
    X <- mvrnorm(n = n, mu = rep(0, p), Sigma = sigma)
  } else if (scenario == "weighted_favoured") {
    # Variables 1 and 2 depend strongly on variable 3; remaining variables are pure noise
    X3 <- rnorm(n)
    X1 <- 1.8 * sin(X3) + 0.5 * rnorm(n)
    X2 <- 1.2 * (X3^2) + 0.5 * rnorm(n)
    
    p_noise <- p - 3
    noise_mat <- matrix(rnorm(n * p_noise), nrow = n, ncol = p_noise)
    X <- cbind(X1, X2, X3, noise_mat)
  } else {
    stop("Unknown scenario.")
  }
  
  colnames(X) <- paste0("V", seq_len(p))
  return(X)
}


# --------------------------------------------------------------------
# Missing data mechanisms
# --------------------------------------------------------------------

# Each function takes the complete matrix X, the column indices `cols` that
# receive missing values, and the target proportions `prop` (one per column).
# Each returns X with NA inserted. Implemented by hand, no mice::ampute().

make_mcar <- function(X, cols, prop) {
  X_na <- X
  n_obs <- nrow(X)
  for (idx in seq_along(cols)) {
    target_col <- cols[idx]
    p_target <- prop[idx]
    mask <- runif(n_obs) < p_target
    X_na[mask, target_col] <- NA
  }
  return(X_na)
}

make_mar <- function(X, cols, prop) {
  X_na <- X
  n_obs <- nrow(X)
  # Missingness in target columns depends on fully observed column V3 (col 3)
  z <- scale(X[, 3])
  
  for (idx in seq_along(cols)) {
    target_col <- cols[idx]
    p_target <- prop[idx]
    # Solve for intercept so that mean probability matches target proportion
    b0 <- uniroot(function(b) mean(plogis(b + 1.5 * z)) - p_target, interval = c(-10, 10))$root
    probs <- plogis(b0 + 1.5 * z)
    mask <- runif(n_obs) < probs
    X_na[mask, target_col] <- NA
  }
  return(X_na)
}

make_mnar <- function(X, cols, prop) {
  X_na <- X
  n_obs <- nrow(X)
  # Missingness depends on the values of the variable itself
  for (idx in seq_along(cols)) {
    target_col <- cols[idx]
    p_target <- prop[idx]
    z <- scale(X[, target_col])
    b0 <- uniroot(function(b) mean(plogis(b + 1.5 * z)) - p_target, interval = c(-10, 10))$root
    probs <- plogis(b0 + 1.5 * z)
    mask <- runif(n_obs) < probs
    X_na[mask, target_col] <- NA
  }
  return(X_na)
}


# --------------------------------------------------------------------
# Evaluation
# --------------------------------------------------------------------

# Imputation quality on the entries that were made missing.
# X_true ... complete data, X_hat ... imputed data, miss ... logical matrix
evaluate_imputation <- function(X_true, X_hat, miss) {
  # Normalized Root Mean Squared Error (NRMSE)
  true_vals <- X_true[miss]
  imp_vals <- X_hat[miss]
  
  rmse <- sqrt(mean((true_vals - imp_vals)^2))
  nrmse <- rmse / sd(true_vals)
  return(nrmse)
}


# --------------------------------------------------------------------
# Simulation loop
# --------------------------------------------------------------------

# Runs one replication for one configuration and returns a data frame row
# with the error of unweighted kNN and hybrid kNN.
run_once <- function(scenario, mechanism) {
  # 1. Generate full dataset
  X_clean <- generate_data(n = n, p = p, scenario = scenario)
  
  # 2. Inject missing values
  if (mechanism == "MCAR") {
    X_miss <- make_mcar(X_clean, cols = incomplete_cols, prop = miss_prop)
  } else if (mechanism == "MAR") {
    X_miss <- make_mar(X_clean, cols = incomplete_cols, prop = miss_prop)
  } else if (mechanism == "MNAR") {
    X_miss <- make_mnar(X_clean, cols = incomplete_cols, prop = miss_prop)
  }
  
  # Logical matrix of missing entries
  miss_mask <- is.na(X_miss)
  
  # 3. Standard unweighted kNN imputation
  X_hat_unweighted <- standard_kNN(X_miss, k = k, dist_metric = "euclidean", aggregate_method = "mean")
  err_unweighted <- evaluate_imputation(X_clean, X_hat_unweighted, miss_mask)
  
  # 4. Hybrid kNN imputation
  hyb_out <- hybrid_kNN(X_miss, k = k, dist_metric = "euclidean", aggregate_method = "mean")
  X_hat_hybrid <- hyb_out$X_hat
  err_hybrid <- evaluate_imputation(X_clean, X_hat_hybrid, miss_mask)
  
  return(data.frame(
    Scenario       = scenario,
    Mechanism      = mechanism,
    Err_Unweighted = err_unweighted,
    Err_Hybrid     = err_hybrid
  ))
}

# Build grid of configurations and run replications
grid <- expand.grid(
  run       = seq_len(n_runs),
  mechanism = mechanisms,
  scenario  = scenarios,
  stringsAsFactors = FALSE
)

cat("Starting simulation loop...\n")
results_list <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  if (i %% 25 == 0 || i == nrow(grid)) {
    cat(sprintf("Iteration %d of %d (%.1f%%)\n", i, nrow(grid), 100 * i / nrow(grid)))
  }
  results_list[[i]] <- run_once(grid[["scenario"]][i], grid[["mechanism"]][i])
}

results <- do.call(rbind, results_list)

if (!dir.exists("results")) dir.create("results")
saveRDS(results, "results/results.rds")


# --------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------

# Reshape to long format for ggplot2
results_long <- rbind(
  data.frame(
    Scenario  = results$Scenario,
    Mechanism = results$Mechanism,
    Method    = "Standard kNN",
    NRMSE     = results$Err_Unweighted
  ),
  data.frame(
    Scenario  = results$Scenario,
    Mechanism = results$Mechanism,
    Method    = "Hybrid kNN",
    NRMSE     = results$Err_Hybrid
  )
)

# Friendly scenario labels
results_long$Scenario <- factor(
  results_long$Scenario,
  levels = c("unweighted_favoured", "weighted_favoured"),
  labels = c("Scenario 1: Equal Relevance (p = 10)", "Scenario 2: Noise Present (p = 10)")
)

results_long[["Mechanism"]] <- factor(results_long[["Mechanism"]], levels = c("MCAR", "MAR", "MNAR"))
results_long[["Method"]] <- factor(results_long[["Method"]], levels = c("Standard kNN", "Hybrid kNN"))

# Generate comparison boxplot
p_sim <- ggplot(results_long, aes(x = Mechanism, y = NRMSE, fill = Method)) +
  geom_boxplot(outlier.size = 0.7, alpha = 0.85, width = 0.6) +
  facet_wrap(~ Scenario, scales = "free_y") +
  scale_fill_manual(values = c("Standard kNN" = "#56B4E9", "Hybrid kNN" = "#E69F00")) +
  theme_bw(base_size = 12) +
  labs(
    title = "Imputation Performance Comparison",
    subtitle = paste0("NRMSE across ", n_runs, " replications per setting (n = ", n, ", k = ", k, ")"),
    y = "NRMSE (Lower is Better)",
    x = "Missing Data Mechanism"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold")
  )

print(p_sim)
ggsave("results/simulation_plot.pdf", plot = p_sim, width = 8, height = 4.5)