# --------------------------------------------------------------------
# Author:
# *Enter your group number here, as well as names and student numbers*
# --------------------------------------------------------------------

## Simulation code comparing standard kNN (unweighted) and hybrid kNN
## (random forest importance weights). Submit this file together with
## hybrid_kNN.R and the report.

library(ggplot2)

source("hybrid_kNN.R")

set.seed(1)


# --------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------

n_runs <- 100                       # replications per configuration
n <- 200                            # number of observations
p <- 10                             # number of variables
k <- 5                              # number of neighbors
miss_prop <- c(0.1, 0.2)            # missing proportion per incomplete variable
scenarios <- c("unweighted_favoured", "weighted_favoured")
mechanisms <- c("MCAR", "MAR", "MNAR")


# --------------------------------------------------------------------
# Data generating processes
# --------------------------------------------------------------------

# Returns an n x p numeric matrix with column names.
# scenario "unweighted_favoured": all variables equally informative
# scenario "weighted_favoured": only some variables are informative
generate_data <- function(n, p, scenario) {
  # *enter your code here*
}


# --------------------------------------------------------------------
# Missing data mechanisms
# --------------------------------------------------------------------

# Each function takes the complete matrix X, the column indices `cols` that
# receive missing values, and the target proportions `prop` (one per column).
# Each returns X with NA inserted. Implemented by hand, no mice::ampute().

make_mcar <- function(X, cols, prop) {
  # *enter your code here*
}

make_mar <- function(X, cols, prop) {
  # *enter your code here*
}

make_mnar <- function(X, cols, prop) {
  # *enter your code here*
}


# --------------------------------------------------------------------
# Evaluation
# --------------------------------------------------------------------

# Imputation quality on the entries that were made missing.
# X_true ... complete data, X_hat ... imputed data, miss ... logical matrix
evaluate_imputation <- function(X_true, X_hat, miss) {
  # *enter your code here*
}


# --------------------------------------------------------------------
# Simulation loop
# --------------------------------------------------------------------

# Runs one replication for one configuration and returns a data frame row
# with the error of unweighted kNN and hybrid kNN.
run_once <- function(scenario, mechanism) {
  # *enter your code here*
}

# results <- do.call(rbind, lapply(...))
# saveRDS(results, "results/results.rds")


# --------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------

# *enter your code here*
