# --------------------------------------------------------------------
# Author:
# *Enter your group number here, as well as names and student numbers*
# --------------------------------------------------------------------

library(ranger)

## Use this code skeleton to implement the hybrid kNN algorithm.  Please submit 
## your implementation together with your report and your simulation code.

## IMPORTANT: Please do not change the function name and make sure that you
##             use the correct input and output as specified below. This will 
##             simplify grading because each group's code will be similarly 
##             structured.  However, feel free to add other arguments to the 
##             function definition as needed.


# Input:
# X ... data matrix
# k ... number of neighbors
# add any other arguments you need

# Output
# A list with the following components:
# X_hat ....... the imputed data set
# weights ..... a list of variable importance weights from the random forest
#               models. That is, each list element contains the weights for 
#               the corresponding variable to be imputed. Use the variable name 
#               as the name of the list element.
# add any other output you want to return as additional components

hybrid_kNN <- function(X, k = 5, dist_metric = c("euclidean", "manhattan"), 
                       aggregate_method = c("mean", "median"), ...) {
  
  dist_metric <- match.arg(dist_metric)
  aggregate_method <- match.arg(aggregate_method)
  
  # Ensure matrix/data.frame format with column names
  X_mat <- as.matrix(X)
  n <- nrow(X_mat)
  p <- ncol(X_mat)
  
  if (is.null(colnames(X_mat))) {
    colnames(X_mat) <- paste0("V", seq_len(p))
  }
  var_names <- colnames(X_mat)
  
  # 1. Standardization (centering and scaling)
  means <- colMeans(X_mat, na.rm = TRUE)
  sds <- apply(X_mat, 2, sd, na.rm = TRUE)
  sds[sds == 0 | is.na(sds)] <- 1 # Guard against zero variance columns
  
  X_scaled <- scale(X_mat, center = means, scale = sds)
  
  # Find columns containing missing values
  vars_with_na <- which(colSums(is.na(X_scaled)) > 0)
  
  weights_list <- list()
  
  if (length(vars_with_na) == 0) {
    return(list(X_hat = X_mat, weights = weights_list))
  }
  
  # 2. Impute variable by variable
  for (target_idx in vars_with_na) {
    target_var <- var_names[target_idx]
    
    missing_rows <- which(is.na(X_scaled[, target_idx]))
    observed_rows <- which(!is.na(X_scaled[, target_idx]))
    
    # Predictor indices and names
    pred_indices <- setdiff(seq_len(p), target_idx)
    pred_names <- var_names[pred_indices]
    
    # --- Fit Random Forest using ranger to extract Variable Importance ---
    # Prepare training data using observed target rows
    train_df <- as.data.frame(X_scaled[observed_rows, , drop = FALSE])
    
    # Simple median imputation for predictors in the RF training set (if missing)
    for (col in pred_names) {
      if (anyNA(train_df[[col]])) {
        med_val <- median(train_df[[col]], na.rm = TRUE)
        train_df[[col]][is.na(train_df[[col]])] <- ifelse(is.na(med_val), 0, med_val)
      }
    }
    
    rf_formula <- as.formula(paste(target_var, "~ ."))
    rf_fit <- ranger::ranger(
      formula = rf_formula,
      data = train_df,
      importance = "permutation",
      num.trees = 100
    )
    
    # Extract importance and enforce non-negativity
    raw_weights <- ranger::importance(rf_fit)
    raw_weights <- pmax(raw_weights, 0)
    
    # If all weights are 0, default to uniform weights
    if (sum(raw_weights) == 0) {
      raw_weights <- setNames(rep(1, length(pred_names)), pred_names)
    }
    
    weights_vector <- raw_weights[pred_names]
    weights_list[[target_var]] <- weights_vector
    
    # --- Compute Weighted Distances & Impute ---
    for (i in missing_rows) {
      row_i <- X_scaled[i, pred_indices, drop = FALSE]
      donor_distances <- numeric(length(observed_rows))
      
      for (d in seq_along(observed_rows)) {
        j <- observed_rows[d]
        row_j <- X_scaled[j, pred_indices, drop = FALSE]
        
        # Mutually observed variables between row i and candidate donor j
        mut_obs <- which(!is.na(row_i) & !is.na(row_j))
        
        if (length(mut_obs) == 0) {
          donor_distances[d] <- Inf
          next
        }
        
        obs_names <- pred_names[mut_obs]
        w_sub <- weights_vector[obs_names]
        sum_w <- sum(w_sub)
        
        # Normalize weights to sum to 1 over mutually observed variables
        if (sum_w > 0) {
          w_norm <- w_sub / sum_w
        } else {
          w_norm <- rep(1 / length(mut_obs), length(mut_obs))
        }
        
        diffs <- as.numeric(row_i[mut_obs] - row_j[mut_obs])
        
        if (dist_metric == "euclidean") {
          donor_distances[d] <- sqrt(sum(w_norm * (diffs^2)))
        } else if (dist_metric == "manhattan") {
          donor_distances[d] <- sum(w_norm * abs(diffs))
        }
      }
      
      valid_donors <- which(is.finite(donor_distances))
      if (length(valid_donors) == 0) {
        # Fallback to column median/mean if no donors share observed variables
        val <- if (aggregate_method == "median") median(X_mat[, target_idx], na.rm = TRUE) else means[target_idx]
        X_scaled[i, target_idx] <- (val - means[target_idx]) / sds[target_idx]
        next
      }
      
      k_eff <- min(k, length(valid_donors))
      nearest_d_indices <- valid_donors[order(donor_distances[valid_donors])[seq_len(k_eff)]]
      nearest_donors_original <- observed_rows[nearest_d_indices]
      
      donor_values <- X_scaled[nearest_donors_original, target_idx]
      
      if (aggregate_method == "mean") {
        X_scaled[i, target_idx] <- mean(donor_values)
      } else if (aggregate_method == "median") {
        X_scaled[i, target_idx] <- median(donor_values)
      }
    }
  }
  
  # 3. Unstandardize back to original scale
  X_hat <- t(apply(X_scaled, 1, function(row) row * sds + means))
  colnames(X_hat) <- var_names
  
  return(list(
    X_hat = X_hat,
    weights = weights_list
  ))
}