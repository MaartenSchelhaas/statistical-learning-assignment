# --------------------------------------------------------------------
# Author:
# *Enter your group number here, as well as names and student numbers*
# --------------------------------------------------------------------

install(ranodmForest)
library(randomForest)
library(ranger)

## Use this code skeleton to implement the hybrid kNN algorithm.  Please submit 
## your implementation together with your report and your simulation code.

## IMPORTANT: Please do not change the function name and make sure that you
##            use the correct input and output as specified below. This will 
##            simplify grading because each group's code will be similarly 
##            structured.  However, feel free to add other arguments to the 
##            function definition as needed.


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

hybrid_kNN <- function(X, k, ...) {
  # Loop through each column index
  for (j in seq_len(ncol(X))) {
    target_col_name <- colnames(X)[j]
    
    # Check if this column actually has missing values
    if (any(is.na(X[, j]))) {
      cat("Imputing column:", target_col_name, "\n")
      
      target_X <- X[, j]
      predictors_X <- X[, -j, drop = FALSE]
      
      # 1. Separate your target variable (X[, j]) and predictors (X[, -j])
      # 2. Train Random Forest on observed cases to get variable importances
      # 3. Compute weighted distances and perform kNN imputation
    }
  }
  
}
