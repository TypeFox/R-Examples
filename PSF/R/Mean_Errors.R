#' Root Mean Square Error Calculation
#'
#' takes difference between Original data and Predicted data as input
#' @param diff is the difference between Original data and Predicted data
#' @return rmse_val as Root Mean Square Error
#' @export
#' @examples
#' ## Generate 100 random numbers within some limits
#' x <- sample(1:7, 100, replace = TRUE)
#' y <- sample(1:4, 100, replace = TRUE)
#' z <- rmse(x - y)
#' z
### RSME/MAE errors ============================================================================
# Function that returns Root Mean Squared Error
rmse <- function(diff)
{
  rmse_val <- sqrt(mean(diff^2))
  return(rmse_val)
}


#' Mean Absolute Error Calculation
#'
#' takes difference between Original data and Predicted data as input
#' @param diff is the difference between Original data and Predicted data
#' @return mae_val as Mean Absolute Error
#' @export
#' @examples
#' ## Generate 100 random numbers within some limits
#' x <- sample(1:7, 100, replace = TRUE)
#' y <- sample(1:4, 100, replace = TRUE)
#' z <- mae(x - y)
#' z
# Function that returns Mean Absolute Error
mae <- function(diff)
{
  mae_val <- mean(abs(diff))
  return(mae_val)
}
#==============================================================================================
