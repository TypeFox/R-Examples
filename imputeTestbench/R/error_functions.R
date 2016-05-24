#' Root Mean Square Error Calculation
#'
#' takes difference between Original data and Predicted data as input
#' @param diff is the difference between Original data and Predicted data
#' @return rmseVal as Root Mean Square Error
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
  rmseVal <- sqrt(mean(diff^2))
  return(rmseVal)
}


#' Mean Absolute Error Calculation
#'
#' takes difference between Original data and Predicted data as input
#' @param diff is the difference between Original data and Predicted data
#' @return maeVal as Mean Absolute Error
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
  maeVal <- mean(abs(diff))
  return(maeVal)
}


#' Mean Absolute Percent Error Calculation
#'
#' takes difference between Original data and Predicted data as input
#' @param diff is the difference between Original data and Predicted data
#' @param dataIn as original input data
#' @return mapeVal as Mean Absolute Error
#' @export
#' @examples
#' ## Generate 100 random numbers within some limits
#' x <- sample(1:7, 100, replace = TRUE)
#' y <- sample(1:4, 100, replace = TRUE)
#' z <- mape((x - y),x)
#' z
# Function that returns Mean Absolute Error
mape <- function(diff, dataIn)
{
  mapeVal <- mean(abs(diff)* 100/dataIn)
  return(mapeVal)
}

#==============================================================================================
