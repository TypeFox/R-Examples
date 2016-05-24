


#' @title cost.mse
#' @description mse
#' @export cost.mse
#' @return Mean Square Error
#' @param y Response
#' @param yhat The predicted value

cost.mse <- function(y, yhat)
{
    mean((y - yhat)^2)
}