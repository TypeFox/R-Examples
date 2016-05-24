

#' @title cost.01
#' @description 0-1 loss function
#' @export cost.01
#' @return 0-1 loss
#' @param y Response
#' @param yhat The predicted value

cost.01 <- function(y, yhat)
{
    mean(y!=yhat)
}
