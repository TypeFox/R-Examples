#' Compute mean absolute percentage error
#' 
#'  Takes in actual and predicted linear response, and returns MAPE value
#'  @param y actual linear response
#'  @param yhat predicted linear response
#'  @details
#'  \code{mape} calculates the mean absolute percentage error in a predicted linear
#'  response.
#'  @return mean absolute percentage error
#'  @author Akash Jain
#'  @seealso \code{\link{actvspred}}, \code{\link{splitdata}}
#'  @examples
#'  # A 'data.frame' with y and yhat
#'  df <- data.frame(y = c(1.5, 2, 3.2),
#'                   yhat = c(3.4, 2.2, 2.7))
#'
#' # Compute mape
#' MAPE <- mape(y = df[, 'y'], yhat = df[, 'yhat'])
#'  @export
mape <- function(y, yhat) {
  if(class(y) != 'integer' && class(y) != 'numeric') {
    stop('Invalid input: y should be numeric or integer vector representing a linear response')
  } else if(class(yhat) != 'integer' && class(yhat) != 'numeric') {
    stop('Invalid input: yhat should be numeric or integer vector of predicted linear response')
  } else if(length(y) != length(yhat)) {
    stop('Invalid input: vectors y and yhat should have the same length')
  } else {
    mape <- mean(abs((y - yhat)/y))
    return(mape)
  }
}