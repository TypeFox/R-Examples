#' Comparison of actual and predicted linear response
#' 
#'  Takes in actual, predicted linear response and quantile value, and returns 
#'  average actual and predicted response in each quantile
#'  @param y actual linear response
#'  @param yhat predicted linear response
#'  @param n quantiles to be created for comparison
#'  @details
#'  \code{actvspred} function divides the data into \code{n} (given as input by the user)
#'  quantiles and computes mean of \code{y} and \code{yhat} for each quantile. All 
#'  \code{NA}'s in \code{n}, \code{y} and \code{yhat} are removed for calculation.
#'  
#'  The function also plots a line chart of average actual response and average predicted
#'  response over \code{n} quantiles. This plot can be used to visualize how close both
#'  the lines are.
#'  @return a data.frame with average actual and predicted response in each quantile
#'  @author Akash Jain
#'  @seealso \code{\link{mape}}, \code{\link{splitdata}}
#'  @examples
#'  # A 'data.frame' with y and yhat
#' df <- data.frame(y = c(1, 2, 3, 6, 8, 10, 15),
#'                  yhat = c(1.2, 2.5, 3.3, 6.9, 9.3, 6.5, 12.3))
#'
#' # Get actual vs predicted table
#' ACTVSPRED <- actvspred(y = df[, 'y'], yhat = df[, 'yhat'], n = 5)
#'  @export
actvspred <- function(y, yhat, n) {
  if(class(y) != 'integer' && class(y) != 'numeric') {
    stop('Invalid input: y should be numeric or integer vector representing a linear response')
  } else if(class(yhat) != 'integer' && class(yhat) != 'numeric') {
    stop('Invalid input: yhat should be numeric or integer vector of predicted linear response')
  } else if(length(y) != length(yhat)) {
    stop('Invalid input: vectors y and yhat should have the same length')
  } else if((class(n) != 'integer' && class(n) != 'numeric') | length(n) != 1) {
    stop('Invalid input: n should be an integer or numeric vector of length 1')
  } else {
    ntile <- as.integer(as.character(cut(yhat, 
                                         breaks = quantile(yhat, probs = seq(0, 1, by = 1/n), na.rm = TRUE), 
                                         labels = 1:n, 
                                         include.lowest = TRUE)))
    data <- data.frame(y = y, yhat = yhat, ntile = ntile)
    dataNomiss <- data[!is.na(data$ntile), ]
    summary <- aggregate(cbind(y, yhat) ~ ntile, FUN = mean, na.rm = TRUE)
    names(summary)[2:3] <- c('avgY', 'avgYhat')
    par(xaxs = 'i', yaxs = 'i')
    plot(x = c(0, summary$ntile), 
         y = c(0, summary$avgYhat), 
         type = 'l',
         ylim = c(0, max(max(summary$avgY), max(summary$avgY))),
         xlab = 'Decile', 
         ylab = 'Average',
         col = 'blue',
         lwd = 2)
    lines(x = c(0, summary$ntile), 
          y = c(0, summary$avgY), 
          type = 'l',
          col = 'red',
          lwd = 2)
    return(summary)
  }
}