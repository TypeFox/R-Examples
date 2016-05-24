#' Gini coefficient of a distribution
#' 
#'  Takes in a distribution and returns gini coefficient
#'  @param y an integer or numeric distribution
#'  @details
#'  To compute the gini coefficient of a distribution, \code{gini} is the right
#'  function. It uses trapezoidal approximation to calculate the area of the curve.
#'  
#'  Lorenz curve is also plotted as output.
#'  @return gini coefficient of the distribution
#'  @author Akash Jain
#'  @seealso \code{\link{auc}}
#'  @examples
#'  # Distribution
#' dist <- c(1, 4, 7, 15, 10)
#'
#' # Gini coefficient
#' GINI <- gini(y = dist)
#'  @export
gini <- function(y) {
  if(class(y) != 'integer' && class(y) != 'numeric') {
    stop('Invalid input: y should be integer or numeric vector')
  } else {
    ginis <- function(data, i) {
      area <- (data$cumperx[i] - data$cumperx[i-1]) * (data$cumpery[i] + data$cumpery[i-1])
      return(area)
    }
    data <- rbind(data.frame(x = 0, y = 0), data.frame(x = 1, y = y))
    data <- data[order(data[, 'y']), ]
    data[, 'cumx'] <- cumsum(data[, 'x'])
    data[, 'cumy'] <- cumsum(data[, 'y'])
    data[, 'cumperx'] <- (data[, 'cumx']/max(data[, 'cumx'])) * 100
    data[, 'cumpery'] <- (data[, 'cumy']/max(data[, 'cumy'])) * 100
    vtGINI <- sapply(2:nrow(data), function(i) ginis(data, i))/10000
    gini <- 1 - sum(vtGINI, na.rm = TRUE)
    par(xaxs = 'i', yaxs = 'i')
    plot(x = data$cumperx, 
         y = data$cumpery, 
         type = 'l', 
         xlab = 'Cumulative proportion', 
         ylab = 'Cumulative share of y', 
         xlim = c(0, 100),
         ylim = c(0, 100),
         col = 'orange',
         lwd = 2)
    abline(a = 0, b = 1, col = 'gray60', lwd = 2)
    return(gini)
  }
}