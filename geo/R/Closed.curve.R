#' Closed smooth curve through positions
#' 
#' A smooth curve is fitted to the positions with periodic spline as
#' implemented in function \code{ps}.
#' 
#' 
#' @param gogn Data frame of positions to be smoothed
#' @param df Degrees of freedom
#' @param n Number of points returned per data point
#' @return Data frame of positions with components \item{lat, lon }{in decimal
#' degrees}
#' @note Missing function \code{spline.des} is required for \code{\link{ps}} to
#' work.
#' @seealso \code{\link{ps}}
#' @keywords manip
#' @export Closed.curve
Closed.curve <-
function(gogn, df = round(nrow(gogn)/2), n = 10)
{
  gogn$index <- c(1:nrow(gogn))
  tmpper <- c(1, nrow(gogn))
  x <- glm(lat ~ ps(index, df = df, period = tmpper), data = gogn)
  y <- glm(lon ~ ps(index, df = df, period = tmpper), data = gogn)
  r <- range(gogn$index)
  pred.frame <- data.frame(index = seq(r[1], r[2], length = nrow(gogn) * n))
  pred.frame$lat <- predict(x, pred.frame)
  pred.frame$lon <- predict(y, pred.frame)
  pred.frame[, c("lat", "lon")]
}

