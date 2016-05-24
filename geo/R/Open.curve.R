#' Open smooth curve through positions
#' 
#' A smooth curve is fitted to the positions with natural spline as implemented
#' in function \code{ns}.
#' 
#' 
#' @param gogn Data frame of positions to be smoothed
#' @param df Degrees of freedom
#' @param n Number of points returned per data point
#' @return Data frame of positions with components \item{lat, lon}{in decimal
#' degrees}
#' @note Function \code{ns} is missing, assume it is some sort of natural
#' spline.
#' @keywords manip
#' @export Open.curve
Open.curve <-
function(gogn, df = nrow(gogn)/2, n = 10)
{
	gogn$index <- c(1:nrow(gogn))
        df <- round(df)
	x <- glm(lat ~ ns(index, df = df), data = gogn)
	y <- glm(lon ~ ns(index, df = df), data = gogn)
	r <- range(gogn$index)
	pred.frame <- data.frame(index = seq(r[1], r[2], length = nrow(gogn) *
		n))
	pred.frame$lat <- predict(x, pred.frame)
	pred.frame$lon <- predict(y, pred.frame)
	pred.frame <- pred.frame[, c("lat", "lon")]
	return(pred.frame)
}

