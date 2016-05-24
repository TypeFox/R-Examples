#' Plot Bi-cross-validation(BCV) Errors
#'
#' Plot the average BCV entrywise MSE against the number of factors tried, 
#' with error bars and the best number of factors picked.
#'
#' The \code{esabcv} object contains the raw BCV result \code{result.list}, 
#' which is a matrix with dimension \code{c(nRepeat, (max.r + 1))} where \code{nRepeat}
#' is the number of BCV repeats and \code{max.r} is the maximum number of factors tried.
#' If either tail of the error curve dominates, then the user has the option to change the start
#' and end rank for plotting.
#'
#' 
#' @param x \code{esabcv} object, typically result of \code{\link{EsaBcv}}.
#' @param start.r the starting number of factors to display in the plot.
#' @param end.r the largest number of factors allowed to display in the plot. 
#' Default is NA, which means to make \code{end.r} as \code{max.r}.
#' @param xlab title for the x axis.
#' @param ylab title for the y axis.
#' @param main title for the plot.
#' @param col.line the line color.
#' @param ... other parameters to be passed through to plotting functions.
#' @return A plot ploting the average BCV entrywise MSE against the number of 
#' factors tried (start.r to \code{max.r + 1}), with error bars (one standard deviation)
#' in grey and selected number of factors marked by a red crossing.
#' @examples
#'
#' \dontrun{
#' data(simdat)
#' result <- EsaBcv(simdat$Y)
#' plot(result)
#' plot(result, start.r = 1)
#'
#' }
#' 
#' @export
plot.esabcv <-function(x, start.r = 0, end.r = NA, xlab = "Number of Factors",
					   ylab = "BCV MSE", main = "Bi-cross-validation Error", col.line = "BLUE", ...) {
  best.r <- x$best.r;
  x <- x$result.list;
  nRepeat <- nrow(x);
  max.r <- ncol(x) - 1;
  if (is.na(end.r) | end.r > max.r)
	  end.r <- max.r;
  means <- colMeans(x);  
  yrange <- range(means[start.r:end.r + 1]);
  plot(start.r:end.r + 1, means[start.r:end.r + 1], type="o", col = col.line, pch=21,  
	   xaxt = "n", xlab = xlab, ylab = ylab, ylim = yrange, ...);
  points(best.r + 1, means[best.r + 1], pch=4, col = "RED", cex = 2);
  axis(1, at = 0:max.r + 1, labels = 0:max.r);
  title(main)
}
