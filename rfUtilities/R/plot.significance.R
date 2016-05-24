#' @title Plot random forests significance
#' @description Plot function for significance object 
#'
#' @param  x        A significance object
#' @param  ...      Additional arguments passed to plot
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'                   
#' @method plot significance
#'
#' @export    	     
plot.significance <- function(x, ...) {
  dots <- as.list(match.call(expand.dots = TRUE)[-1])
  dots[["type"]] <- "n"
  if( x$rf.type == "classification" ) {
  	if (is.null(dots[["xlab"]]) & "xlab" %in% names(dots) == FALSE) dots[["xlab"]] <-  "Error"
  	if (is.null(dots[["ylab"]]) & "ylab" %in% names(dots) == FALSE) dots[["ylab"]] <-  ""
  	if (is.null(dots[["main"]]) & "main" %in% names(dots) == FALSE) dots[["main"]] <- "Distribution of randomized OOB error"
    dots[["x"]] <- stats::density(x$RandOOB)
    dots[["x"]]$y <- dots[["x"]]$y/max(dots[["x"]]$y)
	dots[["xlim"]] <- c(min(c(x$RandOOB,x$test.OOB)), 1) 
      do.call("plot", dots)	
        graphics::polygon(dots[["x"]], col="blue")
        graphics::abline(v = x$test.OOB, col="black", lwd=1.5, lty=2)
  	    graphics::abline(v = stats::quantile(x$RandOOB, p = x$TestQuantile),lwd=1.5, lty=2, col="red") 
        graphics::legend("topright", c("model", "null"), bg="white", col=c("black","red"), 
	                     lty=c(2,2), lwd=c(1.5,1.5) )
							 
  } else if(x$rf.type == "regression" ) {
    if (is.null(dots[["xlab"]]) & "xlab" %in% names(dots) == FALSE) dots[["xlab"]] <-  "Error"
  	if (is.null(dots[["ylab"]]) & "ylab" %in% names(dots) == FALSE) dots[["ylab"]] <-  ""
  	if (is.null(dots[["main"]]) & "main" %in% names(dots) == FALSE) dots[["main"]] <- "Distribution of R-square in randomized models"
	dots[["x"]] <- stats::density(x$RandRsquare)
    dots[["x"]]$y <- dots[["x"]]$y/max(dots[["x"]]$y)
	dots[["xlim"]] <- c(min(x$RandRsquare), 1)
      do.call("plot", dots)  
        graphics::polygon(dots[["x"]], col="blue")
        graphics::abline(v = x$Rsquare, col="black", lwd = 1.5, lty = 2)
  		graphics::abline(v = stats::quantile(x$RandRsquare, p = x$TestQuantile),lwd = 1.5, lty = 2, col = "red") 
        graphics::legend("topright", c("model", "null"), bg="white",  
                         col=c("black","red"), lty=c(2,2), lwd=c(1.5,1.5) )				   
  }
}
