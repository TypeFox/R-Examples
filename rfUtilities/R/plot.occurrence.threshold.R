#' @title Plot occurrence thresholds
#' @description Plot function for occurrence.threshold object 
#' @param  x        A occurrence.threshold object
#' @param  ...      Additional arguments passed to plot
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>                
#' @method plot occurrence.threshold
#' @export    	     
plot.occurrence.threshold <- function(x, ...) {
  if(x$statistic == "delta.ss") { 	
    mdl.lab = "abs-difference of sensitivity-specificity"
	  } 
	else if (x$statistic == "sum.ss") {
      mdl.lab = "sum of sensitivity-specificity"
      }
    else if (x$statistic == "kappa") {
      mdl.lab = "Kappa"  
  }	  
  dots <- as.list(match.call(expand.dots = TRUE)[-1])
  dots[["x"]] <- names(x$thresholds)
  dots[["y"]] <- x$thresholds
  dots[["type"]] <- "l"
    if (is.null(dots[["xlab"]]) & "xlab" %in% names(dots) == FALSE) dots[["xlab"]] <- "probability"  
    if (is.null(dots[["ylab"]]) & "ylab" %in% names(dots) == FALSE) dots[["ylab"]] <- mdl.lab
    if (is.null(dots[["main"]]) & "main" %in% names(dots) == FALSE) dots[["main"]] <- paste0(mdl.lab, " thresholds") 
  do.call("plot", dots)
} 
