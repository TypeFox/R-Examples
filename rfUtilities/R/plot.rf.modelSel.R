#' @title Plot random forests model selection
#' @description Dot plot function for rf.modelSel importance vlaues 
#'
#' @param  x      A rf.modelSel object
#' @param  imp    Plot selected ("sel") or all ("all") importance used in model selection
#' @param  ...    Additional arguments passed to plot
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'                    
#' @method plot rf.modelSel
#'
#' @export    	     
plot.rf.modelSel <- function(x, imp = "sel",  ...) {
  plot.ms <- function(x, ...) {	  
    dots <- as.list(match.call(expand.dots = TRUE)[-1])
      x <- as.matrix( x )    
        ord <- rev(order(x[,1], decreasing=TRUE)[1:nrow(x)])
    dots[["x"]] <- x[ord,1]
    if (is.null(dots[["main"]]) & "main" %in% names(dots) == FALSE) dots[["main"]] <-  lable
    if (is.null(dots[["pch"]]) & "pch" %in% names(dots) == FALSE) dots[["pch"]] <-  20
  	do.call("dotchart", dots)
  }
    if( imp == "sel" ) { imp = x$sel.importance 
	  } else if (imp == "all") { imp = x$importance }
    if (x$s=="mir") {lable = "Row Standardization Variable Importance"} 	
    if (x$s=="se") {lable = "Standardized Error Variable Importance"}
  plot.ms(imp, ...) 
}
