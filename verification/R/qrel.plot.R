qrel.plot <- function(A, ...){

    if(!is.element("quantile", class(A))) {

  	stop("qrel.plot: this function works only on quantile forecasts. \n")

    } else {

	qrelPlotDefault(y.i = A$y.i, obar.i = A$obar.i, prob.y = A$prob.y, ...)

    }
}
