#' @importFrom knitr knit_print
#' @importFrom knitr asis_output
#' @title FlexTable custom printing function for knitr
#'
#' @description FlexTable custom printing function for knitr
#' 
#' @param x a \code{FlexTable} to be printed
#' @param ... further arguments, not used. 
#' @export
knit_print.FlexTable<- function(x, ...){
	asis_output(as.html(x))
}
