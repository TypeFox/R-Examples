#' Print Function For Johnson-Neyman
#'
#' Prints function for objects of class \code{"jn"}
#'
#' @param x An object of class \code{"jn"}.
#' @param \dots Additional arguments (not supported yet).
#' @return none
#'
#' @examples
#' \dontrun{
#' myModel <- lm('DV ~ IV + MOD', data=someData)
#' jnresults <- jn(myModel, dv='DV', iv='IV', mod='MOD')
#' jnresults
#' }
#' @rdname print.jn
#' @export

print.jn <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nConditional effects of ", x$iv, " on ", x$dv, " at values of ", x$mod ,"\n")
  showrect <- cbind(round(x$printsummary$x, digits=4), round(x$printsummary$y, digits=4), round(x$printsummary$se, digits=4), round(x$printsummary$t, digits=4), round(x$printsummary$p, digits=4), round(x$printsummary$llci, digits=4), round(x$printsummary$ulci, digits=4))
  colnames(showrect) <- c(x$mod,'Effect','se','t','p','llci','ulci')
  rownames(showrect) <- rep('',nrow(showrect))
  print(showrect)
}
