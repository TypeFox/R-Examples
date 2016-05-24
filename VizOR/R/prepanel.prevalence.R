##' Adapts y-axis limits to prevalence range and any nonzero box width
##' 
##' When generating a \sQuote{Prevalence Plot} via \code{bwplot(panel=panel.prevalence)},
##' specifying \code{prepanel=prepanel.prevalence} spares one from calculating the
##' y-axis bounds.
##' 
##' @param x Vector of (numeric or difftime) treatment durations
##' @param y The treatment factor
##' @param \dots Used to pass \code{prevalence} and \code{box.width} parameters
##' @author David C. Norris
##' @seealso \code{\link{panel.prevalence}}
##' @keywords hplot
##' @export prepanel.prevalence 
prepanel.prevalence <- function (x, y, ...){
  prevalence <- list(...)$prevalence
  box.width <- list(...)$box.width
  if (is.null(box.width))
    box.width <- 0
  list(ylim=c(0,box.width+signif(max(prevalence), digits=1)))
}
