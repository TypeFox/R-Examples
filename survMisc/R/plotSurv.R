#' @name plotSurv
#' @title Plot a \code{Surv} object.
#'
#' @param x A \code{Surv} object
#' @param l length of arrow. Length is \code{l / nrow(x)}
#' @param ... Additional arguments.
#'  \cr
#' These are passed to as \code{...} (an ellipsis) to the following functions, respectively:
#' \describe{
#'  \item{\code{graphics::arrows}}{for plotting right- or left-censored observations}
#'  \item{\code{graphics::segments}}{for plotting interval-censored observations}
#' }
#' @return A graph (base graphics).
#' The type of graph depends on the \code{type} of the \code{Surv} object.
#' This is given by \code{attr(s, which="type")} :
#'  \item{counting}{Lines with an arrow pointing right if right censored.}
#'  \item{right}{Lines with an arrow pointing right if right censored.}
#'  \item{left}{Lines with an arrow pointing left if left censored.}
#'  \item{interval}{If censored:
#'   \describe{
#'    \item{arrow points right}{right censored}
#'    \item{arrow points left}{left censored}
#'   }
#'   If not censored:
#'    \describe{
#'     \item{lines}{observations of more than one time point}
#'     \item{points}{observation of one time only (i.e. start and end times are the same)}
#'    }
#'  }
#'
#' @keywords plot
#' 
#' @seealso ?graphics::arrows
#' @seealso ?graphics::segments
#'
#' @rdname plotSurv
#' @method plot Surv
#' @aliases plot.Surv
#' @export
#'
#' @examples
#' df0 <- data.frame(t1=c(0, 2, 4, 6, NA, NA, 12, 14),
#'                   t2=c(NA, NA, 4, 6, 8, 10, 16, 18))
#' s5 <- Surv(df0$t1, df0$t2, type="interval2")
#' plot(s5)
#'
plot.Surv <- function(x, l=3, ...){
    stopifnot(inherits(x, "Surv"))
    type1 <- attr(x, which="type")
    l1 <- l / nrow(x)
    ## largest observed time (used for x axis on plot)
    max1 <- ifelse(type1=="interval" | type1=="counting",
                   max(x[,1:2][!is.na(x[, 1:2])]),
                   max(x[,1]))
    ## add jitter to ends of x axis
    j1 <- stats::runif(1, max1 / 50, 2 * max1 / 50)
    maxX <- max1 + j1
    minX <- 0 - j1
    ## blank plot with appropriate axes
    graphics::plot(x=0,
                   type="n",
                   xlim=c(minX, maxX),
                   ylim=c(0, nrow(x)),
                   xlab="time",
                   ylab="observation")
    graphics::grid()
    ## plot based on Surv type
    switch(type1,
           right = plotR(x, l1=l1, ...),
           left = plotL(x, l1=l1, ...),
           counting = plotC(x, l1=l1, ...),
           interval = plotI(x, l1=l1, max1=max1, ...))
    ## add title
    int1 <- ifelse(type1=="Interval",
                   "\n Point = exactly observed event time",
                   "")
    text1 <- " censored survival data\n Arrow = censored observation"
    tit1 <- paste(type1, text1, int1, sep="")
    graphics::title(tit1)
}
### helper functions lines
Seg <- function(x0, x1, y0, y1, ...){
    graphics::segments(x0=x0, x1=x1,
                       y0=y0, y1=y1,
                       ...)
}
## make arrows, default point to right as this is more common
## xpd allows arrows to be on the edge of plot margins
Arr <- function(x0, x1, y0, y1,
                l1, code=2, ...){
    graphics::arrows(x0=x0, x1=x1,
                     y0=y1, y1=y1,
                     length=l1,
                     angle=67.5,
                     code=code,
                     xpd=TRUE,
                     ...)
}
## for right censored data
plotR <- function(x, l1, ...){
    for (i in 1:nrow(x)) {
        if (x[i, "status"]==0) {
            Seg(x0=0, x1=x[i, 1], y0=i, y1=i, ...)
        } else {
            Arr(x0=0, x1=x[i, 1], y0=i, y1=i,
                l1=l1, ...)
        }}}
## for left censored data
plotL <- function(x, l1, ...) {
    for (i in 1:nrow(x)) {
        if (x[i,"status"]==0) {
            Arr(x0=0, x1=x[i, 1], y0=i, y1=i,
                l1=l1, code=1, ...)
        } else {
            Seg(x0=0, x1=x[i, 1], y0=i, y1=i, ...)
        }}}
## for counting format data
plotC <- function(x, l1, ...){
    for (i in 1:nrow(x)) {
        if (x[i,"status"]==0) {
            Seg(x0=x[i, 1], x1=x[i, 2], y0=i, y1=i, ...)
        } else {
            Arr(x0=x[i, 1], x1=x[i, 2], y0=i, y1=i,
                l1=l1, code=1, ...)
        }}}
## for interval format data
plotI <- function(x, l1, max1, ...) {
    for (i in 1:nrow(x)) {
        ## switch based on type of observation
        ## (use point if exact observation)
        switch(x[[i,3]] + 1,
               Arr(x0=x[i, 1], x1=max1, y0=i, y1=i,
                   l1=l1, code=1, ...),
               graphics::points(x[i, 1], i, pch=21, bg=1),
               Arr(x0=0, x1=x[i, 1], y0=i, y1=i,
                   l1=l1, code=1, ...),
               Seg(x0=x[i, 1], x1=x[i, 2], y0=i, y1=i, ...))
    }
}
