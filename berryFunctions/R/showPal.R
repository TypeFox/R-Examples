#' show color palettes
#'
#' Plot examples of the seqential and diverging color palettes in this package.
#' Do not use \code{rainbow}: \url{https://eagereyes.org/basics/rainbow-color-map}
#'
#' @return NULL
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Apr 2016
#' @seealso \code{\link{seqPal}}, \code{\link{divPal}}, package \code{RColorBrewer}
#' @keywords dplot color
#' @export
#' @examples
#' showPal()
#'
#' @param \dots Arguments passed to \code{\link{par}}
#'
showPal <- function(...)
{
op <- par(mfcol=c(8,2), mar=c(0,0,0,0), oma=c(0,0,1.8,0), yaxt="n", xaxt="n", ...)
on.exit(par(op), add=TRUE)
cex <- 3

# Sequential palette -----------------------------------------------------------
plot(rep(1, 12), pch=15, cex=cex, col=seqPal(12))            ; text(  6, 1, "default")
title(main="berryFunctions::seqPal", xpd=NA, line=0.5) 
plot(rep(1,300), pch=15, cex=cex, col=seqPal(300))           ; text(150, 1, "n=300")
plot(rep(1,  7), pch=15, cex=cex, col=seqPal(7))             ; text(  3, 1, "n=7 works, too")
plot(rep(1,300), pch=15, cex=cex, col=seqPal(300, extr=TRUE)); text(150, 1, "extr=TRUE")
plot(rep(1, 12), pch=15, cex=cex, col=seqPal(alpha=0.4))     ; text(  6, 1, "alpha=0.4 (semi-transparency)")
plot(rep(1, 12), pch=15, cex=cex, col=seqPal(reverse=TRUE))  ; text(  6, 1, "rev=TRUE")
plot(rep(1,300), pch=15, cex=cex, col=seqPal(300, yb=TRUE))  ; text(150, 1, "yb=TRUE")
plot(rep(1,300), pch=15, cex=cex, col=seqPal(300, yr=TRUE))  ; text(150, 1, "yr=TRUE")

# Diverging palette ------------------------------------------------------------
plot(rep(1, 12), pch=15, cex=cex, col=divPal(12))            ; text(  6, 1, "default")
title(main="berryFunctions::divPal", xpd=NA, line=0.5)
plot(rep(1,300), pch=15, cex=cex, col=divPal(300))           ; text(150, 1, "n=300")
plot(rep(1,  7), pch=15, cex=cex, col=divPal(7))             ; text(  3, 1, "n=7")
plot(rep(1,300), pch=15, cex=cex, col=seqPal(300,colors=c("darkblue","green","orange")))
text(150, 1, "col=c('darkblue','green','orange'))")
plot(rep(1, 12), pch=15, cex=cex, col=divPal(alpha=0.4))     ; text(  6, 1, "alpha=0.4 (semi-transparency)")
plot(rep(1, 12), pch=15, cex=cex, col=divPal(reverse=TRUE))  ; text(  6, 1, "rev=TRUE")
plot(rep(1,300), pch=15, cex=cex, col=divPal(300, rwb=TRUE)) ; text(150, 1, "rwb=TRUE")
plot(rep(1,300), pch=15, cex=cex, col=divPal(300, ryb=TRUE)) ; text(150, 1, "ryb=TRUE")

}
