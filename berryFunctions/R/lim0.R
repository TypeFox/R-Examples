#' axis limits with one end at zero
#' 
#' Calculates the range needed for ylim or xlim in plot, so that axis
#' starts at zero and is extended by 4\% at the other end
#' 
#' @return Vector with two values: 0 and by 4% contra-extended max (as for xaxs="r")
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 6.6.2013
#' @seealso The \code{\link{extendrange}()} utility in package \pkg{grDevices}
#' @references \code{methods(plot)}, \code{\link[graphics]{plot.default}}.
#'             Actually, I found extendrange via \code{plot.function} in \code{curve}
#' @keywords dplot
#' @export
#' @examples
#' 
#' # basic idea:
#' val <- c(3.2, 1.8, 4.5, 2.8, 0.1, 2.9) # just some numbers
#' plot(val, ylim=lim0(val) ) # you don't even have to set yaxs="i" ;-)
#' 
#' # "normal" plot:
#' plot(val)
#' par("usr")  # -0.076  4.676
#' 
#' # if y-axis is not allowed to go below 0, and we're too lazy to set yaxs="i":
#' plot(val, ylim=lim0(val) )
#' round( par("usr")  , digits=5) # 0.00000 4.66296
#' 
#' # with 0.04 extension as claimed by help page (1/27 in source code = 0.037):
#' plot(val, ylim=lim0(val, f=0.04) )
#' round( par("usr")  , digits=5) # zero is not included on axis anymore
#' 
#' b <- -val
#' plot(b)
#' plot(b, ylim=lim0(b) ) # works with only negative values as well
#' 
#' @param x Numeric. Vector with values
#' @param f Numeric. Extension factor. DEFAULT: 0.04 as in extendrange used eg. by \code{\link{curve}}
#' @param curtail Logical. Should the range returned be trimmed by 4\%? That way, 
#'         plotting doesn't need the default \code{\link{par}} xaxs or yaxs changed. DEFAULT: TRUE
#' 
lim0 <- function(
x,
f=1/27,
curtail=TRUE)
{
if(length(x)==1) x <- c(0,x)
r <- range(as.matrix(x), finite=TRUE)
r2 <- r + c(-f,f) * diff(r) # classical procedure of extendrange
r2[which.min(abs(r2))] <- 0 # set one end to zero
if(curtail) # if par xaxs is "r" as it is by default, first trim the range, so that
r2 + c(f,-f) * diff(r2) # in the plot command, we don't have to change yaxs or xaxs
else
r2
}
