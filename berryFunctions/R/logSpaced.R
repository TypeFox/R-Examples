#' Logarithmically spaced points
#' 
#' Calculates values that are in logarithmic distance from each other
#' e.g. to produce logarithmic interval borders
#' 
#' @return Vector or matrix, depending on base input
#' @note base >1 concentrates points at low values, base<1 at high values.
#'       \code{base} does not relate to \code{base} in \code{\link{log}}!
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Oct 2014
#' @seealso \code{\link{classify}}, \code{\link{log}}
#' @keywords arith
#' @export
#' @examples
#' 
#' logSpaced()
#' logSpaced(base=c(1.1, 1.5, 2), n=6, min=5, max=10)
#' d <- logSpaced(seq(0.8, 1.2, 0.025), main="logarithmically spaced points")
#' 
#' # the default base for the default n (20) will give an approximately equal
#' # bin width across the range on a logarithmic scale:
#' d <- logSpaced()
#' plot(d, rep(1,20), log="x")
#' 
#' @param base Base for calculations, can be a vector to compare several bases. DEFAULT: 1.1708
#' @param n Number of values to be calculated. DEFAULT: 30
#' @param min,max Range where n values are to be distributed, single values each. DEFAULT: 1,n
#' @param plot Should the points be plotted on a line? DEFAULT: TRUE
#' @param pch,las PointCharacter and Label Axis Style. DEFAULT: 3,1
#' @param ylab Y axis label. DEFAULT: "base"
#' @param \dots Further arguments passed to \code{\link{plot}}
#' 
logSpaced <- function(
base=1.1708,
n=20,
min=1, max=n,
plot=TRUE,
pch=3, las=1,
ylab="base",
...
)
{
# input control:
min <- min[1] # in case someone tries to vectorize these
max <- max[1]
n <- n[1]
if(!is.numeric(base)) stop("base must be numeric.")
if(!is.vector(base)) stop("base must be a vector.")
base <- sort(unique(base))
# calculate values (x)
x <- base[1]^(1:n-1)
# map to min-max:
if(all(x==1)) x <- 1:n # if base=1
lmc <- coef(lm(c(min,max)~x[c(1,n)]))
x <- lmc[2]*x + lmc[1]
# plot first value:
if(plot) plot(x, rep(base[1],n), ylim=range(base), pch=pch, las=las, ylab=ylab, ...)
# for all other base values, if this is a vector:
if(length(base)==1) return(x) else
  {
  x <- sapply(base, function(y) y^(1:n-1))
  colnames(x) <- base
  for(i in 1:ncol(x))
    {
    if(all(x[,i]==1)) x[,i] <- 1:n
    lmc <- coef(lm(c(min,max)~x[c(1,n),i]))
    x[,i] <- lmc[2]*x[,i] + lmc[1]
    }
  if(plot) for(i in 2:ncol(x)) points(x[,i], rep(base[i],n), pch=pch)
  return(x)
  }
}

