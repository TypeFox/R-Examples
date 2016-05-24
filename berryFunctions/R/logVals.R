#' Create log-axis values and labels
#' 
#' Create nice values and labels to write at logartihmic axes
#' 
#' @return A list with 
#'        \item{vals}{Values for lines and label positions}
#'        \item{labs}{Formatted values for labels} 
#'        \item{all}{Values for lines}
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Feb 2014
#' @seealso \code{\link{log10}}, \code{\link{logAxis}}, \url{http://r.789695.n4.nabble.com/expression-exponent-labeling-td4661174.html}
#' @keywords aplot dplot
#' @export
#' @examples
#' 
#' # Easiest use: vector with data (logVals automatically finds range):
#' y <- 10^runif(50, -1, 2)
#' plot(y, log="y") # not much control over placement and format of labels
#' plot(y, log="y", yaxt="n")
#' # now do this better, with custom bases:
#' lv <- logVals(y, base=c(1,2,5) )
#' axis(2, lv$vals, lv$labs, las=1)
#' 
#' # Default arguments:
#' lv <- logVals()
#' str(lv) # values, formatted labels, all 10^x values for lines
#' plot(1, ylim=c(1e-3, 1e4), log="y", yaxt="n", yaxs="i")
#' abline(h=lv$all, col=8 )
#' box("plot")
#' axis(2, lv$vals, lv$labs, las=1)
#' lines(seq(0.5, 1.5, len=50), 10^runif(50, -3, 4), col=2)
#' 
#' # Formatting labels:
#' logVals(                )$labs
#' logVals(scient=TRUE     )$labs
#' logVals(exponent=5      )$labs # expression with exponent, see logAxis
#' logVals(big.mark=" "    )$labs
#' logVals(big=".", dec=",")$labs # German style (not recommended)
#' 
#' @param from Lower exponent \emph{OR} vector with data
#' @param to High end
#' @param Range Or give from and to as range
#' @param base Bases to be used, eg. c(1,2,5)
#' @param big.mark Symbol separating thousands, eg. space, comma, dot, etc. see \code{\link{format}} and \code{\link{prettyNum}}
#' @param decimal.mark Character separating comma values, see \code{\link{format}} and \code{\link{prettyNum}}
#' @param scientific See \code{\link{format}}
#' @param exponent Starting at which exponent should \code{labs} be an expression with exponents?
#'        Compare to \code{\link{options}("scipen")}. This is mainly for \code{\link{logAxis}} and only for base 1. DEFAULT: Inf
#' @param expobase1 Should "n * " be appended before 10^exp if n=1? DEFAULT: FALSE
#' @param allbase Base for \code{$all} (for horizontal lines). DEFAULT: 1:9
#' @param \dots Ignored arguments
#' 
logVals <- function(
  from=-7,
  to=7,
  Range,
  base=1,
  big.mark="'",
  decimal.mark=".",
  scientific=FALSE,
  exponent=Inf,
  expobase1=FALSE,
  allbase=1:9,
  ...
  )
{
# Calculate the exponents from vector, if given as first argument:
if( missing(to)  &  NROW(from)>1  )
  {
  rng <- range(log10(from[from>0]), finite=TRUE)
  from <- floor(rng[1])
  to <- ceiling(rng[2])
  }
# or calculate the exponents from range, if given
if( !missing(Range)  )
  {
  from <- floor(Range[1])
  to <- ceiling(Range[2])
  }
# values for lines and labels:
vals <- base*10^rep(floor(from):ceiling(to), each=length(base))
# formatted values for labels:
labs <-  format(vals, big.mark=big.mark, trim=TRUE, scientific=scientific, 
                      drop0trailing=TRUE, decimal.mark=decimal.mark)
# change to expression if value > exponent :
change1 <- abs(log10(vals)) >= exponent  & log10(vals)%%1 ==0 # base=1
change2 <- abs(log10(vals)) >= exponent  & log10(vals)%%1 !=0 # other base
if(expobase1)
  {
  change2 <- change1 | change2        # all should be treated as other bases
  change1 <- rep(FALSE, length(vals))
  }
if(any(change1 | change2)) # only deal with expression if applicable:
  {
  w2c1 <- which(change1)  # w2c= which to change
  w2c2 <- which(change2)
  wnc <- which(!change1 & !change2)
  labs2 <- vector("expression", length(labs))
  for(i in w2c1) labs2[[i]] <- bquote(10^.(floor(log10(vals))[i]))
  for(i in w2c2) labs2[[i]] <- bquote(
         .(rep(base, length.out=length(vals))[i])%.%10^.(floor(log10(vals))[i]))
  for(i in wnc) labs2[[i]] <- labs[i]  # regular labels for all entries not
  labs <- labs2                        # affected by the exponent size rule
  }
# Values for lines:
all <- allbase * 10^rep(floor(from):ceiling(to), each=length(allbase))
# return end result
list(vals=vals, labs=labs, all=all)
}
