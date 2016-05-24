#' Classification into groups
#' 
#' classify continuous values into categories with different methods:\cr
#' - linearly or logarithmically spaced equal intervals,\cr
#' - intervals based on quantiles (equally filled bins),\cr
#' - intervals based on  distance from the mean in normal distributions,\cr
#' - user specified class borders (e.g. for legal or critical limits).
#' 
#' @details 
#' Binning methods are explained very nicely in the link in the section References.\cr
#' \emph{nbins} indicates the number of classes (and thus, colors).\cr \cr
#'
#' \tabular{llll}{
#'  \bold{\code{method}}     \tab |  explanation                         \tab |  meaning of \code{breaks}                             \tab |  default   \cr
#'  ----------               \tab |  -----------                         \tab |  -----------                                          \tab |  -------   \cr
#'  \bold{equalinterval}     \tab |  \emph{nbins} equally spaced classes \tab |  nbins                                                \tab |  100       \cr
#'  \bold{quantile}          \tab |  classes have equal number of values \tab |  the quantiles (or number of them)                    \tab |  0:4/4     \cr
#'  \bold{standarddeviation} \tab |  normal distributions                \tab |  the number of sd in one direction from the mean      \tab |  3         \cr
#'  \bold{logspaced}         \tab |  \emph{nbins} logarithmically spaced \tab |  c(nbins, base),  see \code{\link{logSpaced}}         \tab |  c(100,1.2)\cr
#'  \bold{usergiven}         \tab |  custom breakpoints                  \tab |  personal breakpoint values (including ends of Range) \tab |  none      \cr
#' }
#' The default is set to equalinterval which makes sense for my original intent
#' of plotting lake depth (bathymetry measured at irregularly distributed points) on a linear color scale.\cr
#' This is the workhorse for \code{\link{colPoints}}.\cr
#' 
#' @return list with class numbers (index) and other elements for \code{\link{colPoints}}
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2014
#' @seealso \code{\link{colPoints}}
#' @references See this page on the effect of classification (binning) methods: \cr 
#' \url{http://uxblog.idvsolutions.com/2011/10/telling-truth.html}
#' @keywords classif
#' @export
#' @examples
#' 
#' classify( c(1:10, 20), breaks=12)
#' classify( c(1:10, 20), "q", breaks=0:10/10)
#' classify( c(1:10, 20), "s", sdlab=2 )
#' classify( c(1:10, 20), "s", sdlab=1, breaks=2 )
#' classify( c(1:10, 20), "u", breaks=c(5,27) )
#' classify( c(1:10, 20), "l")
#' 
#' @param x Vector with numeric values
#' @param method Character string (partial matching is performed). Classification method or 
#'        type of binning to compute the class breakpoints. See section Details. DEFAULT: "equalinterval")
#' @param breaks Specification for method, see Details. DEFAULT: different defaults for each method
#' @param Range Ends of color bar for method=equalinterval. DEFAULT: range(x, finite=TRUE)
#' @param sdlab Type of label and breakpoints if \code{method=standarddeviation}. 
#'        1 means \code{-0.5 sd, 0.5 sd}, 2 means \code{-1 sd, mean, 1 sd}, 
#'        3 means actual numbers for type 1, 4 means numbers for type 2.
#' @param quiet Suppress warnings, eg for values outside Range? DEFAULT: FALSE
#' 
classify <- function(
  x, 
  method="equalinterval", 
  breaks, 
  Range=range(x, finite=TRUE), 
  sdlab=1, 
  quiet=FALSE)
{
x <- as.numeric(x)
# error checking:
if(length(Range) != 2) stop("Range must have two values.")
if(diff(Range)==0)
   {
   if(!quiet) warning("The Range values were equal. Range is now extended.")
   Range[1] <- Range[1] -1
   Range[2] <- Range[2] +1
   }
# Partial matching of method:
PossibleValues <- c("equalinterval", "quantile", "logspaced", "standarddeviation", "usergiven")
method <- PossibleValues[pmatch(tolower(method),  PossibleValues)]
if(is.na(method)) stop("method can only be equalinterval, quantile, logspaced, standarddeviation, or usergiven (but the name can be abbreviated).")
# actual work:
if(method=="equalinterval") # --------------------------------------------------
{
if(missing(breaks)) breaks <- 100     # default for breaks           # input control
if(length(breaks)>1) stop("breaks must be a single value if method='equalinterval'.")
nb <- breaks                          # number of bins (classes)
bb <- seqR(Range, length.out=nb+1)    # bin borders
at <- pretty2(bb)                     # position of labels in colPointsLegend / -Hist
la <- at                              # labels -"-
ix <- cut(x, breaks=bb, labels=FALSE, include.lowest=TRUE) # index of class for each value of x
} else if(method=="quantile") # ------------------------------------------------
{
if(missing(breaks)) breaks <- 0:4/4
if(!missing(breaks) & length(breaks)==1) breaks <- seq(0,1, length.out=breaks)
if(any(breaks<0 | breaks>1)) stop("breaks must be between 0 and 1 if method='quantile'.")
nb <- length(breaks) - 1
bb <- unique(quantile(x, probs=breaks))
at <- bb
la <- signif(bb, 2) # rounding
ix <- cut(x, breaks=bb, labels=FALSE, include.lowest=TRUE)
} else if(method=="logspaced") # ------------------------------------------------
{
if(missing(breaks)) breaks <- c(100, 1.2)
if(length(breaks)!=2) stop("breaks must have two values if method='logspaced'.")
nb <- breaks[1]
bb <- logSpaced(base=breaks[2], n=breaks[1], min=0.999*Range[1], max=1.001*Range[2], plot=FALSE)
bb <- signif(bb, 5) # else min is in reality min + 1e-13
bb <- unique(bb) # as logspaced produces a lot of breaks in one region
at <- pretty2(bb)
la <- at
ix <- cut(x, breaks=bb, labels=FALSE, include.lowest=TRUE)
} else if(method=="standarddeviation") # ---------------------------------------
{
if(missing(breaks)) breaks <- 3
breaks <- as.integer(breaks)
if(length(breaks)>1) {breaks <- breaks[1]; if(!quiet) warning("breaks was vector. Only first element is used.")}
if(is.na(breaks)) stop("breaks must be an integer")
if(breaks <0) stop("breaks must be a positive integer >=1 if method='standarddeviation'.")
if(sdlab==2|sdlab==4)
  {
  nb <- 2*breaks
  bb <- mean(x) + (-breaks:breaks)*sd(x)
  }
else
  {
  nb <- 2*breaks - 1
  bb <- mean(x) + (-breaks:(breaks-1)+0.5)*sd(x)
  }
at <- bb
if(sdlab==2)           {la <- paste(-breaks:breaks, "sd"); la[breaks+1] <- "m"} else
if(sdlab==3 | sdlab==4) la <- signif(at, 2) else
                        la <- paste(-breaks:(breaks-1)+0.5, "sd")
ix <- cut(x, breaks=bb, labels=FALSE, include.lowest=TRUE)
} else if(method=="usergiven") # -----------------------------------------------
{
if(missing(breaks)) stop("breaks _must_ be specified if method is 'usergiven'.")
if(length(breaks)==1) {breaks <- c(min(x,na.rm=TRUE), breaks, max(x,na.rm=TRUE))
                       if(!quiet) warning("breaks were expanded by range (x).")}
nb <- length(breaks) - 1
bb <- breaks
at <- bb
la <- signif(breaks, 2)
ix <- cut(x, breaks=bb, labels=FALSE, include.lowest=TRUE)
} else  # ----------------------------------------------------------------------
stop("method went wrong internally. Please tell me! (berry-b@gmx.de).")
# Range Warning:
###
if(any(is.na(ix)))
  {
  ix[ x < min(bb, na.rm=TRUE) ] <- nb+1
  ix[ x > max(bb, na.rm=TRUE) ] <- nb+2
  }
if(min(bb,na.rm=TRUE) > min(x,na.rm=TRUE) | max(bb,na.rm=TRUE) < max(x,na.rm=TRUE) )
  if(!quiet) warning("There are values outside of the range of the given classes.\n",
      "These are given the index ", nb+1, " (lower) and ", nb+2, " (higher).")
# Results
list(nbins=nb, bb=bb, at=at, labels=la, index=ix)
} # Function end
