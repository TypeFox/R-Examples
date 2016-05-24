#' Difference of estimated K functions
#' 
#' \code{kdest} determines the difference in estimated K functions for a set of cases and controls.
#' 
#' This function relies internally on the \code{Kest} and \code{eval.fv} functions from the \code{spatstat} package.  So the arguments are essentially the same as the \code{Kest} function.  See the documentation of the \code{Kdest} for more details about the various arguments.
#' 
#' @param x A \code{ppp} object from the \code{spatstat} package with marks for the case and control groups.
#' @param case The position of the name of the "case" group in levels(x$marks).  The default is 2.
#' @param nsim An non-negative integer.  Default is 0.  The difference in estimated K functions will be calculated for \code{nsim} data sets generated under the random labeling hypothesis.
#' @param level Confidence level of confidence envelopes.  Ignoried if \code{nsim} is 0.
#' @param r Optional. Vector of values for the argument r at which K(r) should be evaluated. Users are advised not to specify this argument; there is a sensible default.
#' @param breaks This argument is for internal use only.
#' @param correction Optional. A character vector containing any selection of the options "none", "border", "bord.modif", "isotropic", "Ripley", "translate", "translation", "none", "good" or "best". It specifies the edge correction(s) to be applied.
#' @param nlarge Optional. Efficiency threshold. If the number of points exceeds nlarge, then only the border correction will be computed (by default), using a fast algorithm.
#' @param domain Optional. Calculations will be restricted to this subset of the window. See Details.
#' @param var.approx	Logical. If TRUE, the approximate variance of Kest(r) under CSR will also be computed.
#' @param ratio	Logical. If TRUE, the numerator and denominator of each edge-corrected estimate will also be saved, for use in analysing replicated point patterns.
#'
#' @return Returns an \code{fv} object.  See documentation for \code{spatstat::Kest}.
#' @author Joshua French
#' @import spatstat
#' @importFrom stats quantile
#' @export
#' @seealso \code{\link[spatstat]{Kest}}
#' @references Waller, L.A. and Gotway, C.A. (2005).  Applied Spatial Statistics for Public Health Data.  Hoboken, NJ: Wiley.  Kulldorff, M. (1997) A spatial scan statistic. Communications in Statistics -- Theory and Methods 26, 1481-1496.
#' @examples 
#' data(grave)
#' kd1 = kdest(grave)
#' plot(kd1, iso ~ r, ylab = "difference", legend = FALSE, main = "")
#' kd2 = kdest(grave, nsim = 9, level = 0.8)
#' plot(kd2)

kdest = function(x, case = 2, nsim = 0, level = 0.95, r=NULL, breaks=NULL, correction=c("border", "isotropic", "Ripley", "translate"), nlarge=3000, domain=NULL, var.approx=FALSE, ratio=FALSE)
{
  if(!is.element("ppp", class(x))) stop("x must be a ppp object")
  if(is.null(x$marks)) stop("x must be marked as cases or controls")
  if(!is.factor(x$marks)) stop("The marks(x) must be a factor")
  nlev = length(levels(x$marks))
  if(case < 1 || case > nlev) stop("case must be an integer between 1 and length(levels(x$marks))")
  if(nsim < 0 | !is.finite(nsim)) stop("nsim must be a non-negative integer")
  if(length(level) != 1) stop("level must have length 1")
  if(level <= 0 | level >= 1) stop("level should be between 0 and 1")
  
  if(nsim == 0)
  {
    out = kd(x, case = case, 
             r = r, breaks = breaks, correction = correction, 
             nlarge = nlarge, domain = domain, 
             var.approx = var.approx, ratio = ratio)
    out = list(out)
  }
  else
  {
    #min/max envelope
    out = spatstat::envelope(x, kd, case = case, nsim = nsim, savefuns = TRUE, 
                   simulate = expression(rlabel(x, permute = TRUE)), 
                   r = r, breaks = breaks, correction = correction, 
                   nlarge = nlarge, domain = domain, 
                   var.approx = var.approx, ratio = ratio)
    #confidence band envelope
    simfuns <- as.data.frame(attr(out, "simfuns"))
    simfuns[,1] <- out$obs
    l = apply(simfuns, 1, stats::quantile, prob  = (1 - level)/2)
    u = apply(simfuns, 1, stats::quantile, prob = 1 - (1-level)/2)
    out = list(out, r = out$r, qlo = l, qhi = u)
  }
  class(out) = "kdenv"
  return(out)
}