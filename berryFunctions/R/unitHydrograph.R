#' Unit Hydrograph
#' 
#' Calculate continuous unit hydrograph with given n and k (in the framework of the linear storage cascade)
#' 
#' @return Vector with the unit hydrograph along t
#' @note The sum under the UH should always be 1 (if t is long enough). This needs yet to be checked...
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, July 2013
#' @seealso \code{\link{lsc}} on how to estimate n and k for a given discharge dataset. \code{deconvolution.uh} in the package hydromad, \url{http://hydromad.catchment.org}
#' @keywords hplot ts
#' @export
#' @examples
#' 
#' Time <- 0:100
#' plot(Time, unitHydrograph(n=2,  k=3, t=Time), type="l", las=1,
#'      main="Unit Hydrograph - linear storage cascade")
#' lines(Time, unitHydrograph(n=2,  k=8, t=Time), col=2)
#' lines(Time, unitHydrograph(n=5.5,k=8, t=Time), col=4)
#' text(c(12, 20, 50), c(0.1, 0.04, 0.025), c("n=2, k=3","n=2, k=8","n=5.5, k=8"),
#'      col=c(1,2,4), adj=0)
#' 
#' # try several parameters (e.g. in Monte Carlo Simulation to estimate
#'   # sensitivity of model towards slight differences/uncertainty in parameters):
#' nreps <- 1e3 # 5e4 eg on faster computers
#' n <- rnorm(nreps, mean=2, sd=0.8); n <- n[n>0]
#' k <- rnorm(nreps, mean=8, sd=1.1); k <- k[k>0]
#' UH <- sapply(1:nreps, function(i) unitHydrograph(n=n[i], k=k[i], t=Time))
#' UHquant <- apply(UH, 1, quantile, probs=0:10/10, na.rm=TRUE)
#' if(interactive()) View(UHquant)
#' 
#' plot(Time, unitHydrograph(n=2, k=8, t=Time), type="l", ylim=c(0, 0.06), las=1)
#' # uncertainty intervals as semi-transparent bands:
#' for(i in 1:5)
#'    polygon(x=c(Time, rev(Time)), y=c(UHquant[i,], rev(UHquant[12-i,])),
#'            col=rgb(0,0,1, alpha=0.3), lty=0)
#' 
#' lines(Time, UHquant[6,], col=4)
#' lines(Time, unitHydrograph(n=2, k=8, t=Time))
#' 
#' # Label a few bands for clarity:
#' points(rep(24,3), UHquant[c(2,5,9),25], pch="+")
#' for(i in 1:3) text(25, UHquant[c(2,5,9)[i],25],
#'                         paste("Q", c(10,40,80)[i], sep=""), adj=-0.1, cex=0.7)
#' # And explain what they mean:
#' Explain <- "Q80: 80% of the 50000 simulations are smaller than this value"
#' legend("topright", bty="n", legend=Explain)
#' 
#' # Some n and k values are cut off at the left, that explains the shift from the
#' # median of simulations relative to the n2k8 line.
#' 
#' @param n Numeric. Number of storages in cascade.
#' @param k Numeric. Storage coefficient [1/s] (resistance to let water run out). High damping = slowly reacting landscape = high soil water absorbtion = high k.
#' @param t Numeric, possibly a vector. Time [s].
#' 
unitHydrograph <- function(
n,
k,
t)
{
if(length(n)>1 | length(k)>1) stop("n and k can only have one single value!
For vectorization, use sapply (see documentation examples).")
t^(n-1) / k^n / gamma(n) * exp(-t/k)  # some say /k^(n-1) for the second term!
}

