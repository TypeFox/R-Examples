#' Nash-Sutcliffe efficiency
#' 
#' Compute Nash-Sutcliffe efficiency from two vectors. Based on eval.NSeff in RHydro Package
#' 
#' Basically does the following, after some input checks: 1 - ( sum((obs - sim)^2) / sum((obs - mean(obs))^2) )
#' 
#' @return Numeric.
#' @note NAs are omitted with warning.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, July 2013
#' @seealso \code{\link{rmse}}, \code{\link{rsquare}}
#' @keywords ts
#' @export
#' @examples
#' 
#' SIM <- dbeta(1:40/40, 3, 10)
#' OBS <- SIM + rnorm(20,,0.2)
#' plot(OBS)
#' lines(SIM)
#' nse(OBS, SIM)
#' 
#' @param obs Vector with observational data.
#' @param sim simulated data.
#' 
nse <- function(
obs,
sim)
{
if(!(is.vector(obs) & is.vector(sim))) stop("Input is not a vector.")
if(length(obs) != length(sim)) stop("Vectors are not of equal length.")
if(any(is.na(obs)|is.na(sim)))
     {
     Na <- which(is.na(obs)|is.na(sim))
     warning(length(Na), " NAs were omitted from ", length(obs), " data points.")
     obs <- obs[-Na] ; sim <- sim[-Na]
     } # end if NA
1 - ( sum((obs - sim)^2) / sum((obs - mean(obs))^2) )
}
