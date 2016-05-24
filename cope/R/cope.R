#' Coverage Probability Excursion (CoPE) Sets 
#'
#' The package cope computes and plots CoPE sets defined in Sommerfeld, Sain and
#' Schwartzman (2015)  for 2D functions. CoPE sets for
#' a real-valued target function \eqn{\mu(s)} on a two-dimensional domain 
#' are designed to bound the excursion set 
#' \eqn{{\mu(s) >= c}} from above and below with a predefined probability.
#' The target 
#' function can be a parameter in spatially
#' indexed linear regression \eqn{Y(s)=X*b(s)+ \epsilon(s)}, where \eqn{s}
#' is the spatial location, \eqn{X} is a known \eqn{n} by \eqn{p} design matrix,
#' \eqn{\epsilon(s)} is an error field and \eqn{Y(s)} is the observed data.
#' 
#' @section Major functions:
#' 
#'  \itemize{
#'  \item{\code{\link{ComputeCope}}}{ Computes CoPE sets for 2D data.}
#'  \item{\code{\link{PlotCope}}}{ Plots CoPE sets.}
#' }
#' 
#' @section Toy example functions:
#' 
#' \itemize{
#'  \item{\code{\link{ToySignal}}}{ Generates a toy signal.}
#'  \item{\code{\link{ToyNoise1}}, \code{\link{ToyNoise2}}, 
#'        \code{\link{ToyNoise3}}}{ Generates realizations of toy noise fields.}
#' }
#' 
#' @examples
#' # An example using the ToySignal and the Toy Noise1 of this package.
#' 
#' # Set sample size.
#' n = 30  
#' # Generate n realizations of the toy noise field.
#' Data = ToyNoise1(n = n)
#' # Add the toy signal to the noise.
#' Data$z = Data$z + rep(ToySignal()$z, n)
#' # Compute the CoPE sets.
#' CopeSet = ComputeCope(Data,level=4/3, mu=ToySignal()$z)
#' # Plot the result.
#' PlotCope(CopeSet)
#' 
#' @references M. Sommerfeld, S. Sain and A. Schwartzman. Confidence regions for 
#'             excursion sets in asymptotically Gaussian
#'             random fields, with an application to climate. Preprint, 2015. 
#' 
#' @import maps
#'
#' @docType package
#' @name cope
NULL