#' Provides tools to calculate the theoretical hydrodynamic response
#' of an aquifer undergoing harmonic straining or pressurization. There are
#' two classes of models here: (1) for sealed wells, based on the model of
#' Kitagawa et al (2011), and (2) for open wells, based on the models of
#' Cooper et al (1965), Hsieh et al (1987), Rojstaczer (1988), and Liu et al
#' (1989). These models treat strain (or aquifer head) as an input to the
#' physical system, and fluid-pressure (or water height) as the output. The
#' applicable frequency band of these models is characteristic of seismic
#' waves, atmospheric pressure fluctuations, and solid earth tides.
#' 
#' @details
#' The following functions provide the primary features of the package:
#' 
#' \code{\link{well_response}} and \code{\link{open_well_response}}, which 
#' take in arguments for well- and aquifer-parameters, and the frequencies
#' at which to calculate the response functions.
#' They both access the constants-calculation routines as necessary, meaning
#' the user need not worry about those functions
#' (e.g., \code{\link{alpha_constants}}).
#' 
#' Helper functions: 
#' \describe{
#' \item{\code{\link{sensing_volume}}}{ can be used to compute the
#' sensing volume of fluid, for the specified well dimensions.}
#' }
#'
#' @section Scientific background:
#'
#' The underlying model is based upon the assumption that fluid flows radially
#' through an homogeneous, isotropic, confined aquifer.
#'
#' The underlying principle is as follows.  When a harmonic wave induces
#' strain in a confined aquifer (one having aquitards above and below it), 
#' fluid flows radially into, and out of a well penetrating the aquifer.
#' The flow-induced drawdown, \eqn{s}, is governed by the following 
#' partial differential equation, expressed in radial coordinates(\eqn{r}):
#' \deqn{
#' \frac{\partial^2 s}{\partial r^2} + \frac{1}{r} 
#' \frac{\partial s}{ \partial r} - \frac{S}{T}\frac{\partial s}{\partial t} = 0
#' }
#' where \eqn{S, T} are the aquifer storativity and transmissivity respectively.
#' 
#' The solution to this PDE, with periodic discharge boundary conditions,
#' gives the amplitude and phase response we wish to calculate.
#' The solution  for an open well was presented by
#' Cooper et al (1965), and subsequently modified by Liu et al (1989).
#' Kitagawa et al (2011) adapted the solution
#' of Hsieh et al (1987) for the case of a sealed well.
#' 
#' These models are applicable to any quasi-static process involving harmonic, 
#' volumetric strain of an aquifer 
#' (e.g. passing Rayleigh waves, or changes in the Earth's tidal potential). 
#' In practice, however, the presence of permeable fractures can violate the
#' assumption of isotropic permeability, which may substantially
#' alter the response by introducing shear-strain coupling. But these
#' complications are beyond the scope of this model.
#' 
#' @docType package
#' @name kitagawa-package
#' @aliases kitagawa
#' @title Spectral response of a water well to harmonic strains at seismic frequencies 
#' @author Andrew J. Barbour <andy.barbour@@gmail.com> 
#' 
#' @import kelvin
#' 
#' @references Abramowitz, M. and Stegun, I. A. (Eds.). "Kelvin Functions." 
#' \eqn{\S 9.9} in Handbook of Mathematical Functions with Formulas, Graphs, 
#' and Mathematical Tables, 9th printing. New York: Dover, pp. 379-381, 1972.
#' 
#' @references Cooper, H. H., Bredehoeft, J. D., Papadopulos, I. S., and Bennett, R. R. (1965),
#' The response of well-aquifer systems to seismic waves, 
#' \emph{J. Geophys. Res.}, \strong{70} (16)
#' 
#' @references Hsieh, P. A., J. D. Bredehoeft, and J. M. Farr (1987),
#' Determination of aquifer transmissivity from Earth tide analysis,
#' \emph{Water Resour. Res.}, \strong{23} (10)
#
#' @references Kitagawa, Y., S. Itaba, N. Matsumoto, and N. Koisumi (2011),
#' Frequency characteristics of the response of water pressure in a closed well to volumetric strain 
#' in the high-frequency domain,
#' \emph{J. Geophys. Res.}, \strong{116}, B08301
#'
#' @references Liu, L.-B., Roeloffs, E., and Zheng, X.-Y. (1989),
#' Seismically Induced Water Level Fluctuations in the Wali Well, Beijing, China,
#' \emph{J. Geophys. Res.}, \strong{94} (B7)
#'
#' @references Roeloffs, E. (1996),
#' Poroelastic techniques in the study of earthquake-related hydrologic phenomena,
#' \emph{Advances in Geophysics}, \strong{37}
#' 
#' @seealso \code{\link{well_response}}, 
#' \code{\link{open_well_response}}, 
#' \code{\link{sensing_volume}}, 
#' \code{\link{wrsp-methods}}
NULL
.kitEnvName = ".kitEnv"
.kitEnv = new.env()
.constants_in = ".kitConstants"
.kitConstants = list(
  radians=list(from.degrees=pi/180, to.degrees=180/pi),
  water=list( density=1000, bulkmod=2.2e9 ),
  gravity=9.80665
)

#
# Datasets, etc
#

#' Access to constants used by default
#' 
#' The response of an aquifer depends on its mechanical
#' and hydrological properties; if these are not known or 
#' specified, these constants are used.
#'
#' @details The function \code{\link{constants}}
#' shows the structure of (optionally),
#' and returns the assumed constants, which do \emph{not} reside in
#' the namespace.
#' 
#' \subsection{Values}{
# The constants here include:
#' \describe{
#' \item{For water: }{Density and bulk modulus}
#' \item{Gravity: }{Standard gravitational acceleration at 6371km radius (Earth)}
#' \item{Conversions: }{Degrees to radians}
#' }
#' }
#' @name kitagawa-constants
#' @seealso \code{\link{well_response}} and \code{\link{open_well_response}}
#' 
#' \code{\link{kitagawa-package}}
#' @family ConstantsCalculators
NULL

#' @rdname kitagawa-constants
#' @param do.str logical; should the structure be printed?
#' @name constants
#' @export
#' @return The constants, invisibly.
#' @examples
#' constants()
#' constants(FALSE) # returns invisibly
constants <- function(do.str=TRUE){
  const <- .kitConstants
  if (do.str) str(const, comp.str = "++++++++\n\t", no.list=TRUE, digits.d = 9)
  return(invisible(const))
}