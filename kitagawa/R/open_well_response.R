#' Spectral response for an open well
#'
#' This is the primary function to calculate the response
#' for an open (exposed to air) well.
#' 
#' @details
#' As opposed to \code{\link{well_response}}, this
#' calculates the theoretical, complex
#' well response for an unsealed (open) well.
#' 
#' The response depends strongly on the physical properties
#' given. Default values are assumed where reasonable--for instance, 
#' the pore-fluid is assumed to be water--but considerable care 
#' should be invested in the choice of
#' parameters, unless the function is used in an optimization scheme.
#' 
#' The responses returned here are,
#' effectively, the amplification of water levels in a well, relative to 
#' the pressure head in the aquifer; or
#' \deqn{Z = \frac{z}{h} \equiv \frac{\rho g z}{p}}
#' If \code{as.pressure=TRUE}, then the responses are scaled by
#' \code{rho*grav} so that they represent water levels relative to
#' aquifer pressure:
#' \deqn{Z = \frac{z}{p}}
#' 
#' Not all parameters need to be given, but should be.  
#' For example, if
#' either \code{rho} or \code{grav} are not specified, they
#' are taken from \code{\link{constants}}.
#' \emph{Parameters which do not end in \code{.} do
#' not need to be specified (they may be excluded); if
#' they are missing, warnings will be thrown.}
#' 
#' @section Models:
#' \subsection{\code{"rojstaczer"}}{
#' Rojstaczer (1988) is based on measurements of water level
#' and strain from volumetric or areal strainmeters.
#' }
#' \subsection{\code{"cooper"}, \code{"hsieh"}, and \code{"liu"}}{
#' Cooper et al (1965), Hsieh et al (1987) and Liu et al (1989) are based
#' on measurements of water level and 
#' displacements from seismometers; these 
#' models are expressed succinctly in Roeloffs (1996).
#' 
#' The sense of the phase shift 
#' for the Liu and Rojstaczer models are reversed from their original presentation, 
#' in order to account for differences in sign convention.
#' }
#'
#' @name open_well_response
#' @export
#' 
#' @param omega  numeric; frequency,  (see \code{freq.units})
#' @param T. numeric; effective aquifer transmissivity \eqn{[m^2/s]}
#' @param S. numeric; well storativity,  \eqn{[unitless]}
#' @param Rs. numeric; the \emph{radius} of the open (screened) section
#' @param rho numeric; fluid density (assumed if missing)
#' @param grav numeric; the local gravitational acceleration (assumed if missing)
#' @param z numeric; From Rojstaczer (1988): the depth from the water table (assumed if missing and if needed)
#' @param Hw numeric; height of water column above confined surface (assumed if missing and if needed)
#' @param Ta numeric; thickness of aquifer (assumed if missing and if needed)
#' @param freq.units character; setup the units of \code{omega}
#' @param model  character; use the response model from either
#'    Rojstaczer (1988),
#'    Liu et al (1989),
#'    Cooper et al (1965), or
#'    Hsieh et al (1987).
#' @param as.pressure logical; should the response be relative to aquifer pressure? (default is aquifer head)
#' 
#' @return An object with class 'owrsp'
#' 
#' @author A. J. Barbour <andy.barbour@@gmail.com>
#'
#' @seealso 
#' \code{\link{owrsp-methods}} for a description of the class 'owrsp' and its methods, and
#' \code{\link{kitagawa-package}} for references and more background.
#' @family WellResponseFunctions
#' 
#' @examples
#' OWR <- open_well_response(1:10,1,1)
#' plot(OWR)
#' OWR <- open_well_response(1/(1:200),1,1,Ta=100,Hw=10,model="liu",freq.units="Hz")
#' plot(OWR)
open_well_response <- function(omega, T., S.,
                               Rs.=(8/12)*(1200/3937),
                               rho, grav, z, Hw, Ta,
                               freq.units=c("rad_per_sec","Hz"),
                               model=c("rojstaczer","liu","cooper","hsieh"),
                               as.pressure=TRUE) UseMethod("open_well_response")
#' @rdname open_well_response
#' @method open_well_response default
#' @S3method open_well_response default
open_well_response.default <- function(omega, T., S., 
                                       Rs.=(8/12)*(1200/3937),
                                       rho, grav, z, Hw, Ta,
                                       freq.units=c("rad_per_sec","Hz"),
                                       model=c("rojstaczer","liu","cooper","hsieh"),
                                       as.pressure=TRUE){
  # Pick a model
  model <- match.arg(model)
  # Enforce units of omega to be radians/sec
  freq.units <- match.arg(freq.units)
  fc <- switch(freq.units, rad_per_sec=1, Hz=2*pi)
  omega <- fc*omega
  #
  # Setup constants
  const <- kitagawa::constants(FALSE)
  if (missing(rho)){
    rho <- const$water$density
  }
  if (missing(grav)){
    grav <- const$gravity
  }
  rhog <- rho*grav
  #
  # Diffusiv time in unified framework
  Dtau. <- omega_constants(omega, c.type="diffusivity_time", S.=S., T.=T.)
  # e.g., Rojstaczer 1988 Eq 11:
  #   Qp is
  #     z^2 omega / 2 / D
  #   and we want sqrt(Q')
  #   Dtau is 
  #     sqrt( omega / 2 / D.) == sqrt (omega * S / 2 / T)
  #   so
  #     sqrt(Qp) == z * sqrt(omega / 2 / D) == z * Dtau
  #
  if (model=="rojstaczer"){
    #
    Zunits <- ifelse(as.pressure, "P/E", "Z/E")
    #
    # Rojstaczer 1988
    # Eq A3 - A4
    #
    if (missing(z)){
      z <- 1
      warning("Depth from the water table 'z' not given. using default")
    }
    sQp <- z * Dtau.
    exptau <- exp(-sQp)
    #
    A. <- (exptau*cos(sQp) - 1)
    B. <- -1*exptau*sin(sQp)
    # Gain = [P / pg As epsilon] == [z / As epsilon]
    wellresp <- complex(real=A., imaginary=B.)
    # fix a -1 sign convention in Rojstaczer (relative to Kitagawa)
    amp <- Mod(wellresp)
    phs <- -1*Arg(wellresp)
    wellresp <- complex(modulus=amp, argument=phs)
    #
  } else if (model %in% c("liu","cooper","hsieh")){
    #    
    Zunits <- ifelse(as.pressure, "Z/P", "Z/H")
    # 
    # Calc various constants needed
    #
    Alpha. <- omega_constants(omega, c.type="alpha", S.=S., T.=T., Rs.=Rs.)
    KelPhiPsi <- alpha_constants(Alpha., c.type="Phi") #alpha, kel0, kel1, phi, psi
    Kel0. <- KelPhiPsi[,2] # zero order Kelvin function (complex)
    Phi. <- KelPhiPsi[,4] 
    Psi. <- KelPhiPsi[,5]
    kei. <- Im(Kel0.)
    ker. <- Re(Kel0.)
    omegsq <- omega**2 # used a few times
    #
    # Hw is the height of the water column above the upper limit of the aquifer
    if (missing(Hw)){
      Hw <- 1
      warning("water column height 'Hw' not given. using default")
    }
    # Ta is the aquifer thickness
    if (missing(Ta)){
      Ta <- 1
      warning("aquifer thickness 'Ta' not given. using default")
    }
    #
    if (model=="liu"){
      #
      # from Liu et al (1989), expressed in Roeloffs 1996 eq 22
      #
      # Eq A13
      U. <- (Ta/T.)*Kel0.
      # Eq A16 (Beta)
      onei <- complex(real=0,imaginary=1)
      gamma <- sqrt(2*onei*omega/(Rs.**2 * grav * U.))
      expgam <- exp(-1*gamma*Ta)
      exp2gam <- exp(-2*gamma*Ta)
      A. <- -1 * omegsq / grav * (Hw + (1 - expgam)/(1 + expgam)/gamma)
      B. <- -1 * onei * omega * U. * Rs.**2 * gamma * expgam / (1 - exp2gam)
      # Eq A20 -- x / h
      wellresp <- 1 / (A. + B. + 1)
      # fix a -1 sign convention in Liu (relative to Hsieh/Cooper)
      amp <- Mod(wellresp)
      phs <- -1*Arg(wellresp)
      wellresp <- complex(modulus=amp, argument=phs)
      #
    } else if (model=="hsieh"){
      #
      # from Hsieh et al (1987), Eq 12-16
      #
      # R = rc**2 * omega / 2 T
      R. <- Rs.**2 * (Dtau.**2 / S.)
      U1. <- Psi. * ker.  +  Phi. * kei.
      U2. <- Phi. * ker.  -  Psi. * kei.
      A. <- 1  -  R. * U1.
      B. <-       R. * U2.
      # Eq 12 -- x / h
      wellresp <- 1 / complex(real=A., imaginary=B.)
      #
    } else if (model=="cooper"){
      #
      # from Cooper et al (1965) eq 28, also 
      # expressed in Roeloffs 1996 eq 20 or Liu eq 3
      #
      # the effective height of the water column in the well
      He. <- Hw + 3*Ta/8 
      cT. <- omega * Rs.**2 / 2 / T.
      A. <- 1  -  cT. * kei.  -  omegsq * He. / grav
      B. <- cT. * ker.
      # the amplification of water level in the well relative to 
      # pressure head in the aquifer
      # Eq 28 -- A == x / h == rho * g x / p
      wellresp <- 1 / complex(real=A., imaginary=B.)
    }
  }
  #
  # optionally scale to Z/P from Z/H, or p/E from Z/E
  rhog <- ifelse(as.pressure, 1/rhog, 1)
  wellresp <- wellresp * rhog
  #
  omega <- omega/fc
  toret <- list(
    Aquifer=list(Transmiss=T., Storativ=S., Diffusiv=T./S., Thickness=ifelse(missing(Ta), 1, Ta)),
    Well=list(ScreenRad=Rs., Z=ifelse(missing(z), 1, z), Zw=ifelse(missing(Hw), 1, Hw)),
    Fluid=list(Density=rho),
    Omega=list(Units=freq.units),
    Gravity=grav,
    Model=list(Model=model, Pressure=as.pressure),
    Response=cbind(omega, wellresp),
    Response.units=Zunits
  )
  class(toret) <- "owrsp"
  return(toret)
}