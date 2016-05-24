#' Spectral response for a sealed well
#'
#' This is the primary function to calculate the response
#' for a sealed well.
#' 
#' @details
#' The response depends strongly on the physical properties
#' given. Default values are assumed where reasonable--for instance, 
#' the pore-fluid is assumed to be water--but considerable care 
#' should be invested in the choice of
#' parameters, unless the function is used in an optimization scheme.
#' 
#' Assumed values are:
#' \tabular{rlrl}{
#' \code{Avs}  \tab 1 \tab \tab amplification factor for volumetric strain\cr
#' \code{Aw}   \tab 1 \tab \tab amplification factor for water well\cr
#' }
#' 
#' The responses returned here are,
#' effectively, the amplification of water levels in a well, relative to 
#' the aquifer strain; or
#' \deqn{Z = \frac{z}{\epsilon} \equiv \frac{p}{\rho g \epsilon}}
#' If \code{as.pressure=TRUE}, then the responses are scaled by
#' \code{rho*grav} so that they represent water pressure relative to
#' aquifer strain:
#' \deqn{Z = \frac{p}{\epsilon}}
#' 
#' Not all parameters need to be given, but should be.  
#' For example, if
#' either \code{rho} or \code{grav} are not specified, they
#' are taken from \code{\link{constants}}.
#' \emph{Parameters which do not end in \code{.} do
#' not need to be specified (they may be excluded); if
#' they are missing, warnings will be thrown.}
#'
#' @name well_response
#' @export
#' 
#' @param omega  frequency, (see \code{freq.units})
#' @param T.     effective aquifer transmissivity \eqn{[m^2/s]}
#' @param S.     well storativity,  \eqn{[unitless]}
#' @param Vw.    well volume,	 \eqn{[m^3]}
#' @param Rs.    radius of screened portion,  \eqn{[m]}
#' @param Ku.    undrained bulk modulus,  \eqn{[Pa]}
#' @param B.     Skempton's coefficient,  \eqn{[unitless, bounded]}
#' @param Avs   amplification factor for volumetric strain \eqn{E_{kk,obs}/E_{kk}},  \eqn{[]}
#' @param Aw    amplification factor of well volume change for \eqn{E_{kk}},  \eqn{[]}
#' @param rho   fluid density \eqn{[kg/m^3]}
#' @param Kf    bulk modulus of fluid,  \eqn{[Pa]}
#' @param grav  local gravitational acceleration \eqn{[m/s^2]}
#' @param freq.units  set the units of \code{omega}
#' @param as.pressure logical; should the response for water pressure? (default is water heights)
#'
#' @return An object with class 'wrsp'
#' 
#' @author A. J. Barbour <andy.barbour@@gmail.com>
#'
#' @seealso 
#' \code{\link{wrsp-methods}} for a description of the class 'wrsp' and its methods,
#' \code{\link{sensing_volume}} to easily estimate the volume \code{Vw.}, and
#' \code{\link{kitagawa-package}} for references and more background.
#' @family WellResponseFunctions
#' 
#' @examples
#' #### dummy example
#' well_response(1:10, T.=1, S.=1, Vw.=1, Rs.=1, Ku.=1, B.=1)
#' 
#' #### a more physically realistic calculation:
#' # Physical params applicable for B084 borehole
#' # (see: http://pbo.unavco.org/station/overview/B084/ for details)
#' #
#' Rc <- 0.0508   # m, radius of water-sensing (2in)
#' Lc <- 146.9    # m, length of grouted region (482ft)
#' Rs <- 3*Rc     # m, radius of screened region (6in)
#' Ls <- 9.14     # m, length of screened region (30ft)
#' #
#' # calculate the sensing volume for the given well parameters
#' Volw <- sensing_volume(Rc, Lc, Rs, Ls) # m**3, ~= 1.8
#' #
#' Frqs <- 10**seq.int(from=-4,to=0,by=0.1) # log10-space
#' head(Rsp <- well_response(omega=Frqs, T.=1e-6, S.=1e-5, 
#' Vw.=Volw, Rs.=Rs, Ku.=40e9, B.=0.2, freq.units="Hz"))
#' 
#' # Access plot.wrsp:
#' plot(Rsp)
#'
well_response <- function(omega, T., S., Vw., Rs., Ku., B., 
                          Avs, Aw, rho, Kf, grav,
                          freq.units=c("rad_per_sec","Hz"),
                          as.pressure=TRUE) UseMethod("well_response")

#' @rdname well_response
#' @method well_response default
#' @S3method well_response default
well_response.default <- function(omega, T., S., Vw., Rs., Ku., B., 
                                  Avs, Aw, rho, Kf, grav,
                                  freq.units=c("rad_per_sec","Hz"),
                                  as.pressure=TRUE){
    # Enforce units of omega to be radians/sec
    freq.units <- match.arg(freq.units)
    fc <- switch(freq.units, rad_per_sec=1, Hz=2*pi)
    omega <- fc*omega
    #
    # Setup constants
    defA <- 1
    if (missing(Avs)){
      Avs <- defA
    }
    if (missing(Aw)){
      Aw <- defA
    }
    const <- kitagawa::constants(FALSE)
    if (missing(rho)){
      rho <- const$water$density
    }
    if (missing(Kf)){
      Kf <- const$water$bulkmod
    }
    if (missing(grav)){
      grav <- const$gravity
    }
    rhog <- rho*grav
    #
    # Alpha function
    Alpha. <- omega_constants(omega, c.type="alpha", S.=S., T.=T., Rs.=Rs.)
    #   A1, and A2 (functions of Phi and Psi, calculated internally)
    Amat <- alpha_constants(Alpha., c.type="A")
    stopifnot(ncol(Amat)==7)
    #  A1,2 are in Mod(A.[,6:7]) 
    A12 <- matrix(Mod(Amat[,6:7]),ncol=2)  # is complex, but imag is zero, so == abs
    stopifnot(ncol(A12)==2)
    A1 <- A12[,1]
    A2 <- A12[,2]
    #
    TVFRG <- 2 * pi * T. / omega / Vw. / rhog
    #
    #check Skemptons coefficient: it should be bound to [0,1]
    .in0to1(B.)
    #
    XX. <- Ku. * B. / Aw * TVFRG - A2
    rNum. <- XX. * XX. + A1 * A1
    #
    YY. <- Kf * TVFRG  -  A2
    rDen. <- YY. * YY. + A1 * A1
    ##
    ## complex response EQ 17
    cNum <- complex(real=(Ku. * B. / Aw * TVFRG - A2), imaginary=A1)
    cDen <- complex(real=(Kf * TVFRG  -  A2), imaginary=A1)
    #
    # z / epsilon == p / epsilon * rho * g
    # but is there a mistake?
    # / rhog
    wellresp <- -1 * Kf * Aw / Avs  *  cNum / cDen 
    ##
    # optionally scale to p/epsilon
    rhog <- ifelse(as.pressure, rhog, 1)
    Zunits <- ifelse(as.pressure, "P/E", "Z/E")
    wellresp <- wellresp * rhog
    ##
    omega <- omega/fc
    toret <- list(Aquifer=list(Transmiss=T., Storativ=S., Diffusiv=T./S.),
                  Well=list(Volume=Vw., ScreenRad=Rs.),
                  Solid=list(BulkModU=Ku.,Skemp=B.),
                  Fluid=list(BulkMod=Kf, Density=rho),
                  Amplification=list(Ekk=Avs, Well=Aw),
                  Omega=list(Units=freq.units),
                  Gravity=grav,
                  Model=list(Model="kitagwawa", Pressure=as.pressure),
                  Response=cbind(omega=omega, wellresp=wellresp),
                  Response.units=Zunits
    )
    class(toret) <- "wrsp"
    return(toret)
}
