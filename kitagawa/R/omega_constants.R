#' Calculate any constants that depend on angular frequency \eqn{\omega}
#' 
#' @description This function accesses the appropriate method to calculate the
#' \eqn{\omega}-dependent constant associated with the choice of \code{c.type}.
#' 
#' \emph{This function is not likely to be needed by the user.}
#' 
#' @details
#' \subsection{What is \code{"omega"}?}{
#' The name is in reference to
#'  radial frequency \eqn{\omega}, which is defined as
#' \deqn{\omega \equiv 2 \pi / \tau}
#' where \eqn{\tau} is the period of oscillation.
#' }
#' \subsection{What is the \code{"alpha"} calculation?}{
#'  
#'  The parameter \eqn{\alpha} is defined as
#'  \deqn{\alpha \equiv r_w \sqrt{\omega S / T}}
#'  where \eqn{r_w} is the radius of the well,
#'  where \eqn{S} is the storativity, and \eqn{T} is
#'  transmissivity.  See the parameter \code{...} for details
#'  on how to pass these values.
#'  
#'  This definition is common to many papers on the topic.  For example,
#'  it corresponds to  \strong{Equation 12} in Kitagawa et al (2011).
#' Because the computation of \eqn{\alpha} depends also on physical
#' properties, additional parameters can be
#' passed through (e.g. the transmissivity).  
#' }
#' \subsection{What is the \code{"diffusivity_time"} calculation?}{
#' This is a similar calculation to \code{\link{omega_norm}}. It uses
#' the effective hydraulic diffusivity \eqn{D}, which is defined in
#' this case as the ratio of transmissivity to storativity:
#' \deqn{D \equiv \frac{T}{S}}
#' }
#' 
#' @section Warnings Issued:
#' 
#' In the case \code{c.type='alpha'}, the 
#' parameters \code{S.}, \code{T.},  and \code{Rs.} should
#' be passed; otherwise, values are assumed to ensure the 
#' calculation does not fail, and the results are essentially meaningless.
#' 
#' Warnings will be issued if any necessary parameters are missing, indicating
#' default values were used.
#'
#' @name omega_constants
#' @export
#' 
#' @param omega   frequency,  \eqn{[rad/sec]}
#' @param c.type  the constant to calculate 
#' @param ...     additional params passed to calculator.  In the case of 
#'   \code{ctype="alpha"}, set 
#'   \code{S., T., Rs.}; and, in the case of
#'   \code{ctype="diffusivity_time"}, set
#'   \code{D.} or \code{S., T.}.
#'
#' @return Values of the constant repesented by \code{c.type} for the given
#' parameters
#' 
#' @author A. J. Barbour <andy.barbour@@gmail.com>
#'
#' @seealso \code{\link{alpha_constants}}, \code{\link{well_response}}, and
#' \code{\link{kitagawa-package}} for references and more background.
#' @family ConstantsCalculators
#'  
#' @examples
#' # alpha
#' omega_constants() # default is alpha, but will give warnings about S., T., Rs.
#' omega_constants(T.=1,S.=1,Rs.=1)  # 0, no warnings
#' omega_constants(1:10)  # sequence, with warnings about S., T., Rs.
#' omega_constants(1:10,T.=1,S.=1,Rs.=1) # sequence, no warnings
#' # diffusivity time
#' omega_constants(c.type="diffusivity_time", D.=1)  # 0, no warnings
#' omega_constants(c.type="diff", D.=1)  # 0, no warnings (arg matching)
#' omega_constants(c.type="diff")  # 0, warnings about S., T. because no D.
#' omega_constants(c.type="diff", S.=1)  # 0, warnings about T. because no D. or S.
omega_constants <-
function(omega=0, c.type=c("alpha","diffusivity_time"), ...) UseMethod("omega_constants")

#' @rdname omega_constants
#' @method omega_constants default
#' @S3method omega_constants default
omega_constants.default <-
  function(omega=0, c.type=c("alpha","diffusivity_time"), ...){
    #
    # switch constants-calculation method
    c.type <- match.arg(c.type)
    c.meth <- switch(c.type, 
                     alpha=".wc_alpha",
                     diffusivity_time=".wc_difftime")
    #
    cdef <- 1 # default constant, if S T or Rs are missing
    # here are the methods available:
    #
    # 1) Alpha
    #
    .wc_alpha.default <- function(omega, S., T., Rs.){
        # Kitagawa equation 12
        # storativity S, transmissivity T, radius of screened portion Rs
        if (missing(S.)){
          warning("storativity was missing, used defaults")
          S. <- cdef
        }
        if (missing(T.)){
          warning("tranmissivity was missing, used defaults")
          T. <- cdef
        }
        if (missing(Rs.)){
          Rs. <- cdef
          warning("radius was missing, used defaults")
        }
        alpha <- sqrt( omega * S. / T. ) * Rs.
        # ensure not null
        .nullchk(alpha)
        return(alpha)
    } # end .wc_alpha
    #
    # 2) Diffusivity time 
    #
    .wc_difftime.default <- function(omega, D., S., T.){
      # From Rojstaczer response
      if (missing(D.)){
        if (missing(S.) & missing(T.)){
          warning("diffusivity was missing, checked for S. and T.")
        }
        # diffusivity is just T/S, so check for those
        if (missing(S.)){
          warning("storativity was missing, used default")
          S. <- cdef
        }
        if (missing(T.)){
          warning("tranmissivity was missing, used default")
          T. <- cdef
        }
        D. <- T./S.
      }
      dtau <- sqrt( omega / 2 / D.)
      # ensure not null
      .nullchk(dtau)
      return(dtau)
    } # end .wc_difftime
    #
    # do the calculation with the method of choice
    c.calc <- function(...) UseMethod(c.meth)
    toret <- c.calc(omega,...)
    return(toret)
  }
