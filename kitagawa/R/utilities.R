#' @title General utility functions
#' 
#' @description
#' General utility functions
#' 
#' @details
#' \code{\link{.nullchk}} quickly checks for \code{NULL} and \code{NA},
#' and raises an error if \code{TRUE}; 
#' \emph{This function is not likely to be needed by the user.}
#' 
#' \code{\link{.in0to1}} checks if values are numeric and in [0,1] (inclusive).
#' 
#' \code{\link{is.wrsp}} and \code{\link{is.owrsp}} report whether an object 
#' has S3 class 'wrsp' or 'owrsp', respectively.  Such an object
#' would be returned by, for example, \code{\link{well_response}}.
#' 
#' @name kitagawa-utilities
#' @docType methods
#'
#' @seealso \code{\link{kitagawa-package}}
#' 
#' @author A. J. Barbour <andy.barbour@@gmail.com>
#' 
#' @param X  something to be checked (vector, scalar, wrsp object, ...)
#' 
#' @examples
#' \dontrun{
#' .nullchk(1:10) # OK
#' .nullchk(NULL) # error
#' .nullchk(c(1:10,NULL)) # error
#' .nullchk(NA) # error
#' .nullchk(c(1:10,NA)) # error
#' 
#' .in0to1(1:10) # error
#' .in0to1(NULL) # error
#' .in0to1(c(1:10,NULL)) # error
#' .in0to1(NA) # error
#' .in0to1(c(1:10,NA)) # error
#' .in0to1(c(1,NA)) # error
#' 
#' is.wrsp(1) # FALSE
#' }
NULL

#' @name .nullchk
#' @rdname kitagawa-utilities
#' @export
.nullchk <- function(X){
  stopifnot(!is.null(X) & !(NA %in% X))
}

#' @name .in0to1
#' @rdname kitagawa-utilities
#' @export
.in0to1 <- function(X){
  X <- as.numeric(X)
  stopifnot((X >= 0) & (X <= 1))
}

#' @name is.wrsp
#' @rdname kitagawa-utilities
#' @export
is.wrsp <- function(X) inherits(X, "wrsp")

#' @name is.owrsp
#' @rdname kitagawa-utilities
#' @export
is.owrsp <- function(X) inherits(X, "owrsp")

#' Dimensionless frequency from diffusivity and depth
#' @details
#' Dimensionless frequency \eqn{Q} is defined as \deqn{Q=\frac{z^2 \omega}{2 D}}
#' where 
#' \eqn{z} is the well depth,
#' \eqn{\omega} is the angular frequency
#' and \eqn{D} is the hydraulic diffusivity
#' 
#' @param omega numeric; angular frequency
#' @param Diffusiv numeric; hydraulic diffusivity
#' @param z numeric; depth
#' @param invert logical; should \code{omega} be taken as normalized
#' frequency?
#' @return \code{\link{omega_norm}} returns dimensionless frequency, unless \code{invert=TRUE}
#' where it will assume \code{omega} is dimensionless frequency, and return radial frequency.
#' 
#' @author A. J. Barbour <andy.barbour@@gmail.com>
#' @family utilities
#' @seealso \code{\link{open_well_response}}, \code{\link{kitagawa-package}}
#' @export
omega_norm <- function(omega, Diffusiv, z, invert=FALSE){
  # From Rojstaczer 1988 Eq 11
  # omega <- Q * 2 * Diffus / z**2
  # and thus
  # Q <- omega * z**2 / 2 / Diffus
  stopifnot(Diffusiv>0)
  #
  z2 <- z*z
  if (invert){
    Q <- omega
    frq <- Q * 2 * Diffusiv / z2
  } else {
    # Diffusiv time in unified framework
    Dtau. <- omega_constants(omega, c.type="diffusivity_time", D.=Diffusiv)
    #   Qp is
    #     z^2 omega / 2 / D
    #   and we want sqrt(Q')
    #   Dtau is 
    #     sqrt( omega / 2 / D.) == sqrt (omega * S / 2 / T)
    #   so
    #     sqrt(Qp) == z * sqrt(omega / 2 / D) == z * Dtau
    frq <- Dtau.**2  * z2
  }
  return(frq)
}

#' Calculate volume of fluids in the sensing region of the borehole.
#'
#' This function calculates the volume of fluid in the screened section, 
#' namely \strong{Equation 2} in Kitagawa et al (2011).
#' 
#' Although typical scientific boreholes with water-level sensors are 
#' drilled very deeply, pore-fluids are only allowed to flow through
#' a relatively short section, known as the "screened" section.  The
#' calculation assumes two pairs of radii and lengths: one for the cemented (grout)
#' section, and another for the screened section.
#' 
#' The volume calculated is
#' \deqn{
#' \pi R_C^2 (L_C - L_S) + \pi R_S^2 L_S
#' }
#' where 
#' \eqn{R} and \eqn{L} denote radius and length respectively, and subscripts
#' \eqn{C} and \eqn{S} denote the cemented and screened sections respectively.
#' 
#' This calculation assumes the measurement is for a sealed well.
#'
#' @name sensing_volume
#' @export
#' 
#' @param rad_grout   radius of the grouting  \eqn{[m]}
#' @param len_grout   length of the grouting  \eqn{[m]}
#' @param rad_screen  radius of the screened interval  \eqn{[m]}
#' @param len_screen  length of the screened interval  \eqn{[m]}
#' 
#' @return scalar, with units of \eqn{[m^3]}
#' 
#' @author A. J. Barbour <andy.barbour@@gmail.com>
#' @family utilities
#' @seealso \code{\link{well_response}}
#' 
#' @examples
#' #### dummy example
#' sensing_volume(1, 1, 1, 1)
#' #
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
#' sensing_volume(Rc, Lc, Rs, Ls) # m**3, ~= 1.8
sensing_volume <- function(rad_grout, len_grout, rad_screen, len_screen){
  Rc. <- rad_grout
  Lc. <- len_grout
  Rs. <- rad_screen
  Ls. <- len_screen
  Vw. <- Rc. * Rc. * (Lc. - Ls.) + Rs. * Rs. * Ls.
  return(pi * Vw.)
}

#' Add proper logarithm ticks to a plot axis.
#' 
#' @details 
#' This uses \code{\link{pretty}} with \code{n==5}, and assumes 
#' that the data along the axis \code{ax} has 
#' \emph{already} been transformed into its logarithm.
#' \emph{Only integer exponents are labeled.}
#' 
#' The functions 
#'  \code{\link{log_ticks}},
#'  \code{\link{log2_ticks}}, and
#'  \code{\link{log10_ticks}} are wrapper functions.
#' 
#' Set the \code{axt} parameter (e.g. \code{xaxt}) to \code{'n'} 
#' in the original plot command to prevent adding default tick marks.
#' 
#' @rdname logticks
#' @export
#' @param ax numeric; the axis number to add tick-marks to
#' @param n.minor numeric; the number of minor ticks to display
#' @param t.lims numeric; the upper and lower tick limits (in log space)
#' @param t.ratio numeric; the ratio of minor to major tick lengths.
#' @param major.ticks numeric; the axis limits.
#' @param base numeric; the base of the logarithm (somewhat experimental)
#' @param ... additional parameters passed to the \code{axis} call for the major ticks.
#' 
#' @author A. J. Barbour <andy.barbour@@gmail.com>
#' @references 
#' This was modified from a post on StackOverflow:
#' \url{http://stackoverflow.com/questions/6955440/displaying-minor-logarithmic-ticks-in-x-axis-in-r}
#' @family PlotUtilities
#' @examples
#' x <- 10^(0:8)
#' y <- 1:9
#' plot(log10(x),y,xaxt="n",xlab="x",xlim=c(0,9))
#' logticks()
logticks <- function(ax=1, n.minor=9, t.lims, t.ratio=0.5, major.ticks=NULL, base=c("ten","ln","two"), ...){
  # axis limits
  lims <- par("usr")
  # x or y axis
  if (ax %in% c(1,3)) lims <- lims[1:2] else lims[3:4]
  # prettify
  if (is.null(major.ticks)) major.ticks <- unique(as.integer(pretty(lims, n=5)))
  
  if(missing(t.lims)) t.lims <- range(major.ticks)
  
  ml <- t.lims[1]
  mu <- t.lims[2]
  
  major.ticks <- major.ticks[major.ticks >= ml & major.ticks <= mu]
  
  # setup logarithm base and function
  base <- match.arg(base)
  LOG <- switch(base, ten=log10, ln=log, two=log2)
  base <- switch(base, two=2, ln="e", ten=10)
  
  # create tick label expressions
  tick.labels <- sapply(major.ticks, function(i) as.expression(bquote(.(base) ^ .(i))) )
  # add major ticks
  axis(ax, at=major.ticks, labels=tick.labels, ...)
  
  if (base=="e") base <- exp(1)
  n <- n.minor+2
  minors <- LOG(pretty(base ^ major.ticks[1:2], n)) - major.ticks[1]
  minors <- minors[-c(1,n)]
  
  minor.ticks <- c(outer(minors, major.ticks,`+`))
  minor.ticks.loc <- minor.ticks[minor.ticks > ml & minor.ticks < mu]
  
  # add minor ticks
  axis(ax, at=minor.ticks.loc, tcl=par("tcl")*t.ratio, labels=FALSE)
  
  return(invisible(base))
}
#' @rdname logticks
#' @export
log_ticks <- function(...) logticks(base="ln", ...)
#' @rdname logticks
#' @export
log2_ticks <- function(...) logticks(base="two", ...)
#' @rdname logticks
#' @export
log10_ticks <- function(...) logticks(base="ten", ...)