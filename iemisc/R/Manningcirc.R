#' Circular cross-section using the Gauckler-Manning-Strickler equation
#'
#' Manningcirc and Manningcircy solve for a missing variable for a circular
#' cross-section.
#'
#' The Manningcirc function solves for one missing variable in the Gauckler-
#' Manning equation for a circular cross-section and uniform flow. The
#' possible inputs are Q, n, Sf, y, and d. If y or d are not initially known,
#' then Manningcircy can solve for y or d to use as input in the Manningcirc
#' function.
#'
#' The Manningcircy function solves for one missing variable in the Gauckler-
#' Manning equation for a circular cross-section and uniform flow. The possible
#' inputs are y, d, y_d (ratio of y/d), and theta.
#'
#'
#'
#'
#' Gauckler-Manning-Strickler equation is expressed as
#'
#' \deqn{V = \frac{K_n}{n}R^\frac{2}{3}\sqrt{S}}
#'
#' \describe{
#'	\item{\emph{V}}{the velocity (m/s or ft/s)}
#'	\item{\emph{n}}{Manning's roughness coefficient (dimensionless)}
#'	\item{\emph{R}}{the hydraulic radius (m or ft)}
#'	\item{\emph{S}}{the slope of the channel bed (m/m or ft/ft)}
#'	\item{\emph{\eqn{K_n}}}{the conversion constant -- 1.0 for SI and
#'        3.2808399 ^ (1 / 3) for English units -- m^(1/3)/s or ft^(1/3)/s}
#' }
#'
#'
#'
#'
#' This equation is also expressed as
#'
#' \deqn{Q = \frac{K_n}{n}\frac{A^\frac{5}{3}}{P^\frac{2}{3}}\sqrt{S}}
#'
#' \describe{
#'	\item{\emph{Q}}{the discharge [m^3/s or ft^3/s (cfs)] is VA}
#'	\item{\emph{n}}{Manning's roughness coefficient (dimensionless)}
#'	\item{\emph{P}}{the wetted perimeter of the channel (m or ft)}
#'	\item{\emph{A}}{the cross-sectional area (m^2 or ft^2)}
#'	\item{\emph{S}}{the slope of the channel bed (m/m or ft/ft)}
#'	\item{\emph{\eqn{K_n}}}{the conversion constant -- 1.0 for SI and
#'        3.2808399 ^ (1 / 3) for English units -- m^(1/3)/s or ft^(1/3)/s}
#' }
#'
#'
#'
#'
#' Other important equations regarding the circular cross-section follow:
#' \deqn{R = \frac{A}{P}}
#'
#' \describe{
#'	\item{\emph{R}}{the hydraulic radius (m or ft)}
#'	\item{\emph{A}}{the cross-sectional area (m^2 or ft^2)}
#'	\item{\emph{P}}{the wetted perimeter of the channel (m or ft)}
#' }
#'
#'
#'
#'
#' \deqn{A = \left(\theta - \sin \theta\right) \frac{d^2}{8}}
#'
#' \describe{
#'	\item{\emph{A}}{the cross-sectional area (m^2 or ft^2)}
#'	\item{\emph{d}}{the diameter of the cross-section (m or ft)}
#'	\item{\emph{\eqn{\theta}}}{see the equation defining this parameter}
#' }
#'
#'
#'
#'
#' \deqn{\theta = 2 \arcsin\left[1 - 2\left(\frac{y}{d}\right)\right]}
#' \describe{
#'	\item{\emph{\eqn{\theta}}}{see the equation defining this parameter}
#'	\item{\emph{y}}{the flow depth (normal depth in this function) [m or ft]}
#'	\item{\emph{d}}{the diameter of the cross-section (m or ft)}
#' }
#'
#'
#'
#'
#' \deqn{d = 1.56 \left[\frac{nQ}{K_n\sqrt{S}}\right]^\frac{3}{8}}
#' \describe{
#'	\item{\emph{d}}{the initial diameter of the cross-section [m or ft]}
#'	\item{\emph{Q}}{the discharge [m^3/s or ft^3/s (cfs)] is VA}
#'	\item{\emph{n}}{Manning's roughness coefficient (dimensionless)}
#'	\item{\emph{S}}{the slope of the channel bed (m/m or ft/ft)}
#'	\item{\emph{\eqn{K_n}}}{the conversion constant -- 1.0 for SI and
#'        3.2808399 ^ (1 / 3) for English units -- m^(1/3)/s or ft^(1/3)/s}
#' }
#'
#' Note: This will only provide the initial conduit diameter, check the design
#'       considerations to determine your next steps.
#'
#'
#'
#'
#' \deqn{P = \frac{\theta d}{2}}
#'
#' \describe{
#'	\item{\emph{P}}{the wetted perimeter of the channel (m or ft)}
#'	\item{\emph{\eqn{\theta}}}{see the equation defining this parameter}
#'	\item{\emph{d}}{the diameter of the cross-section (m or ft)}
#' }
#'
#'
#'
#'
#' \deqn{B = d \sin\left(\frac{\theta}{2}\right)}
#'
#' \describe{
#'	\item{\emph{B}}{the top width of the channel (m or ft)}
#'	\item{\emph{\eqn{\theta}}}{see the equation defining this parameter}
#'	\item{\emph{d}}{the diameter of the cross-section (m or ft)}
#' }
#'
#'
#'
#'
#' Assumptions: uniform flow and prismatic channel
#'
#' Note: Units must be consistent
#'
#'
#' @param Q numeric vector that contains the discharge value [m^3/s or ft^3/s],
#'   if known.
#' @param n numeric vector that contains the Manning's roughness coefficient n,
#'   if known.
#' @param d numeric vector that contains the diameter value (m or ft),
#'   if known.
#' @param Sf numeric vector that contains the the bed slope (m/m or ft/ft),
#'   if known.
#' @param y numeric vector that contains the flow depth (m or ft), if known.
#' @param y_d numeric vector that contains the filling ration (y/d), if known.
#' @param theta numeric vector that contains the angle theta (radians), if
#'        known.
#' @param units character vector that contains the system of units [options are
#'   \code{SI} for International System of Units and \code{Eng} for English units
#'   (United States Customary System in the United States and Imperial Units in
#'   the United Kingdom)]
#'
#' @return the missing parameter (Q, n, or Sf) & theta, area (A), wetted
#'   perimeter (P), top width (B), and R (hydraulic radius) as a \code{\link[base]{list}}
#'   for the Manningcirc function.
#' @return the missing parameter (d or y) & theta, area (A), wetted
#'   perimeter (P), top width (B), and R (hydraulic radius) as a \code{\link[base]{list}}
#'   for the Manningcircy function.
#'
#'
#' @source
#' r - Better error message for stopifnot? - Stack Overflow answered by Andrie on Dec 1 2011. See \url{http://stackoverflow.com/questions/8343509/better-error-message-for-stopifnot}.
#'
#'
#' @references
#' \enumerate{
#'    \item Terry W. Sturm, \emph{Open Channel Hydraulics}, 2nd Edition, New York City, New York: The McGraw-Hill Companies, Inc., 2010, page 8, 36, 102, 120, 123-125, 153-154.
#'    \item Dan Moore, P.E., NRCS Water Quality and Quantity Technology Development Team, Portland Oregon, "Using Mannings Equation with Natural Streams", August 2011, \url{http://www.wcc.nrcs.usda.gov/ftpref/wntsc/H&H/xsec/manningsNaturally.pdf}.
#'    \item Gilberto E. Urroz, Utah State University Civil and Environmental Engineering, CEE6510 - Numerical Methods in Civil Engineering, Spring 2006,  "Solving selected equations and systems of equations in hydraulics using Matlab", August/September 2004, \url{http://ocw.usu.edu/Civil_and_Environmental_Engineering/Numerical_Methods_in_Civil_Engineering/}.
#'    \item Tyler G. Hicks, P.E., \emph{Civil Engineering Formulas: Pocket Guide}, 2nd Edition, New York City, New York: The McGraw-Hill Companies, Inc., 2002, page 423, 425.
#'    \item Wikimedia Foundation, Inc. Wikipedia, 26 November 2015, “Manning formula”, \url{https://en.wikipedia.org/wiki/Manning_formula}.
#' }
#'
#' @encoding UTF-8
#'
#'
#'
#'
#' @seealso \code{\link{Manningtrap}} for a trapezoidal cross-section, \code{\link{Manningrect}} for a
#'   rectangular cross-section, \code{\link{Manningtri}} for a triangular cross-section,
#'   and \code{\link{Manningpara}} for a parabolic cross-section.
#'
#'
#'
#'
#' @examples
#' library(iemisc)
#' library(iemiscdata)
#' \dontrun{
#' # The equations used in this function were solved analytically after
#' # following these steps:
#' library(rSymPy) # review the package to determine its system dependencies
#' Q <- Var("Q")
#' n <- Var("n")
#' k <- Var("k")
#' Sf <- Var("Sf")
#' y <- Var("y")
#' d <- Var("d")
#' theta <- Var("theta")
#'
#' # Simplify with rSymPy
#' eqsimp <- sympy("expr = n*Q*((theta*d)/2)**(2/3)-k*(theta-sin(theta))*(d**2/8)**(5/3)*sqrt(Sf)")
#' # eqsimp is "Q*n - k*Sf**(1/2)*d**2*(theta - sin(theta))/8"
#' # This is the equation that was used to solve for the missing variables
#' }
#'
#'
#'
#'
#' # Example 4.1 from Sturm (page 124-125)
#' Manningcircy(y_d = 0.8, d = 2, units = "Eng")
#' y <- Manningcircy(y_d = 0.8, d = 2, units = "Eng")
#' # defines all list values within the object named y
#' y$y # gives the value of y
#'
#'
#'
#'
#'
#' @name Manningcirc
NULL

#' @export
#' @rdname Manningcirc
Manningcirc <- function (Q = NULL, n = NULL, Sf = NULL, y = NULL, d = NULL, units = c("SI", "Eng")) {

checks <- c(Q, n, Sf, y, d)
units <- units

if (length(checks) < 4) {

stop("There are not at least 4 known variables. Try again with at least 4 known variables.")
# Source 1 / only process enough known variables and provide a stop warning if not enough

} else {

if (any(checks == 0)) {

stop("Either Q, n, Sf, y, or d is 0. None of the variables can be 0. Try again.")
# Source 1 / only process with a non-zero value for Q, n, Sf, y, and d and provide a stop warning if Q, n, Sf, y, or d = 0

} else {

if (units == "SI") {

   k <- 1

} else if (units == "Eng") {

   k <- 3.2808 ^ (1 / 3)

} else if (all(c("SI", "Eng") %in% units == FALSE) == FALSE) {

stop("Incorrect unit system. Try again.")
# Source 1 / only process with a specified unit and provide a stop warning if not

}

if (missing(Q)) {

theta <- 2 * acos(1 - (2 * (y / d)))
A <- (theta - sin(theta)) * (d ^ 2 / 8)
P <- ((theta * d) / 2)
B <- d * sin(theta / 2)
R <- A / P

Q <- (k / n) * sqrt(Sf) * A * R ^ (2/3)

return(list(Q = Q, theta = theta, A = A, P = P, B = B, R = R))

} else if (missing(n)) {

theta <- 2 * acos(1 - (2 * (y / d)))
A <- (theta - sin(theta)) * (d ^ 2 / 8)
P <- ((theta * d) / 2)
B <- d * sin(theta / 2)
R <- A / P

n <- (k / Q) * sqrt(Sf) * A * R ^ (2/3)

return(list(n = n, theta = theta, A = A, P = P, B = B, R = R))

} else if (missing(Sf)) {

theta <- 2 * acos(1 - (2 * (y / d)))
A <- (theta - sin(theta)) * (d ^ 2 / 8)
P <- ((theta * d) / 2)
B <- d * sin(theta / 2)
R <- A / P

Sf <- ((Q * n) / ((k * A * R ^ (2/3)))) ^ 2

return(list(Sf = Sf, theta = theta, A = A, P = P, B = B, R = R))
}
}
}
}


#' @export
#' @rdname Manningcirc
Manningcircy <- function (y = NULL, d = NULL, y_d = NULL, theta = NULL, Sf = NULL, Q = NULL, units = c("SI", "Eng")) {

checks <- c(y, d, y_d, theta, Sf, Q)
units <- units

if (length(checks) < 1) {

stop("There are not at least 1 known variables. Try again with at least 1 known variables.")
# Source 1 / only process enough known variables and provide a stop warning if not enough

} else {

if (any(checks == 0)) {

stop("Either y, d, theta, y_d, Sf, or Q is 0. None of the variables can be 0. Try again.")
# Source 1 / only process with a non-zero value for y, d, theta, y_d, Sf, and Q and provide a stop warning if y, d, theta, y_d, Sf, or Q = 0

} else {

if (units == "SI") {

   k <- 1

} else if (units == "Eng") {

   k <- 3.2808 ^ (1 / 3)

} else if (all(c("SI", "Eng") %in% units == FALSE) == FALSE) {

stop("Incorrect unit system. Try again.")
# Source 1 / only process with a specified unit and provide a stop warning if not

}

if (missing(y) & missing(y_d)) {

y <- (d / 2) * (1 - cos(theta / 2))

theta <- 2 * acos(1 - (2 * (y / d)))
A <- (theta - sin(theta)) * (d ^ 2 / 8)
P <- ((theta * d) / 2)
B <- d * sin(theta / 2)
R <- A / P

return(list(y = y, theta = theta, A = A, P = P, B = B, R = R))

} else if (missing(y) & missing(y_d) & missing(theta)) {

rh <- (n * Q) / (k * sqrt(Sf))

f <- function (theta) ((theta - sin(theta)) * (d ^ 2 / 8)) * (((theta - sin(theta)) * (d ^ 2 / 8) / ((theta * d) / 2)) ^ (2 / 3)) - rh

theta <- uniroot(f, c(-1000, 1000))
theta <- theta$root

y <- (d / 2) * (1 - cos(theta / 2))

A <- (theta - sin(theta)) * (d ^ 2 / 8)
P <- ((theta * d) / 2)
B <- d * sin(theta / 2)
R <- A / P

return(list(y = y, theta = theta, A = A, P = P, B = B, R = R))

} else if (missing(d)) {

if (missing(Q) & missing(Sf))

stop("Q and Sf are needed to compute d. Try again with a value for Q and Sf.")
# Source 1 / only process enough known variables and provide a stop warning if not enough

d <- 1.56 * ((n * Q) / (k * sqrt(Sf))) ^ (3 / 8)

theta <- 2 * acos(1 - (2 * (y / d)))
A <- (theta - sin(theta)) * (d ^ 2 / 8)
P <- ((theta * d) / 2)
B <- d * sin(theta / 2)
R <- A / P

return(list(d = d, theta = theta, A = A, P = P, B = B, R = R))

} else if (missing(theta) & missing(y)) {

theta <- 2 * acos(1 - (2 * (y_d)))
y <- (d / 2) * (1 - cos(theta / 2))

A <- (theta - sin(theta)) * (d ^ 2 / 8)
P <- ((theta * d) / 2)
B <- d * sin(theta / 2)
R <- A / P

return(list(theta = theta, y = y, A = A, P = P, B = B, R = R))

} else if (missing(theta)) {

theta <- 2 * acos(1 - (2 * (y / d)))

A <- (theta - sin(theta)) * (d ^ 2 / 8)
P <- ((theta * d) / 2)
B <- d * sin(theta / 2)
R <- A / P

return(list(theta = theta, A = A, P = P, B = B, R = R))
}
}
}
}
