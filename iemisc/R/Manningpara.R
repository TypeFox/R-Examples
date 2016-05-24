#' Parabolic cross-section for the Gauckler-Manning-Strickler equation
#'
#' This function solves for one missing variable in the Gauckler-Manning-
#' Strickler equation for a parabolic cross-section and uniform flow.
#' This function is NOT operable yet.
#'
#'
#'
#'
#' Gauckler-Manning-Strickler equation is expressed as
#'
#' \deqn{V = \frac{K_n}{n}R^\frac{2}{3}S^\frac{1}{2}}
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
#' \deqn{Q = \frac{K_n}{n}\frac{A^\frac{5}{3}}{P^\frac{2}{3}}S^\frac{1}{2}}
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
#' Other important equations regarding the parabolic cross-section follow:
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
#' \deqn{A = \frac{2}{3}By}
#'
#' \describe{
#'	\item{\emph{A}}{the cross-sectional area (m^2 or ft^2)}
#'	\item{\emph{y}}{the flow depth (normal depth in this function) [m or ft]}
#'	\item{\emph{B}}{the top width of the channel (m or ft)}
#' }
#'
#'
#'
#'
#' \deqn{P = \left(\frac{B}{2}\right)\left[\sqrt{\left(1 + x^2\right)} + \left(\frac{1}{x}\right) \ln \left(x + \sqrt{\left(1 + x^2\right)}\right)\right]}
#'
#' \describe{
#'	\item{\emph{P}}{the wetted perimeter of the channel (m or ft)}
#'	\item{\emph{B}}{the top width of the channel (m or ft)}
#'	\item{\emph{x}}{4y/b (dimensionless)}
#' }
#'
#'
#'
#'
#' \deqn{x = \frac{4y}{B}}
#'
#' \describe{
#'	\item{\emph{x}}{4y/b (dimensionless)}
#'	\item{\emph{B}}{the top width of the channel (m or ft)}
#'	\item{\emph{y}}{the flow depth (normal depth in this function) [m or ft]}
#' }
#'
#'
#'
#'
#' \deqn{B = B_1 \left(\sqrt{\frac{y}{y_1}}\right)}
#'
#' \describe{
#'	\item{\emph{B}}{the top width of the channel (m or ft)}
#'	\item{\emph{y}}{the flow depth (normal depth in this function) [m or ft]}
#'	\item{\emph{\eqn{B_1}}}{the "bank-full width" (m or ft)}
#'	\item{\emph{\eqn{y_1}}}{the "bank-full depth" (m or ft)}
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
#' @param m numeric vector that contains the "cross-sectional side slope of m:1
#'   (horizontal:vertical)", if known.
#' @param Sf numeric vector that contains the the bed slope (m/m or ft/ft),
#'   if known.
#' @param y numeric vector that contains the flow depth (m or ft), if known.
#' @param B1 numeric vector that contains the the "bank-full width", if known.
#' @param y1 numeric vector that contains the "bank-full depth", if known.
#' @param units character vector that contains the system of units [options are
#'   \code{SI} for International System of Units and \code{Eng} for English units
#'   (United States Customary System in the United States and Imperial Units in
#'   the United Kingdom)]
#'
#' @return the missing parameter (Q, n, m, Sf, B1, y1, or y) & area (A), wetted
#'   perimeter (P), top width (B), and R (hydraulic radius) as a \code{\link[base]{list}}.
#'
#'
#' @source
#' r - Better error message for stopifnot? - Stack Overflow answered by Andrie on Dec 1 2011. See \url{http://stackoverflow.com/questions/8343509/better-error-message-for-stopifnot}.
#'
#'
#' @references
#' \enumerate{
#'    \item Terry W. Sturm, \emph{Open Channel Hydraulics}, 2nd Edition, New York City, New York: The McGraw-Hill Companies, Inc., 2010, page 8, 36, 102, 120, 153-154.
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
#'   and \code{\link{Manningcirc}} for a circular cross-section.
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
#' m <- Var("m")
#' k <- Var("k")
#' Sf <- Var("Sf")
#' y <- Var("y")
#'
#' # Simplify with rSymPy
#' eqsimp <- sympy("expr = n*Q*(2*y*sqrt(1+m**2))**(2/3)-k*(m*y**2)**(5/3)*sqrt(Sf)")
#' # eqsimp is "Q*n - k*m*Sf**(1/2)*y**2"
#' # This is the equation that was used to solve for the missing variables
#' }
#'
#'
#' # Modified Exercise 4.1 from Sturm (page 153)
#' Manningpara(Q = 3000, m = 3, Sf = 0.002, n = 0.025, units = "Eng")
#' # Q = 3000 cfs, m = 3, Sf = 0.002 ft/ft, n = 0.025, units = English units
#' # This will solve for y since it is missing and y will be in ft
#'
#'
#' # Modified Exercise 4.5 from Sturm (page 154)
#' Manningpara(Q = 950, m = 2, Sf = 0.022, n = 0.023, units = "SI")
#' # Q = 950 m^3/s, m = 2, Sf = 0.022 m/m, n = 0.023, units = SI units
#' # This will solve for y since it is missing and y will be in m
#'
#' @export
Manningpara <- function (Q = NULL, n = NULL, m = NULL, Sf = NULL, y = NULL, B1 = NULL, y1 = NULL, units = c("SI", "Eng")) {

checks <- c(Q, n, m, Sf, y, B1, y1)
units <- units

if (length(checks) < 4) {

stop("There are not at least 4 known variables. Try again with at least 4 known variables.")
# Source 1 / only process enough known variables and provide a stop warning if not enough

} else {

if (any(checks == 0)) {

stop("Either Q, n, m, Sf, or y is 0. None of the variables can be 0. Try again.")
# Source 1 / only process with a non-zero value for Q, n, m, Sf, and y and provide a stop warning if Q, n, m, Sf, or y = 0

} else {

if (units == "SI") {

   k <- 1

} else if (units == "Eng") {

   k <- 3.2808399 ^ (1 / 3)

} else if (all(c("SI", "Eng") %in% units == FALSE) == FALSE) {

stop("Incorrect unit system. Try again.")
# Source 1 / only process with a specified unit and provide a stop warning if not

}

if (missing(Q)) {

B <- B1 * sqrt(y / y1)
x <- (4 * y) / B
A <- (2 / 3) * B * y
P <- (B / 2) * (sqrt(1 + x ^ 2) + (1 / x) * log(x + (sqrt(1 + x ^ 2))))
R <- A / P

Q <- (k / n) * sqrt(Sf) * y ^ 2 * m

return(list(Q = Q, A = A, P = P, B = B, R = R))

} else if (missing(n)) {

B <- B1 * sqrt(y / y1)
x <- (4 * y) / B
A <- (2 / 3) * B * y
P <- (B / 2) * (sqrt(1 + x ^ 2) + (1 / x) * log(x + (sqrt(1 + x ^ 2))))
R <- A / P

n <- (k / Q) * m * sqrt(Sf) * y ^ 2

return(list(n = n, A = A, P = P, B = B, R = R))

} else if (missing(m)) {

B <- B1 * sqrt(y / y1)
x <- (4 * y) / B
A <- (2 / 3) * B * y
P <- (B / 2) * (sqrt(1 + x ^ 2) + (1 / x) * log(x + (sqrt(1 + x ^ 2))))
R <- A / P

return(list(m = m, A = A, P = P, B = B, R = R))

} else if (missing(y)) {

y <- sqrt((Q * n) / (k * m * sqrt(Sf)))

B <- B1 * sqrt(y / y1)
x <- (4 * y) / B
A <- (2 / 3) * B * y
P <- (B / 2) * (sqrt(1 + x ^ 2) + (1 / x) * log(x + (sqrt(1 + x ^ 2))))
R <- A / P

return(list(y = y, A = A, P = P, B = B, R = R))

} else if (missing(Sf)) {

B <- B1 * sqrt(y / y1)
x <- (4 * y) / B
A <- (2 / 3) * B * y
P <- (B / 2) * (sqrt(1 + x ^ 2) + (1 / x) * log(x + (sqrt(1 + x ^ 2))))
R <- A / P

Sf <- ((Q * n) / (k * m * y ^ 2)) ^ 2

return(list(Sf = Sf, A = A, P = P, B = B, R = R))
}
}
}
}
