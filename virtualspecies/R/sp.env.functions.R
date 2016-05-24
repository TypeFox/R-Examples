#' Linear function
#' @description A simple linear function of the form
#' \deqn{ax+b}{a*x+b}
#' @param x a numeric value or vector
#' @param a a numeric value or vector
#' @param b a numeric value or vector
#' @return a numeric value or vector resulting from the function
#' @export
#' @seealso \code{\link{logisticFun}}, \code{\link{quadraticFun}}
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#'
#' Maintainer: Boris Leroy \email{leroy.boris@@gmail.com}
#' @examples
#' x <- 1:100
#' y <- linearFun(x, a = 0.5, b = 0)
#' plot(y ~ x, type = "l")
linearFun <- function(x, a, b) {a * x + b}

#' Logistic function
#' 
#' @description A simple logistic function of the form
#' \deqn{\frac{1}{{1 + e^{\frac{x - \beta}{\alpha}}}}}{
#' 1 / (1 + exp((x - \beta)/\alpha))}
#' @param x a numeric value or vector
#' @param alpha a numeric value or vector
#' @param beta a numeric value or vector
#' @return a numeric value or vector resulting from the function
#' @export
#' @seealso \code{\link{linearFun}}, \code{\link{quadraticFun}}
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#'
#' Maintainer: Boris Leroy \email{leroy.boris@@gmail.com}
#' @details
#' The value of \code{beta} determines the 'threshold' of the logistic curve 
#' (i.e. the inflexion point).
#' 
#' The value of \code{alpha} determines the slope of the curve (see examples):
#' \itemize{
#' \item{\code{alpha} very close to 0 will result in a threshold-like response.}
#' \item{Values of \code{alpha} with the same order of magnitude as the range of
#' \code{x} (e.g., the range of\code{x} / 10) will result in a 
#' logistic function.} 
#' \item{\code{alpha} very far from 0 will result in a linear function.}
#' }
#' @examples
#' x <- 1:100
#' y <- logisticFun(x, alpha = -10, b = 50)
#' plot(y ~ x, type = "l")
#' 
#' # The effect of alpha:
#' y1 <- logisticFun(x, alpha = -0.01, b = 50)
#' y2 <- logisticFun(x, alpha = -10, b = 50)
#' y3 <- logisticFun(x, alpha = -1000, b = 50)
#' 
#' par(mfrow = c(1, 3))
#' plot(y1 ~ x, type = "l", main = expression(alpha %->% 0))
#' plot(y2 ~ x, type = "l", main = expression(alpha %~~% range(x)/10))
#' plot(y3 ~ x, type = "l", main = expression(alpha %->% infinity))
logisticFun <- function(x, alpha, beta) {1 / (1 + exp((x - beta)/alpha))}

#' Quadratic function
#' 
#' @description A simple quadratic function of the form
#' \deqn{ax^2+bx+c}{
#' a*x^2+b*x+c}
#' @param x a numeric value or vector
#' @param a a numeric value or vector
#' @param b a numeric value or vector
#' @param c a numeric value or vector
#' @return a numeric value or vector resulting from the function
#' @export
#' @seealso \code{\link{linearFun}}, \code{\link{quadraticFun}}
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#'
#' Maintainer: Boris Leroy \email{leroy.boris@@gmail.com}
#' @examples
#' x <- 1:100
#' y <- quadraticFun(x, a = 2, b = 2, c = 3)
#' plot(y ~ x, type = "l")
quadraticFun <- function(x, a, b, c) {a * x^2 + b * x + c}


#' Normal function defined by extremes
#' 
#' @description A modified version of the normal function based on three parameters:
#' \itemize{
#' \item{the mean}
#' \item{the absolute difference between the mean and extreme values}
#' \item{the percentage of area under the curve between the specified extreme values}
#' }
#' 
#' See the example for an easier understanding.
#' @param x a numeric value or vector. The input environmental variable.
#' @param mean a numeric value or vector. The optimum (mean) of the normal curve
#' @param diff a numeric value or vector. The absolute difference between the mean and extremes.
#' @param prob a numeric value or vector. The percentage of the area under the curve between the 
#' chosen extreme values
#' @return a numeric value or vector resulting from the function
#' @export
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}, Florian David
#'
#' Maintainer: Boris Leroy \email{leroy.boris@@gmail.com}
#' @examples
#' # Let's define the response of a species to temperature which
#' #  - has an optimum at 20 degrees C
#' #  - occurs 99% of the time between 13 and 27 degrees C.
#' # In that case, mean = 20, diff = 7, and prob = 0.99
#' 
#' # First, we generate an arbitrary temperature variable 
#' # between 0 and 30 degrees C
#' temp <- seq(0, 30, length = 1000)
#' 
#' 
#' # Then, we calculate the response to this variable with the chosen values
#' response <- custnorm(x = temp, mean = 20, diff = 7, prob = .99)
#' 
#' plot(response ~ temp, type = "l")
custnorm <- function(x, mean, diff, prob)
{
  prob <- prob + (1 - prob)/2
  sd <- - diff / qnorm(p = 1 - prob)
  dnorm(x, mean = mean, sd = sd)
}

#' Beta response function
#' 
#' @description Generation of a beta response curve (see references) according to the equation:
#' \deqn{\frac{1}{{1 + e^{\frac{x - \beta}{\alpha}}}}}{
#' P = k (x - p1)^\alpha (p2 - x)^\gamma}
#' k is automatically estimated to have a maximum value of P equal to 1.
#' @param x a numeric value or vector. The input environmental variable.
#' @param p1 a numeric value or vector. Lower tolerance bound for the species
#' @param p2 a a numeric value or vector. Upper tolerance bound for the species
#' @param alpha a numeric value or vector. Parameter controlling the shape of the curve (see details)
#' @param gamma a numeric value or vector. Parameter controlling the shape of the curve (see details)
#' @return a numeric value or vector resulting from the function
#' @details
#' p1 and p2 can be seen as the upper and lower critical threshold of the curve.
#' \code{alpha} and \code{gamma} control the shape of the curve near p1 and p2, respectively.
#' When \code{alpha} = \code{gamma}, the curve is symetric. Low values of \code{alpha} and \code{gamma} 
#' result in smooth (< 1) to plateau (< 0.01) curves. Higher values result in 
#' peak (> 10) curves. 
#' 
#' When \code{alpha} < \code{gamma}, the curve is skewed to the right.
#' When \code{gamma} < \code{alpha}, the curve is skewed to the left.
#' @export
#' @references
#' Oksanen, J. & Minchin, P.R. (2002). Continuum theory revisited: what shape 
#' are species responses along ecological gradients? \emph{Ecological Modelling}
#' \bold{157}:119-129.
#' @seealso \code{\link{linearFun}}, \code{\link{quadraticFun}}, \code{\link{custnorm}}
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#'
#' Maintainer: Boris Leroy \email{leroy.boris@@gmail.com}
#' @examples
#' temp <- seq(-10, 40, length = 100)
#' # A curve similar to a thermal performance curve
#' P <- betaFun(x = temp, p1 = 0, p2 = 35, alpha = 0.9, gamma = 0.08)
#' plot(P ~ temp, type = "l")
betaFun <- function(x, p1, p2, alpha, gamma) 
{
  k <- 1/((alpha * (p2 - p1) / (alpha + gamma))^alpha) / ((gamma * (p2 - p1) / (alpha + gamma))^gamma)
  ifelse(x > p1 & x < p2, k * ((x - p1)^alpha) * (p2 - x)^gamma, 0)
}

#' Huisman-Olff-Fresco response function
#' 
##' @description A Huisman-Olff-Fresco response function:
#' \deqn{P = \frac{1}{{1 + e^{a + b x}}} \frac{1}{1 + e^{c - dx}}}{
#' P = (1 / (1 + exp(a + bx))) * (1 / (1 + exp(c -dx)))}
#' @param x a numeric value or vector
#' @param a a numeric value or vector
#' @param b a numeric value or vector
#' @param c a numeric value or vector
#' @return a numeric value or vector resulting from the function
#' @seealso \code{\link{linearFun}}, \code{\link{quadraticFun}}
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#'
#' Maintainer: Boris Leroy \email{leroy.boris@@gmail.com}
#' @examples
#' temp <- seq(-10, 40, length = 100)
#' # A curve similar to a thermal performance curve
#' P <- HOFFun(x = temp, a = -200, b = 10, c = 10, d = 0.1)
#' plot(P ~ temp, type = "l")
#' 
#' 
# .HOFFun <- function(x, a, b, c, d)
# {
#   if (a == 0 & b == 0)
#   {
#     stop("a and b can't both be set to zero")
#   } else if (c == 0 & d == 0)
#   {
#     stop("c and d can't both be set to zero")
#   }
#   M/(1 + exp(a + b * x)) * 1/(1 + exp(c - d * x))
#   (1 / (1 + exp(a + b * x))) * (1 / (1 + exp(c - d * x)))
# }


# Functions useful for the PCA approach

.f <- function(x, co) x %*% co

.pca.coordinates <- function(x, pca, na.rm)
{
  x <- sweep(x, 2L, pca$cent, check.margin=FALSE)
  x <- sweep(x, 2L, pca$norm, "/", check.margin=FALSE)
  x1 <- apply(x, 1, .f, co = pca$c1[, 1])
  x2 <- apply(x, 1, .f, co = pca$c1[, 2])
  return(cbind(x1, x2))
}

.prob.gaussian <- function(x, means, sds)
{
  dnorm(x[1], mean = means[1], sd = sds[1]) * dnorm(x[2], mean = means[2], sd = sds[2])
}


.thermalFun <- function(Pmax, Tb, To, rho, sigma)
{
  Pmax * exp(-exp(rho * (Tb - To) - 6) - sigma * (Tb - To)^2)
}
