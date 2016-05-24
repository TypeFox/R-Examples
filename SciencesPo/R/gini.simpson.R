#' @encoding UTF-8
#' @title Gini-Simpson Index
#'
#' @description Computes the Gini/Simpson coefficient. \code{NA}s from the data are omitted.
#'
#' @param x A data.frame, a matrix-like, or a vector.
#' @param na.rm A logical value to deal with NAs.
#' @param \dots Additional arguements (currently ignored)
#'
#' @details The Gini-Simpson quadratic index is a classic measure of diversity, widely used by social scientists and ecologists. The Gini-Simpson is also known as Gibbs-Martin index in sociology, psychology and management studies, which in turn is also known as the Blau index. The Gini-Simpson index is computed as \eqn{1 - \lambda = 1 - \sum_{i=1}^R p_i^2 = 1 - 1/{}^2D}.
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}.
#' @keywords Diversity, Concentration, Inequality
#' @importFrom stats na.omit
#' @seealso \code{\link{politicalDiversity}}.
#' @examples
#' # generate a vector (of incomes)
#' x <- as.table(c(69,50,40,22))
#' rownames(x) <- c("AB","C","D","E")
#' gini.simpson(x)
#'
#' @rdname gini.simpson
#' @export
`gini.simpson` <- function(x, na.rm=TRUE, ...) UseMethod("gini.simpson")

#' @rdname gini.simpson
#' @export
`gini.simpson` <- function(x, na.rm = TRUE){
  # reference: Sachs, Angewandte Statistik, S. 57
  if(na.rm) x <- na.omit(x)
  x <- as.table(x)
  ptab <- prop.table(x)
  idx <-sum(ptab*(1-ptab))
  print(idx, digits = max(3, getOption("digits") - 3))
}##--end of gini.simpson
NULL




#' @encoding UTF-8
#' @title Weighted Gini Index
#'
#' @description Computes the unweighted and weighted Gini index of a distribution.
#'
#' @param x A data.frame, a matrix-like, or a vector.
#' @param weights A vector containing weights for \code{x}.
#' @param \dots Additional arguements (currently ignored)
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#' .
#' @keywords Diversity, Concentration
#' @seealso \code{\link{gini.simpson}}.
#' @examples
#' # generate a vector (of incomes)
#' x <- c(778, 815, 857, 888, 925, 930, 965, 990, 1012)
#' # compute Gini index
#' gini(x)
#'
#' gini(c(100,0,0,0))
#'
#' @rdname gini
#' @export
`gini` <- function(x, weights, ...) UseMethod("gini")


#' @export
#' @rdname gini
`gini` <- function(x, weights = rep(1, length = length(x)), ...){
  # TODO : add unbiased estiamtes see: gini(c(100,0,0,0))
  ox <- order(x)
  x <- x[ox]
  weights <- weights[ox]/sum(weights)
  p <- cumsum(weights)
  nu <- cumsum(weights * x)
  n <- length(nu)
  nu <- nu/nu[n]
  idx <- sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
  print(idx, digits = max(3, getOption("digits") - 3))
  # print(paste("Gini Index:", idx*100," Gini Coefficient:", idx))
}
NULL

# http://stats.stackexchange.com/questions/68940/how-is-the-weighted-gini-criterion-defined



# Original Zeileis:
# Gini <- function(x)
# {
#   n <- length(x)
#   x <- sort(x)
#   G <- sum(x * 1:n)
#   G <- 2*G/(n*sum(x))
#   G - 1 - (1/n)
# }

# other:
# http://rss.acs.unt.edu/Rdoc/library/reldist/html/gini.html
# http://finzi.psych.upenn.edu/R/library/dplR/html/gini.coef.html


# Gini <- function(x, n = rep(1, length(x)), unbiased = TRUE, conf.level = NA, R = 1000, type = "bca", na.rm = FALSE) {

#  x <- rep(x, n)    # same handling as Lc
#  if(na.rm) x <- na.omit(x)
#  if (any(is.na(x)) || any(x < 0)) return(NA_real_)
#
#  i.gini <- function (x, unbiased = TRUE){
#    n <- length(x)
#    x <- sort(x)
#
#    res <- 2 * sum(x * 1:n) / (n*sum(x)) - 1 - (1/n)
#    if(unbiased) res <- n / (n - 1) * res
#
#    # limit Gini to 0 here, if negative values appear, which is the case with
#    # Gini( c(10,10,10))
#    return( pmax(0, res))
#
#    # other guy out there:
#    #     N <- if (unbiased) n * (n - 1) else n * n
#    #     dsum <- drop(crossprod(2 * 1:n - n - 1, x))
#    #     dsum / (mean(x) * N)
#    # is this slower, than above implementation??
#  }
#
#  if(is.na(conf.level)){
#    res <- i.gini(x, unbiased = unbiased)
#
#  } else {
#    # adjusted bootstrap percentile (BCa) interval
#    boot.gini <- boot(x, function(x, d) i.gini(x[d], unbiased = unbiased), R=R)
#    ci <- boot.ci(boot.gini, conf=conf.level, type=type)
#    res <- c(gini=boot.gini$t0, lwr.ci=ci[[4]][4], upr.ci=ci[[4]][5])
#  }
#
#  return(res)
#
#}




#' @encoding UTF-8
#' @title The Lorenz Curve
#'
#' @description Computes the (empirical) ordinary and generalized Lorenz curve of a vector.
#'
#' @param x A vector of non-negative values.
#' @param n A vector of frequencies of the same length as \code{x}.
#' @param plot A logical. If TRUE the Lorenz curve will be plotted.
#' @param \dots Additional arguements (currently ignored)
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}.
#'
#' @details The Gini coefficient ranges from a minimum value of zero, when all individuals are equal, to a theoretical maximum of one in an infinite population in which every individual except one has a size of zero. It has been shown that the sample Gini coefficients originally defined need to be multiplied by n/(n-1) in order to become unbiased estimators for the population coefficients.
#'
#' @keywords Diversity, Concentration
#' @seealso \code{\link{gini}}, \code{\link{gini.simpson}}.
#' @examples
#' # generate a vector (of incomes)
#' x <- c(778, 815, 857, 888, 925, 930, 965, 990, 1012)
#' # compute Lorenz values
#' lorenz(x)
#' # generate some weights:
#' wgt <- runif(n=length(x))
#' # compute the lorenz with especific weights
#' lorenz(x, wgt)
#'
#' @rdname lorenz
#' @export
`lorenz` <- function(x, n = rep(1, length(x)), plot = FALSE, ...) UseMethod("lorenz")

#' @export
#' @rdname lorenz
`lorenz` <- function(x, n = rep(1, length(x)), plot = FALSE, ...)
{
  ina <- !is.na(x)
  n <- n[ina]
  x <- as.numeric(x)[ina]
  k <- base::length(x)
  o <- base::order(x)
  x <- x[o]
  n <- n[o]
  x <- n * x
  p <- base::cumsum(n)/sum(n)
  L <- base::cumsum(x)/sum(x)
  p <- c(0, p)
  L <- c(0, L)
  L2 <- L * base::mean(x)/mean(n)
  idx <- list(p, L, L2)
  names(idx) <- c("p", "L", "L.general")
  class(idx) <- "lorenz"
  if (plot)
    graphics::plot(idx)
  print(idx, digits = max(3, getOption("digits") - 3))
}##-end of Lorenz


# gini <- function(x, unbiased = TRUE, na.rm = FALSE){
#  if (!is.numeric(x)){
#    warning("'x' is not numeric; returning NA")
#    return(NA)
#  }
#  if (!na.rm && any(na.ind <- is.na(x)))
#    stop("'x' contain NAs")
#  if (na.rm)
#    x <- x[!na.ind]
#  n <- length(x)
#  mu <- mean(x)
#  N <- if (unbiased) n * (n - 1) else n * n
#  ox <- x[order(x)]
#  dsum <- drop(crossprod(2 * 1:n - n - 1,  ox))
#  dsum / (mu * N)
#}
