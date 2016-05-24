
#' Random sample generated from an specified gaussian mixture model.
#' 
#' @param n Sample size
#' @param Pi A vector indicating the mixing proportions
#' @param Mu A two dimensional array where second component indicates the mean 
#' of each gaussian component
#' @param S A three dimensional array where third component indicates the variance 
#' of each gaussian component
#' @param labels A logical indicating whether or not a label shoud be returned indicating
#' the component from where observation has been generated (default FALSE)
#' @return A matrix with \code{n} row and columns given by the dimension of \code{Mu} and \code{S}. 
#' If \code{labels = T} another column is included indicating the component from where the observation 
#' was generated.
#' @examples
#' Pi = c(0.5, 0.3, 0.2)
#' Mu = array(c(## Mu first component
#'              5, 5,
#'              ## Mu second component
#'              1, 1,
#'              ## Mu third component
#'              0, 0), dim = c(2,3))
#' S = array(c(## Sigma first component
#'             1, 0,
#'             0, 1,
#'             ## Sigma second component
#'             0.2, 0,
#'             0, 0.2,
#'             ## Sigma third component
#'             0.05, 0,
#'             0, 0.05), dim = c(2, 2, 3))
#' X = rmixnorm(100, Pi = Pi, Mu = Mu, S = S, labels = TRUE)
#' plot(X[,1:2], col=X[,3])
#' @export
rmixnorm <- function(n, Pi, Mu, S, labels = F) {
  z <- sample(x = 1:length(Pi), size = n, prob = Pi, replace = TRUE)
  rmn <- matrix(0, nrow = n, ncol = nrow(Mu))
  for (i in 1:length(Pi)) {
    n_z <- sum(z == i)
    if (n_z != 0) {
      if (ncol(Mu) == 1) 
        rmn[z == i, ] <- stats::rnorm(n_z, mean = Mu[, i], sd = sqrt(S[, , i])) else rmn[z == i, ] <- mvtnorm::rmvnorm(n_z, mean = Mu[, i], sigma = S[, , i])
    }
  }
  if (labels) 
    rmn <- cbind(rmn, z)
  rmn
}
#' Random sample generated from an specified gaussian mixture model. 
#' 
#' The parameters are defined from the parameters obtained using other 
#' packages (' \code{Mclust}, 
#' \code{rmixmod})
#' 
#' @param n sample size
#' @param solution solution comming from packages \code{Mclust} or \code{rmixmod}
#' @param ... arguments passed to function \code{\link{rmixnorm}}
#' @export
#' @examples
#' require(mclust)
#' mod1 = Mclust(iris[,1:4])
#' rmixnorm_solution(10, mod1)
rmixnorm_solution <- function(n, solution, ...) {
  if ("MixmodCluster" %in% methods::is(solution)) {
    if (solution@dataType == "quantitative") {
      return(rmixnorm_rmixmod(n, solution, ...))
    }
  }
  if ("Mclust" %in% methods::is(solution)) {
    return(rmixnorm_mclust(n, solution, ...))
  }
  stop("Not recognized format for solution")
}
rmixnorm_mclust <- function(n, mclust_solution, ...) {
  func_pi <- mclust_solution$parameters$pro
  func_mean <- mclust_solution$parameters$mean
  func_sigma <- mclust_solution$parameters$variance$sigma
  rmixnorm(n, func_pi, func_mean, func_sigma, ...)
}
rmixnorm_rmixmod <- function(n, mixmod_solution, ...) {
  func_pi <- mixmod_solution@bestResult@parameters@proportions
  func_mean <- t(mixmod_solution@bestResult@parameters@mean)
  func_sigma <- array(do.call("c", mixmod_solution@bestResult@parameters@variance), c(NCOL(mixmod_solution@bestResult@parameters@mean),
                                                                                      NCOL(mixmod_solution@bestResult@parameters@mean),
                                                                                      NROW(mixmod_solution@bestResult@parameters@mean)))
  #func_sigma <- do.call("abind", list(mixmod_solution@bestResult@parameters@variance, along = 3))
  rmixnorm(n, func_pi, func_mean, func_sigma, ...)
}

#' Density function of specified gaussian mixture model.
#' 
#' @param x vector/matrix where density function is evaluated
#' @param Pi a vector indicating the mixing proportions
#' @param Mu a two dimensional array where second component indicates the mean 
#' of each gaussian component
#' @param S a three dimensional array where third component indicates the variance 
#' of each gaussian component
#' @param part subcomposition where x shoud be evaluated. Take into an account that
#' if x has dimensions K, K components must be selected by \code{part}
#' @param closure are probabilities Pi summing up to 1 (default TRUE)
#' @export
dmixnorm <- function(x, Pi, Mu, S, part = 1:length(Pi), closure = TRUE) {
  # z = sample(x=1:length(Pi), size=n, prob=Pi, replace=T) rmn = matrix(0, nrow=n, ncol=nrow(Mu))
  if (is.vector(x)) {
    dmn <- 0
  } else {
    dmn <- rep(0, times = nrow(x))
  }
  for (i in part) {
    if (ncol(Mu) == 1) 
      dmn <- dmn + Pi[i] * stats::dnorm(x, mean = Mu[, i], sd = sqrt(S[, , i])) else dmn <- dmn + Pi[i] * mvtnorm::dmvnorm(x, mean = Mu[, i], sigma = S[, , i])
  }
  if(closure){
    dmn/sum(Pi[part])
  }else{
    dmn
  }
}

#' Density function of specified gaussian mixture model.
#' 
#' The parameters are defined from the parameters obtained using other 
#' packages (' \code{Mclust}, 
#' \code{rmixmod})
#' 
#' @param x vector/matrix where density function is evaluated
#' @param solution solution comming from packages \code{Mclust} or \code{rmixmod}
#' @param ... arguments passed to function \code{\link{dmixnorm}}
#' @export
#' @examples
#' require(mclust)
#' mod1 = Mclust(iris[,1:4])
#' rmixnorm_solution(10, mod1)
dmixnorm_solution <- function(x, solution, ...) {
  if ("MixmodCluster" %in% methods::is(solution)) {
    if (solution@dataType == "quantitative") {
      return(dmixnorm_rmixmod(x, solution, ...))
    }
  }
  if ("Mclust" %in% methods::is(solution)) {
    return(dmixnorm_mclust(x, solution, ...))
  }
  stop("Not recognized format for solution")
}
dmixnorm_mclust <- function(x, mclust_solution, ...) {
  func_pi <- mclust_solution$parameters$pro
  func_mean <- mclust_solution$parameters$mean
  func_sigma <- mclust_solution$parameters$variance$sigma
  dmixnorm(x, func_pi, func_mean, func_sigma, ...)
}
dmixnorm_rmixmod <- function(x, mixmod_solution, ...) {
  K = 
  func_pi <- mixmod_solution@bestResult@parameters@proportions
  func_mean <- t(mixmod_solution@bestResult@parameters@mean)
  func_sigma <- array(do.call("c", mixmod_solution@bestResult@parameters@variance), c(NCOL(mixmod_solution@bestResult@parameters@mean),
                                                                                      NCOL(mixmod_solution@bestResult@parameters@mean),
                                                                                      NROW(mixmod_solution@bestResult@parameters@mean)))
  #func_sigma <- do.call("abind", list(mixmod_solution@bestResult@parameters@variance, along = 3))
  dmixnorm(x, func_pi, func_mean, func_sigma, ...)
}
## functions
dmixnorm_solution_func <- function(solution, ...) {
  if ("MixmodCluster" %in% methods::is(solution)) {
    if (solution@dataType == "quantitative") {
      return(function(x) dmixnorm_rmixmod(x, solution, ...))
    }
  }
  if ("Mclust" %in% methods::is(solution)) {
    return(function(x) dmixnorm_mclust(x, solution, ...))
  }
  stop("Not recognized format for solution")
}

get_order <- function(fittedMean, fittedVariance, mainMean, mainVariance) {
  K <- ncol(mainMean)
  res <- matrix(0, nrow = K, ncol = K)
  for (i in 1:K) {
    for (j in 1:K) {
      res[i, j] <- norm(mainVariance[, , j] - fittedVariance[, , i], type = "2") + norm(mainMean[, j] - fittedMean[, 
        i], type = "2")
    }
  }
  return(apply(res, 1, which.min))
}
#' CLR evaluated on gaussian mixture model posterioris of X
#' 
#' @param X dataframe where density function is evaluated
#' @param Pi a vector indicating the mixing proportions
#' @param Mu a two dimensional array where second component indicates the mean 
#' of each gaussian component
#' @param Sigma a three dimensional array where third component indicates the variance 
#' of each gaussian component
#' @export
clr_mixnorm <- function(X, Pi, Mu, Sigma) {
  log.dnorm <- lapply(1:length(Pi), function(i) log(Pi[i]) + mvtnorm::dmvnorm(X, mean = Mu[, i], sigma = Sigma[, 
    , i], log = T))
  log.dnorm.mean <- Reduce("+", log.dnorm)/length(log.dnorm)
  data.frame(do.call(cbind, lapply(log.dnorm, function(comp) comp - log.dnorm.mean)))
} 
