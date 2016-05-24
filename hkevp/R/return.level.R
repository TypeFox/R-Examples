#' @title The associated return level
#' @author Quentin Sebille
#'
#'
#' @description 
#' Computation of the associated return level with given period and GEV parameters.
#' 
#' 
#' 
#' @param period
#' An integer greater than 1.
#' Gives the period \eqn{T} over which the return level is computed. See details.
#' 
#' @param mu
#' A numerical value or a vector of real values.
#' GEV location parameter. Must be of length one or same length as \code{sigma} and/or \code{xi}.
#' @param sigma
#' A numerical value or a vector of real values.
#' GEV scale parameter. Must be of length one or same length as \code{mu} and/or \code{xi}.
#' @param xi
#' A numerical value or a vector of real values.
#' GEV shape parameter. Must be of length one or same length as \code{mu} and/or \code{sigma}.
#'
#'
#'
#' @details
#' The \eqn{T}-year return level is a typical value of interest in Extreme Value Theory. It represents the value that is expected to be exceeded once over \eqn{T} years. Given the parameters \eqn{\mu}, \eqn{\sigma} and \eqn{\xi} of the GEV distribution associated to the yearly maxima, we can compute the associated \eqn{T}-return level \eqn{y_T} by:
#' \deqn{y_T := \mu + \frac{\sigma}{\xi} \left[ \log\left(\frac{T}{T-1}\right)^{-\xi} -1 \right] ~.}
#'
#'
#'
#' @return 
#' This function returns either a numerical value or a numerical vector. If \code{mu}, \code{sigma}, \code{xi} are of length one each, it returns a single numerical value. Otherwise, it returns a vector with length corresponding to the given GEV parameters.
#'
#'
#' @examples
#' return.level(period = 100, mu = 1, sigma = 1, xi = 1)
#' return.level(period = 200, mu = 1:10, sigma = 1, xi = 0)
#' 
#' 
return.level <- function(period, mu, sigma, xi) {
  ## Errors
  if ((period != round(period)) || (period <= 1)) stop("'period' must be an integer greater than 1!")
  LENGTHS <- c(length(mu), length(sigma), length(xi))
  if ( !(length(unique(LENGTHS[LENGTHS != 1])) %in% c(0,1)) ) stop("GEV parameters must have same length (or length one)!")
  
  ## Initialisation of the result table
  mu <- as.vector(mu)
  sigma <- as.vector(sigma)
  xi <- as.vector(xi)
  gev <- cbind(mu, sigma, xi)
  result <- rep(NA, nrow(gev))
  
  ## Separation of two cases : xi == 0 and xi != 0
  isxi0 <- gev[,3] == 0
  result[isxi0] <- gev[isxi0,1] - gev[isxi0,2]*log(log(period/(period - 1)))
  result[!isxi0] <- gev[!isxi0,1] + gev[!isxi0,2]*(log(period/(period - 1)) ^ (-gev[!isxi0,3]) - 1)/gev[!isxi0,3]
  
  return(result)
}
