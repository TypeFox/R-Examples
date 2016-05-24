#internal function representing the probability that an event will be occuring x seconds into the future, given that no event is currently occurring

p_0 <- function(x,phi,zeta) phi * (1- exp(-1 * x * zeta / (phi * (1- phi))))


#internal function representing the probability that an event will be occuring x seconds into the future, given that an event is currently occurring

p_1 <- function(x,phi,zeta) phi + (1 - phi) * exp(-1 * x * zeta / (phi * (1- phi)))


#calculates a vector of probabilities that an event is occuring at each U time point, given the interval record up to that time

PIRpsi <- function(phi, zeta, U, c, d) {
  psi <- vector(mode = "numeric", length = length(U))
  psi[1] <- phi
  for(i in 2:length(U)){
    psi[i] <- ((psi[i-1] * p_1(c + d,phi,zeta) + (1 - psi[i-1]) * (p_0(c + d,phi,zeta) - p_0(d, phi, zeta) * exp(-1 * zeta * c / (1-phi)))) / (1 - (1-psi[i-1]) * exp(-1 * zeta * c/(1-phi))))^U[i-1] * p_0(d,phi,zeta)^(1-U[i-1]) 
  }
  return(psi)
}


#' @title Calculate log-likelihood
#' 
#' @description Calculates the log-likelihood of within-session PIR data
#' 
#' @param phi  value for prevalence
#' @param zeta value for incidence
#' @param U a vector containing interval-level PIR data
#' @param c the length of the active interval
#' @param d the length of the recording interval
#' 
#' @details \code{phi} must be a value between 0 and 1, inclusive of 0 and 1. \
#' code{zeta} must be a value greater than 0.
#' 
#' The vector \code{U} should only contain values of 1 or 0.
#' 
#' \code{c} must be some positive value in whatever time units the observation took place in (typically seconds). \code{d} must be some non-negative value - a 
#' \code{d} of zero represents a PIR observation where no time was set aside for recording.
#' 
#' @return The value of the log-likelihood
#'
#' 
#' @author Daniel Swan <dswan@@utexas.edu>
#' @export

PIR_loglik <- function(phi, zeta, U, c, d){
  
  if (phi < 0 | phi > 1) stop("Values for phi must be between 0 and 1, inclusive")
  if (zeta < 0) stop("Values for zeta must be greater than 0")
  if (length(U[which(U != 0 & U != 1)]) > 0) stop("Values for U must be 0 or 1")
  if (c <= 0) stop("c must be a positive value")
  if (d < 0)  stop("d cannot be negative")
  
  psi <- PIRpsi(phi, zeta, U, c, d)
  loglik <- sum(U * log(1 - (1 - psi) * exp(-1 * zeta * c / (1-phi))) + 
                  (1 - U) * (log(1 - psi) - zeta * c / (1-phi)))
  return(loglik)
}


