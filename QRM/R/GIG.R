## Copyright (C) 2013 Marius Hofert, Bernhard Pfaff
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


##
## Generalised Inverse Gaussian
##
## Random variates
rGIG <- function(n, lambda, chi, psi, envplot = FALSE, messages = FALSE){
   if((chi < 0) | (psi < 0))
     stop("Invalid parameters for GIG with chi, psi both negative")
   if((chi == 0) & (lambda <= 0))
     stop("Invalid parameters for GIG with chi=0 and lambda <=0")
   if((psi == 0) & (lambda >= 0))
     stop("Invalid parameters for GIG with psi=0 and lambda>=0")
   if((chi == 0) & (lambda > 0))
     return(rgamma(n, shape = lambda, rate = (psi/2)))
   if((psi == 0) & (lambda < 0))
     return(1/rgamma(n, shape = ( - lambda), rate = (chi/2)))
   message <- NULL
   if(abs(lambda) < 1)
     message <- paste(message, "Not necessarily efficient rejection method", "\n")
   neglambda <- FALSE
   if(lambda < 0){
     neglambda = TRUE
     lambda <- abs(lambda)
     tmp <- c(chi, psi)
     chi <- tmp[2]
     psi <- tmp[1]
   }
   #Local function:
   efunc <- function(x, lambda, chi, psi){
     (x^(lambda - 1)) * exp( - (chi/x + psi * x)/2)
   }
   #Local function:
   calcmode <- function(lambda, chi, psi){
     if(psi > 0){
       return(((lambda - 1) + sqrt(chi * psi + (1 - lambda)^2))/psi)
     } else if(psi == 0){
       return(chi/(2 - 2 * lambda))
     } else {
       stop("Problem in mode function")
     }
   }
   # Local function:
   objectiveOneParam <- function(theta,lambda, chi, psi, themode){
     if(lambda >= 1) stop("lambda must be less than 1 to call objectiveOneParam\n")
     if(length(theta)!= 1) stop("theta must be a single number in objectiveOneParam\n")
     if(theta <= 0){
       out <- NA
     } else {
       Delta1 <- (exp(themode * theta) - 1)/theta
       Delta2 <- (2 * exp(( - themode * psi)/2))/psi
       xL <- calcmode(lambda, chi, psi + 2 * theta)
       xH <- chi/(2 - 2 * lambda)
       S1 <- efunc(xL, lambda, chi, psi + 2 * theta)
       S2 <- efunc(xH, lambda, chi, 0)
       out <- Delta1 * S1 + Delta2 * S2
     }
     out
   }
   #Local function:
   objectiveTwoParams <- function(theta, lambda, chi, psi, themode){
     if(lambda < 1) stop("lambda must be at least equal to 1 to call objectiveTwoParam\n")
     if(length(theta)!= 2) stop("theta must be vector of length 2 in objectiveTwoParams\n")
     if((theta[1] <= 0) | (theta[2] <= 0)){
       out <- NA
     } else if((psi - 2 * theta[2]) < 0) {
       out <- NA
     } else {
       Delta1 <- (exp(themode * theta[1]) - 1)/theta[1]
       Delta2 <- exp( - themode * theta[2])/theta[2]
       xL <- calcmode(lambda, chi, psi + 2 * theta[1])
       xH <- calcmode(lambda, chi, psi - 2 * theta[2])
       S1 <- efunc(xL, lambda, chi, psi + 2 * theta[1])
       S2 <- efunc(xH, lambda, chi, psi - 2 * theta[2])
       out <- Delta1 * S1 + Delta2 * S2
     }
     out
  }
  #Call internal function:
  themode <- calcmode(lambda, chi, psi)
  if(lambda < 1)
  {
     lambdaLOW <- TRUE
     theta <- 0.01
     optimout <- optim(c(0.25), objectiveOneParam, gr = NULL, lambda = lambda, chi = chi, psi = psi, themode = themode, method = "BFGS")
     spar <- optimout$par[1]
     ppar <- psi/2
  } else {
     lambdaLOW <- FALSE
     theta <- c(0.01, psi/4)
     optimout <- optim(c(0.01, 0.25), objectiveTwoParams, gr = NULL, lambda = lambda, chi = chi, psi = psi, themode = themode, method = "BFGS")
     spar <- optimout$par[1]
     ppar <- optimout$par[2]
   }
   if(optimout$convergence != 0){
     message <- paste(optimout$message, "Problems finding optimal s and p (use option envplot for reassurance)","\n")
     print(message)
   }
   xL <- calcmode(lambda, chi, psi + 2 * spar)
   xH <- calcmode(lambda, chi, psi - 2 * ppar)
   S1 <- efunc(xL, lambda, chi, psi + 2 * spar)
   S2 <- efunc(xH, lambda, chi, psi - 2 * ppar)
   Delta1 <- (exp(themode * spar) - 1)/spar
   Delta2 <- exp( - themode * ppar)/ppar
   k <- 1/((Delta1/S2) + (Delta2/S1))
   k1 <- k/S2
   k2 <- k/S1
   rpar <- k1 * Delta1
   if(envplot){
     xdat <- seq(from = 0.01, to = themode * 20, length = 1000)
     envelope2 <- (xdat <= themode) * exp(spar * xdat) * S1 +
			(xdat > themode) * exp( - ppar * xdat) * S2
     envelope <- (xdat <= themode) * exp(spar * xdat) * k1 + (xdat >
			themode) * exp( - ppar * xdat) * k2
     ydat <- efunc(xdat, lambda, chi, psi)
     yr <- range(ydat, envelope, envelope2)
     plot(xdat, ydat, ylim = yr, type = "l")
     abline(v = themode)
     lines(xdat, envelope, col = 2)
     lines(xdat, envelope2, lty = 2, col = 2)
   }
   xsim <- rep(0, n)
   new = n
   tmp <- rgig(as.integer(n),
               as.double(rpar),
               as.double(spar),
               as.double(ppar),
               as.double(k1),
               as.double(k2),
               as.double(lambda),
               as.double(chi),
               as.double(psi),
               as.double(S1),
               as.double(S2))
   if(messages){
     efficiency <- n / new
     message <- paste(message, "Efficiency", round(efficiency * 100, 1),"\n")
     cat(message)
   }
   if(neglambda){
     return(1 / tmp)
   } else {
     return(tmp)
   }
}
## Fitting
EGIG <- function(lambda,chi,psi,k=1){
  if ((chi[1]>0) & (psi[1]>0)){
    term1 <- k*log(chi/psi)/2
    term2 <- log(besselK(x = sqrt(chi*psi), nu = lambda + k, expon.scaled = FALSE))
    term3 <- log(besselK(x = sqrt(chi * psi), nu = lambda, expon.scaled = FALSE))
    out <- exp(term1 + term2 - term3)
  } else if ((chi[1]==0) & (lambda>0)){
    alpha <- lambda
    beta <- psi/2
    out <- gamma(k + alpha) / (gamma(alpha) * beta^k)
  } else if ((psi[1]==0) & (lambda <0)){
    alpha <- -lambda
    beta <- chi / 2
    out <- (gamma(alpha - k) * beta^k) / gamma(alpha)
  } else {
    stop("These GIG parameters are not allowed")
  }
  out
}
## Fitting log-GIG
ElogGIG <- function(lambda, chi, psi){
 if ((chi[1] == 0) & (lambda > 0)){
    alpha <- lambda
    beta <- psi / 2
    out <- psi(alpha) - log(beta)
  } else if ((psi[1] == 0) & (lambda < 0)){
    alpha <- -lambda
    beta <- chi / 2
    out <- log(beta) - psi(alpha)
  } else {
    stop("Log Moment of general GIG not implemented")
  }
  out
}
