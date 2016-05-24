### Calculate the hessian of the hyperbolic distribution for
### different parameterizations
### CYD 23/04/10
### CYD 22/04/10
### CYD 01/04/10
### 1:4 param = hyperbChangePars 1:4;
### 5 = log scale of number 2 and 4 elements of param
hyperbHessian <- function(x, param, hessianMethod = c("exact", "tsHessian"),
                          whichParam = 1:5) {
  if (hessianMethod == "exact") {
    hessian <- matrix(nrow = 4, ncol = 4)
    n <- length(x)
    if (whichParam == 1) {
      mu <- param[1]
      delta <- param[2]
      hyperbPi <- param[3]
      zeta <- param[4]

      hessian[1,1] <- -zeta*sqrt(1 + hyperbPi^2)*delta*
                      sumX(x, mu, delta, r = 0, k = 3)
      hessian[1,2] <- hessian[2,1] <-
                      -zeta/delta^2*(sqrt(1 + hyperbPi^2)*
                      sumX(x, mu, delta, r = 3, k = 3) -
                      2*sqrt(1 + hyperbPi^2)*sumX(x, mu, delta, r = 1, k = 1)
                      - n*hyperbPi)
      hessian[1,3] <- hessian[3,1] <-
                      -zeta/delta*(hyperbPi/sqrt(1 + hyperbPi^2)*
                      sumX(x, mu, delta, r = 1, k = 1) + n)
      hessian[1,4] <- hessian[4,1] <-
                      -1/delta*(sqrt(1 + hyperbPi^2)*
                      sumX(x, mu, delta, r = 1, k = 1) + n*hyperbPi)
      hessian[2,2] <- n/delta^2 + zeta/delta^3*(sqrt(1 + hyperbPi^2)*
                      sumX(x, mu, delta, r = 4, k = 3) -
                      3*sqrt(1 + hyperbPi^2)*sumX(x, mu, delta, r = 2, k = 1) -
                      2*hyperbPi*sumX(x, mu, delta, r = 1, k = 0))
      hessian[2,3] <- hessian[3,2] <-
                      zeta/delta^2*(hyperbPi/sqrt(1 + hyperbPi^2)*
                      sumX(x, mu, delta, r = 2, k = 1) +
                      sumX(x, mu, delta, r = 1, k = 0))
      hessian[2,4] <- hessian[4,2] <-
                      1/delta^2*(sqrt(1 + hyperbPi^2)*
                      sumX(x, mu, delta, r = 2, k = 1) +
                      hyperbPi*sumX(x, mu, delta, r = 1, k = 0))
      hessian[3,3] <- n/(1 + hyperbPi^2)^2*(hyperbPi^2 - 1) -
                      zeta/((1 + hyperbPi^2)^(7/2)*delta)*
                      sumX(x, mu, delta, r = 0, k = -1)*
                      (2*hyperbPi^2 + hyperbPi^4 + 1)
      hessian[3,4] <- hessian[4,3] <-
                      -1/delta*(sumX(x, mu, delta, r = 1, k = 0) +
                      (hyperbPi/sqrt(1 + hyperbPi ^2))*
                      sumX(x, mu, delta, r = 0, k = -1))
      hessian[4,4] <- -n*(1 + 1/zeta^2 -
                      besselRatio(x = zeta, nu = 1, orderDiff = -1)/zeta -
                      (besselRatio(x = zeta, nu = 1, orderDiff = -1))^2)

    } else if (whichParam == 2) {
      mu <- param[1]
      delta <- param[2]
      alpha <- param[3]
      beta <- param[4]
      tau <- delta*sqrt(alpha^2 - beta^2)
      brTau <- besselRatio(x = tau, nu = 1, orderDiff = -1)
      hessian[1,1] <- -alpha*delta^2*sumX(x, mu, delta, r = 0, k = 3)
      hessian[1,2] <- hessian[2,1] <-
                      alpha*delta*sumX(x, mu, delta, r = 1, k = 3)
      hessian[1,3] <- hessian[3,1] <-
                      -sumX(x, mu, delta, r = 1, k = 1)
      hessian[1,4] <- hessian[4,1] <- -n
      hessian[2,2] <- n*(tau/delta^2*brTau +
                      (alpha^2 - beta^2)*(brTau^2 - 1)) +
                      alpha*delta^2*sumX(x, mu, delta, r = 0, k = 3) -
                      alpha*sumX(x, mu, delta, r = 0, k = 1)
      hessian[2,3] <- hessian[3,2] <-
                      -n*delta*alpha + 2*alpha*n*delta/tau*brTau +
                      n*delta*alpha*brTau^2 -
                      delta*sumX(x, mu, delta, r = 0, k = 1)
      hessian[2,4] <- hessian[4,2] <-
                      n*(delta*beta - 2*beta*delta/tau*brTau -
                      beta*delta*brTau^2)
      hessian[3,3] <- (n*(-alpha^6*delta^2 + (-1 + delta^2*beta^2)*
                      alpha^4 - 4*alpha^2*beta^2 + beta^4))/
                      (alpha^2*(alpha^2 - beta^2)^2) +
                      n*delta^2/tau*brTau +
                      n*delta^2*alpha^2/(alpha^2 - beta^2)*brTau^2
      hessian[3,4] <- hessian[4,3] <-
                      n*beta*alpha*(4 + alpha^2*delta^2 - beta^2*delta^2)/
                      (alpha^2 - beta^2)^2 -
                      n*delta^2*beta*alpha/(alpha^2 - beta^2)*brTau^2
      hessian[4,4] <- (-2*n*(alpha^2 + beta^2) - n*beta^2*delta^2*
                      (alpha^2 - beta^2))/(alpha^2 - beta^2)^2 -
                      n*delta^2/tau*brTau +
                      n*delta^2*beta^2/(alpha^2 - beta^2)*brTau^2

    } else if (whichParam == 3 | whichParam == 4) {
      stop("The exact hessian formula is not available for this parameterization yet. Please use method tsHessian instead.")
    } else if (whichParam == 5) {
      mu <- param[1]
      delta <- exp(param[2])
      hyperbPi <- param[3]
      zeta <- exp(param[4])

      hessian[1,1] <- -zeta*sqrt(1 + hyperbPi^2)*delta*
                      sumX(x, mu, delta, r = 0, k = 3)
      hessian[1,2] <- hessian[2,1] <-
                      -zeta/delta*(sqrt(1 + hyperbPi^2)*
                      sumX(x, mu, delta, r = 3, k = 3) -
                      2*sqrt(1 + hyperbPi^2)*sumX(x, mu, delta, r = 1, k = 1)
                      - n*hyperbPi)
      hessian[1,3] <- hessian[3,1] <-
                      -zeta/delta*(hyperbPi/sqrt(1 + hyperbPi^2)*
                      sumX(x, mu, delta, r = 1, k = 1) + n)
      hessian[1,4] <- hessian[4,1] <-
                      -zeta/delta*(sqrt(1 + hyperbPi^2)*
                      sumX(x, mu, delta, r = 1, k = 1) + n*hyperbPi)
      hessian[2,2] <- zeta/delta*(sqrt(1 + hyperbPi^2)*
                      sumX(x, mu, delta, r = 4, k = 3) -
                      2*sqrt(1 + hyperbPi^2)*sumX(x, mu, delta, r = 2, k = 1) -
                      hyperbPi*sumX(x, mu, delta, r = 1, k = 0))
      hessian[2,3] <- hessian[3,2] <-
                      zeta/delta*(hyperbPi/sqrt(1 + hyperbPi^2)*
                      sumX(x, mu, delta, r = 2, k = 1) +
                      sumX(x, mu, delta, r = 1, k = 0))
      hessian[2,4] <- hessian[4,2] <-
                      zeta/delta*(sqrt(1 + hyperbPi^2)*
                      sumX(x, mu, delta, r = 2, k = 1) +
                      hyperbPi*sumX(x, mu, delta, r = 1, k = 0))
      hessian[3,3] <- n/(1 + hyperbPi^2)^2*(hyperbPi^2 - 1) -
                      zeta/((1 + hyperbPi^2)^(7/2)*delta)*
                      sumX(x, mu, delta, r = 0, k = -1)*
                      (2*hyperbPi^2 + hyperbPi^4 + 1)
      hessian[3,4] <- hessian[4,3] <-
                      -zeta/(sqrt(1 + hyperbPi^2)*delta)*
                      (sqrt(1 + hyperbPi^2)*sumX(x, mu, delta, r = 1, k = 0) +
                      hyperbPi*sumX(x, mu, delta, r = 0, k = -1))
      hessian[4,4] <- -n*zeta^2 +
                      2*n*zeta*besselRatio(x = zeta, nu = 1, orderDiff = -1) +
                      n*zeta^2*
                      (besselRatio(x = zeta, nu = 1, orderDiff = -1))^2 -
                      sqrt(1 + hyperbPi^2)*zeta/delta*
                      sumX(x, mu, delta, r = 0, k = -1) -
                      hyperbPi*zeta/delta*sumX(x, mu, delta, r = 1, k = 0)
    }

  } else {
    if (whichParam == 1) {
      llfuncH <- function(param) {
        llparam <- param
        KNu <- besselK(llparam[4], nu = 1)
        hyperbDens <- (2*llparam[2]* sqrt(1 + param[3]^2)*KNu)^(-1)*
                      exp(-llparam[4]* (sqrt(1 + param[3]^2)*
                      sqrt(1 + ((x - param[1])/llparam[2])^2) -
                      param[3]*(x - param[1])/llparam[2]))
        return(sum(log(hyperbDens)))
      }
    } else if (whichParam == 3) {
      llfuncH <- function(param) {
        llparam <- hyperbChangePars(3, 2, param = param)
        return(sum(log(dhyperb(x = x, param = llparam))))
      }
    } else if (whichParam == 4) {
              llfuncH <- function(param) {
                llparam <- hyperbChangePars(4, 2, param = param)
                return(sum(log(dhyperb(x = x, param = llparam))))
              }
            } else if (whichParam == 2) {
              llfuncH <- function(param) {
                llparam <- param
                return(sum(log(dhyperb(x = x, param = llparam))))
              }

            } else if (whichParam == 5) {
              llfuncH <- function(param) {
                KNu <- besselK(exp(param[4]), nu = 1)
                hyperbDens <- (2*exp(param[2])* sqrt(1 + param[3]^2)*KNu)^(-1)*
                              exp(-exp(param[4])* (sqrt(1 + param[3]^2)*
                              sqrt(1 + ((x - param[1])/exp(param[2]))^2) -
                              param[3]*(x - param[1])/exp(param[2])))
                return(sum(log(hyperbDens)))
              }
            }

    hessian <- tsHessian(param = param, fun = llfuncH)
  }
  if (whichParam == 1) {
    rownames(hessian) <- colnames(hessian) <- c("mu","delta","hyperbPi","zeta")
  } else if (whichParam == 2) {
    rownames(hessian) <- colnames(hessian) <- c("mu","delta","alpha","beta")
  } else if (whichParam == 3) {
    rownames(hessian) <- colnames(hessian) <- c("mu","delta","phi","gamma")
  } else if (whichParam == 4) {
    rownames(hessian) <- colnames(hessian) <- c("mu","delta","xi","chi")
  } else if (whichParam == 5) {
    rownames(hessian) <- colnames(hessian) <-
      c("mu","log(delta)","hyperbPi","log(zeta)")
  }

  return(hessian)
}




sumX <- function(x, mu, delta, r, k) {
  sum((mu - x)^r / (delta^2 + (x - mu)^2)^(k/2))
}


################################################################################

