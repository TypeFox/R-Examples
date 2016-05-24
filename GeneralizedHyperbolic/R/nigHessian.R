### Calculate the hessian of the Normal Inverse Gaussian distribution for
### different parameterizations
### CYD 23/04/10
### CYD 22/04/10
### CYD 01/04/10
### 1:4 param = hyperbChangePars 1:4;
### 5 = log scale of number 2 and 4 elements of param
nigHessian <- function(x, param, hessianMethod = c("tsHessian", "exact"),
                          whichParam = 1:5, ...) {
  if (hessianMethod == "exact") {
    stop("Exact hessian not implemented yet. Use method tsHessian instead.")
  }

  else {
    if (whichParam == 1) {
      llfuncH <- function(param) {
        llparam <- param
        mu <- llparam[1]
        delta <- llparam[2]
        hyperbPi <- llparam[3]
        zeta <- llparam[4]
        KNu <- besselK(zeta, nu = -1/2)
        KNu2 <- besselK((zeta*sqrt(1 + hyperbPi^2)/delta)*
                sqrt(delta^2 + (x - mu)^2), nu = -1)
        nigDens <- sqrt(1 + hyperbPi^2)*zeta^(1/2)*
                   (sqrt(2*pi)*KNu)^(-1)*(delta^2 + (x - mu)^2)^(-1/2)*
                   KNu2*exp((zeta*hyperbPi/delta)*(x - mu))
        return(sum(log(nigDens)))
      }
    } else if (whichParam == 3) {
      llfuncH <- function(param) {
        llparam <- hyperbChangePars(3, 2, param = param)
        return(sum(log(dnig(x = x, param = llparam))))
      }
    } else if (whichParam == 4) {
      llfuncH <- function(param) {
        llparam <- hyperbChangePars(4, 2, param = param)
        return(sum(log(dnig(x = x, param = llparam))))
      }
    } else if (whichParam == 2) {
      llfuncH <- function(param) {
        llparam <- param
        return(sum(log(dnig(x = x, param = llparam))))
      }
      
    } else if (whichParam == 5) {
      llfuncH <- function(param) {
        mu <- param[1]
        delta <- exp(param[2])
        hyperbPi <- param[3]
        zeta <- exp(param[4])
        KNu <- besselK(zeta, nu = -1/2)
        KNu2 <- besselK((zeta*sqrt(1 + hyperbPi^2)/delta)*
                        sqrt(delta^2 + (x - mu)^2), nu = -1)
        nigDens <- sqrt(1 + hyperbPi^2)*zeta^(1/2)*
                   (sqrt(2*pi)*KNu)^(-1)*(delta^2 + (x - mu)^2)^(-1/2)*
                     KNu2*exp((zeta*hyperbPi/delta)*(x - mu))
        return(sum(log(nigDens)))
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


################################################################################

