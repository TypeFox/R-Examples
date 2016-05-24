ssbayes <- function(k, nu, q0, deltastar, eta, zeta, alpha0, beta0, xi, prec="known", crit="1"){
  
  prec <- match.arg(prec, choices=c("known", "unknown"))
  crit <- match.arg(crit, choices=c("1", "2"))
  
  if(length(q0)!=(k + 1)){stop("q0 must have k+1 elements.")}
  
  if(prec=="known"){
    
    if(crit=="1"){
      
      rho <- 1 / (1 + sqrt(k))
      Sigma <- matrix(c(rep(rho, k * k)), k, k, byrow=TRUE) + diag(1 - rho, k)
      V1 <- ((qnorm(eta, 0, 1) + qmvnorm(zeta, sigma=Sigma)$quantile) / deltastar)^2
      nopt <- (1 + 1/sqrt(k)) * V1/nu - q0
      nopt[1] <- (1 + sqrt(k)) * V1/nu - q0[1]
      
    }else{
      
      V2 <- ((qnorm(eta, 0, 1) + qnorm(zeta, 0, 1)) / deltastar)^2
      nopt <- (1 + 1/sqrt(k)) * V2/nu - q0
      nopt[1] <- (1 + sqrt(k)) * V2/nu - q0[1]
      
    }    
    
  }else{
    
    if(crit=="1"){
      
      Vn1 <- function(n, k, q0, eta, zeta, deltastar, alpha0, beta0, xi){
        rho <- 1 / (1 + sqrt(k))
        Sigma <- matrix(c(rep(rho, k * k)), k, k, byrow=TRUE) + diag(1 - rho, k)
        nint <- ceiling(n)
        alpha1 <- alpha0 + nint/2
        beta1 <- beta0 / (1 - qbeta(xi, nint/2, alpha0))
        Vn <- ((deltastar / (qmvt(zeta, tail="lower.tail", df=(2 * alpha1), sigma=Sigma)$quantile +
                               qt(eta, df=(2 * alpha1))))^2 * alpha1/beta1)^(-1)
        return(Vn)
      }
      
      fun1 <- function(n, k, q0, eta, zeta, deltastar, alpha0, beta0, xi){
        diff <- n - (sqrt(k) + 1)^2 * Vn1(n, k, q0, eta, zeta, deltastar, alpha0, beta0, xi) + sum(q0)
        return(diff)
      }
      
      optn <- uniroot(fun1, c(1, 10000), k, q0, eta, zeta, deltastar, alpha0, beta0, xi)$root
      nopt <- (1 + 1/sqrt(k)) * Vn1(optn, k, q0, eta, zeta, deltastar, alpha0, beta0, xi) - q0
      nopt[1] <- (1 + sqrt(k)) * Vn1(optn, k, q0, eta, zeta, deltastar, alpha0, beta0, xi) - q0[1]
      
    }else{
      
      Vn2 <- function(n, k, q0, eta, zeta, deltastar, alpha0, beta0, xi){
        nint <- ceiling(n)
        alpha1 <- alpha0 + nint/2
        beta1 <- beta0 / (1 - qbeta(xi, nint/2, alpha0))
        Vn <- ((deltastar / (qt(zeta, df=(2 * alpha1)) + qt(eta, df=(2 * alpha1))))^2 * alpha1/beta1)^(-1)
        return(Vn)
      }
      
      fun2 <- function(n, k, q0, eta, zeta, deltastar, alpha0, beta0, xi){
        diff <- n - (sqrt(k) + 1)^2 * Vn2(n, k, q0, eta, zeta, deltastar, alpha0, beta0, xi) + sum(q0)
        return(diff)
      }
      
      optn <- uniroot(fun2, c(1, 10000), k, q0, eta, zeta, deltastar, alpha0, beta0, xi)$root
      nopt <- (1 + 1/sqrt(k)) * Vn2(optn, k, q0, eta, zeta, deltastar, alpha0, beta0, xi) - q0
      nopt[1] <- (1 + sqrt(k)) * Vn2(optn, k, q0, eta, zeta, deltastar, alpha0, beta0, xi) - q0[1]
      
    }
    
  }

  if(min(nopt) < 0){warning("There are sample sizes calculated as <0 and therefore set to 0.")}
  
  nopt <- as.matrix(sapply(nopt, function(x) max(0, ceiling(x))))
  rownames(nopt) <- c("Control:", paste("Group ", LETTERS[1:k], ":", sep=""))
  colnames(nopt) <- ""
  return(nopt)
  
}
