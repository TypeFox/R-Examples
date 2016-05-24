marginBitsD <- function(params, eta1, eta2, X3, X4, dat, VC, precision, univariate=TRUE, bd=bd) {
    
    sigma <- nu <- df2.sigma <- df2.nu <- dF2.sigma <- dF22.sigma <- dF2.nu <- dF22.nu <- d2f2sigma2 <- NULL
    d2f2delta2sigma <- d2f2nu2 <- d2f2delta2nu <- d2f2nusigma <- d2f2sigmanu <- d2F2dsigma2 <- NULL 
    d2F22dsigma2 <- d2F2ddelta2dsigma <- d2F22ddelta2dsigma <- d2F2dnu2 <- d2F22dnu2 <- d2F2ddelta2dnu <- NULL 
    d2F22ddelta2dnu <- d2F2dsigmadnu <- d2F22dsigmadnu <- d2f2nusigma <- d2F2dnudsigma <- d2F22dnudsigma <- NULL
    F1 <- F2 <- F22 <- dF1 <- d2F1ddelta1delta1 <- dF2 <- dF22 <- dF2.sigma <- dF22.sigma <- dF2.nu <- NULL
    dF22.nu <- d2F2ddelta22 <- d2F22ddelta22 <- d2F2dsigma2 <- d2F22dsigma2 <- d2F2ddelta2dsigma <- NULL 
    d2F22ddelta2dsigma <- d2F2dnu2 <- d2F22dnu2 <- d2F2ddelta2dnu <- d2F22ddelta2dnu <- d2F2dnudsigma <- NULL
    d2F22dnudsigma <- d2F2dsigmadnu <- d2F22dsigmadnu <- NULL


  etasqv <- etanu <- NULL
  
  eps <- sqrt(.Machine$double.eps)


  i1 <- dat[,1]
  i0 <- 1-i1
  ind <- i1==0 

  if(univariate==FALSE) {
  
    F1 <- pnorm(-eta1)        
    F1 <- ifelse(F1>0.00000001,F1,0.00000001)
    F1 <- ifelse(F1<0.99999999,F1,0.99999999)
    dF1 <- -dnorm(-eta1)
    d2F1ddelta1delta1 <- -dF1*eta1
  
  }
  
  #fd.prec <- 10^(-7) 
  fd.prec <- eps
  
 if(VC$margins[2]=="P") {
    i2 <- dat[,2]
    mu <- exp(eta2) 
    f2 <- dpois(i2, mu)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- ppois(i2, mu)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
   
   
    # df2, dF2, dF22  
    gamma <- exp(-exp(eta2))*exp(eta2)^(i2)/(dpois(i2, exp(eta2)))
    df2 <- as.vector(1/gamma*(-exp(-exp(eta2))*exp(eta2)^(i2)+exp(-exp(eta2))*i2*exp(eta2)^(i2-1))*exp(eta2))
    
   if(univariate==FALSE) {
  
      df2.int.P <- function(y, mu) { 
        gamma <- exp(-mu)*mu^(y)/(dpois(y, mu))
        as.vector(1/gamma*(-exp(-mu)*mu^(y)+exp(-mu)*y*mu^(y-1))*mu) 
      }
    
      ly <- length(i2)
      cdf.f <- rep(0, ly)
      nmu <- rep(exp(eta2), length = ly) 
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        allval <- seq(0, y.y)
        pdfall <- df2.int.P(allval, mu = mm)
        cdf.f[i] <- sum(pdfall)
      }  
    
      dF2 <- cdf.f
      dF22 <- dF2-df2
    
   }
    
    # Hessian derivative components
    
    gamma.plus <- exp(-exp(eta2+fd.prec))*exp(eta2+fd.prec)^(i2)/(dpois(i2, exp(eta2+fd.prec)))
    df2.plus <- as.vector(1/gamma.plus*(-exp(-exp(eta2+fd.prec))*exp(eta2+fd.prec)^(i2)+exp(-exp(eta2+fd.prec))*i2*exp(eta2+fd.prec)^(i2-1))*exp(eta2+fd.prec))
    d2f2delta22 <- (df2.plus - df2)/fd.prec
    
    if(univariate==FALSE) {
   
      ly <- length(i2)
      cdf.f.plus <- rep(0, ly)
      nmu.plus <- rep(exp(eta2 + fd.prec), length = ly) 
    
      for (i in 1:ly) {
      
        y.y <- i2[i]
        mm <- nmu.plus[i]
        allval <- seq(0, y.y)
        pdfall <- df2.int.P(allval, mu = mm)
        cdf.f.plus[i] <- sum(pdfall)
      
      }  
    
      dF2.plus <- cdf.f.plus
      d2F2ddelta22 <- (dF2.plus - dF2) / fd.prec
    
    
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    }
    
    
    
    
  } else if (VC$margins[2]=="NB") {
    i2 <- dat[,2]
    
    if(univariate==TRUE) {
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
    }
      
    sigma <- exp(sigma.star)
    sigma <- ifelse(sigma>10^-8, sigma, 10^-8)
    sigma <- ifelse(sigma<10^7, sigma, 10^7)
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^7, mu, 10^7)
    f2 <- dNBI(i2, sigma = sigma,  mu = mu)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- pNBI(i2, sigma = sigma, mu = mu)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
    
    df2 <- f2*(i2*(mu*sigma)^(-1)*sigma*mu-(i2+1/sigma)*(mu*sigma+1)^(-1)*sigma*mu)
    df2 <- as.vector(df2)
    
    df2.sigma <- f2*(digamma(i2+1/sigma)*(-1/sigma^2)+i2*(mu*sigma)^(-1)*mu-(digamma(1/sigma)*(-1/sigma^2)+(-1/sigma^2)*log(mu*sigma+1)+(i2+1/sigma)*(1/(mu*sigma+1))*mu))
    
    
    if(univariate==FALSE) {
    
      df2.mu.int.NB<-function(y, mu, sigma, nu) {
        f2 <- dNBI(y, sigma = sigma,  mu = mu)
        f2 <- ifelse(f2>precision, f2, precision)
        df2 <- f2*(y*(mu*sigma)^(-1)*sigma*mu-(y+1/sigma)*(mu*sigma+1)^(-1)*sigma*mu) 
       df2
      }
    
      df2.sigma.int.NB<-function(y, mu, sigma, nu) {
        f2 <- dNBI(y, sigma = sigma,  mu = mu)
        f2 <- ifelse(f2>precision, f2, precision)
        df2 <- f2*(digamma(y+1/sigma)*(-1/sigma^2)+y*(mu*sigma)^(-1)*mu-(digamma(1/sigma)*(-1/sigma^2)+(-1/sigma^2)*log(mu*sigma+1)+(y+1/sigma)*(1/(mu*sigma+1))*mu)) 
        df2
      }
    
    
      ly <- length(i2)
      cdf.f <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        pdfall <- df2.mu.int.NB(allval, mm, ms)
        cdf.f[i] <- sum(pdfall)
      }  
    
      dF2 <- cdf.f
      dF22 <- dF2 - df2
    
    
      ly <- length(i2)
      cdf.f.sigma <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        pdfall <- df2.sigma.int.NB(allval, mm, ms)
        cdf.f.sigma[i] <- sum(pdfall)
     } 
    
      dF2.sigma <- cdf.f.sigma
      dF22.sigma <- dF2.sigma - df2.sigma
    
    }
    
    # Hessian derivative components
    
    
    mu.plus <- ifelse(exp(eta2+fd.prec)>0.0001, exp(eta2+fd.prec), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    
    f2.plus <- dNBI(i2, sigma = sigma,  mu = mu.plus)
    f2.plus <- ifelse(f2.plus>precision, f2.plus, precision)
    df2.plus <- f2.plus*(i2*(mu.plus*sigma)^(-1)*sigma*mu.plus-(i2+1/sigma)*(mu.plus*sigma+1)^(-1)*sigma*mu.plus)
    df2.plus <- as.vector(df2.plus)
    d2f2delta22 <- (df2.plus - df2) / fd.prec
    
    
    sigma.plus <- sigma+fd.prec
    f2.plus.sigma <- dNBI(i2, sigma = sigma.plus,  mu = mu)
    f2.plus.sigma <- ifelse(f2.plus.sigma>precision, f2.plus.sigma, precision)
    df2.sigma.plus <- f2.plus.sigma*(digamma(i2+1/sigma.plus)*(-1/sigma.plus^2)+i2*(mu*sigma.plus)^(-1)*mu-(digamma(1/sigma.plus)*(-1/sigma.plus^2)+(-1/sigma.plus^2)*log(mu*sigma.plus+1)+(i2+1/sigma.plus)*(1/(mu*sigma.plus+1))*mu))
    d2f2sigma2 <- (df2.sigma.plus - df2.sigma) / fd.prec
    
    
    df2.dsigma.plus <- f2.plus.sigma*(i2*(mu*sigma.plus)^(-1)*sigma.plus*mu-(i2+1/sigma.plus)*(mu*sigma.plus+1)^(-1)*sigma.plus*mu)
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec
    
    if(univariate==FALSE) {
    
      ly <- length(i2)
      cdf.f.plus <- rep(0, ly)
      nmu <- rep(mu.plus, length = ly) 
      nsigma <- rep(sigma, length = ly)
    
      for (i in 1:ly) {
         y.y <- i2[i]
         mm <- nmu[i]
         ms <- nsigma[i]
         allval <- seq(0, y.y)
         pdfall <- df2.mu.int.NB(allval, mm, ms)
         cdf.f.plus[i] <- sum(pdfall)
        }  
      
      dF2.mu.plus <- cdf.f.plus
      d2F2ddelta22 <- (dF2.mu.plus - dF2) / fd.prec
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    
      ly <- length(i2)
      cdf.f.sigma.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma.plus, length = ly)
    
      for (i in 1:ly) {
         y.y <- i2[i]
         mm <- nmu[i]
         ms <- nsigma[i]
         allval <- seq(0, y.y)
         pdfall <- df2.sigma.int.NB(allval, mm, ms)
         cdf.f.sigma.plus[i] <- sum(pdfall)
        } 
    
      dF2.sigma.plus <- cdf.f.sigma.plus
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
    
      ly <- length(i2)
      cdf.f.dsigma.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma.plus, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        pdfall <- df2.mu.int.NB(allval, mm, ms)
        cdf.f.dsigma.plus[i] <- sum(pdfall)
      } 
    
      dF2.dsigma.plus <- cdf.f.dsigma.plus
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma   
    
    }
    
  } else if (VC$margins[2]=="D") {
    i2 <- dat[,2]
    
    if(univariate==TRUE) {
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
      if (is.null(VC$X3)) nu.star <- params[(VC$X2.d2 + 2)]
      if (!is.null(VC$X3)) nu.star <- etanu <- X4 %*% params[(VC$X2.d2 + VC$X3.d2 + 1):(VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]        
      if(is.null(VC$X3))  nu.star <- params[(VC$X1.d2+VC$X2.d2+2)]
      if(!is.null(VC$X3)) nu.star <- etanu <- X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)] 
    }
      
    nu <- plogis(nu.star)
    nu <- ifelse(nu<precision, precision, nu)
    nu <- ifelse(nu>(1-precision), 1-precision, nu)
    sigma <- exp(sigma.star)
    sigma <- ifelse(sigma>precision, sigma, precision)
    sigma <- ifelse(sigma<10^7, sigma, 10^7)
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^7, mu, 10^7)  
    f2 <- dDEL(i2, mu=mu, sigma=sigma, nu=nu)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- pDEL(i2, mu=mu, sigma=sigma, nu=nu)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
    
    
    
    bigS <- getFromNamespace("tofyDEL2", "gamlss.dist")
    
    
    S_prime_mu <- (exp(bigS(i2, mu+fd.prec, sigma, nu))-exp(bigS(i2, mu, sigma, nu)))/fd.prec
    S_prime_eta2 <- mu*S_prime_mu
    S <- exp(bigS(i2, mu, sigma, nu))
    df2 <- as.vector(f2*(-mu*nu+(-1/sigma)*(1+mu*sigma*(1-nu))^(-1)*mu*sigma*(1-nu)+S_prime_eta2/S))
    
    S_prime_sigma <- (exp(bigS(i2, mu, sigma+fd.prec, nu))-exp(bigS(i2, mu, sigma, nu)))/fd.prec
    d_prime <- ((1/sigma^2)*log(1+mu*sigma*(1-nu))-(1/sigma)*(1/(1+mu*sigma*(1-nu)))*mu*(1-nu))*(1+mu*sigma*(1-nu))^(-1/sigma)
    df2.sigma <- f2*(d_prime*(1+mu*sigma*(1-nu))^(1/sigma)+S_prime_sigma/S)
    
    S_prime_nu <- (exp(bigS(i2, mu, sigma, nu+fd.prec))-exp(bigS(i2, mu, sigma, nu)))/fd.prec
    df2.nu <- f2*(-mu+((-1/sigma)*(1+mu*sigma*(1-nu))^(-1)*(-mu*sigma)+S_prime_nu/S))
    
    
    if(univariate==FALSE) {
    
      df2.mu.int.D <- function(y, mu, sigma, nu) {
        f2 <- dDEL(y, mu=mu, sigma=sigma, nu=nu)
        f2 <- ifelse(f2>precision, f2, precision)
      
        S_prime_mu <- (exp(bigS(y, mu+fd.prec, sigma, nu))-exp(bigS(y, mu, sigma, nu)))/fd.prec
        S_prime_eta2 <- mu*S_prime_mu
        S <- exp(bigS(y, mu, sigma, nu))
      
        df2 <- f2*(-mu*nu+(-1/sigma)*(1+mu*sigma*(1-nu))^(-1)*mu*sigma*(1-nu)+S_prime_eta2/S)
        df2
      
      }
    
      df2.sigma.int.D <-function(y, mu, sigma, nu) {
        f2 <- dDEL(y, mu=mu, sigma=sigma, nu=nu)
        f2 <- ifelse(f2>precision, f2, precision)
      
        S_prime_sigma <- (exp(bigS(y, mu, sigma+fd.prec, nu))-exp(bigS(y, mu, sigma, nu)))/fd.prec
        d_prime <- ((1/sigma^2)*log(1+mu*sigma*(1-nu))-(1/sigma)*(1/(1+mu*sigma*(1-nu)))*mu*(1-nu))*(1+mu*sigma*(1-nu))^(-1/sigma)
        S <- exp(bigS(y, mu, sigma, nu))
      
        df2 <- f2*(d_prime*(1+mu*sigma*(1-nu))^(1/sigma)+S_prime_sigma/S)
        df2
      
      }
    
      df2.nu.int.D <- function(y, mu, sigma, nu) {
        f2 <- dDEL(y, mu=mu, sigma=sigma, nu=nu)
        f2 <- ifelse(f2>precision, f2, precision)
      
        S_prime_nu <- (exp(bigS(y, mu, sigma, nu+fd.prec))-exp(bigS(y, mu, sigma, nu)))/fd.prec
        S <- exp(bigS(y, mu, sigma, nu))
      
        df2 <- f2*(-mu+((-1/sigma)*(1+mu*sigma*(1-nu))^(-1)*(-mu*sigma)+S_prime_nu/S))
        df2
      }
    
    
      ly <- length(i2)
      cdf.f <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
      nnu <- rep(nu, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        pdfall <- df2.mu.int.D(allval, mm, ms, mn)
        cdf.f[i] <- sum(pdfall)
      }  
    
      dF2<-cdf.f
      dF22<-dF2-df2
    
    
      ly <- length(i2)
      cdf.f.sigma <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
      nnu <- rep(nu, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        pdfall <- df2.sigma.int.D(allval, mm, ms, mn)
        cdf.f.sigma[i] <- sum(pdfall)
      } 
    
      dF2.sigma <- cdf.f.sigma
      dF22.sigma <- dF2.sigma-df2.sigma
    
    
      ly <- length(i2)
      cdf.f.nu <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
      nnu <- rep(nu, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn <- nnu[i]
        allval <- seq(0, y.y)
        pdfall <- df2.nu.int.D(allval, mm, ms, mn)
        cdf.f.nu[i] <- sum(pdfall)
      } 
    
      dF2.nu <- cdf.f.nu
      dF22.nu <- dF2.nu-df2.nu
    
    }
    
    
    # Hessian derivative components
    
    #fd.prec2 <- 10^-3
    fd.prec2 <- sqrt(eps)
    
    #delta2delta2
    mu.plus <- ifelse(exp(eta2+fd.prec2)>0.0001, exp(eta2+fd.prec2), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    f2.plus <- dDEL(i2, sigma = sigma,  mu = mu.plus, nu=nu)
    f2.plus <- ifelse(f2.plus>precision, f2.plus, precision)
    S_prime_mu.plus <- (exp(bigS(i2, mu.plus+fd.prec, sigma, nu))-exp(bigS(i2, mu.plus, sigma, nu)))/fd.prec
    S_prime_eta2.plus <- mu.plus*S_prime_mu.plus
    S.plus <- exp(bigS(i2, mu.plus, sigma, nu))
    df2.plus <- as.vector(f2.plus*(-mu.plus*nu+(-1/sigma)*(1+mu.plus*sigma*(1-nu))^(-1)*mu.plus*sigma*(1-nu)+S_prime_eta2.plus/S.plus))
    df2.plus <- as.vector(df2.plus)
    d2f2delta22 <- (df2.plus - df2) / fd.prec2
    
    #sigma^2
    sigma.plus <- sigma + fd.prec2
    f2.plus.sigma <- dDEL(i2, sigma = sigma.plus,  mu = mu, nu = nu)
    f2.plus.sigma <- ifelse(f2.plus.sigma>precision, f2.plus.sigma, precision)
    S_prime_sigma <- (exp(bigS(i2, mu, sigma.plus+fd.prec, nu))-exp(bigS(i2, mu, sigma.plus, nu)))/fd.prec
    d_prime <- ((1/sigma.plus^2)*log(1+mu*sigma.plus*(1-nu))-(1/sigma.plus)*(1/(1+mu*sigma.plus*(1-nu)))*mu*(1-nu))*(1+mu*sigma.plus*(1-nu))^(-1/sigma.plus)
    S.plus.sigma <- exp(bigS(i2, mu, sigma.plus, nu))
    df2.sigma.plus <- f2.plus.sigma*(d_prime*(1+mu*sigma.plus*(1-nu))^(1/sigma.plus)+S_prime_sigma/S.plus.sigma)
    d2f2sigma2 <- (df2.sigma.plus - df2.sigma) / fd.prec2
    
    #delta2sigma
    S_prime_mu.dsigma <- (exp(bigS(i2, mu+fd.prec, sigma.plus, nu))-exp(bigS(i2, mu, sigma.plus, nu)))/fd.prec
    S_prime_eta2.dsigma <- mu*S_prime_mu.dsigma
    df2.dsigma.plus <- as.vector(f2.plus.sigma*(-mu*nu+(-1/sigma.plus)*(1+mu*sigma.plus*(1-nu))^(-1)*mu*sigma.plus*(1-nu)+S_prime_eta2.dsigma/S.plus.sigma))
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec2
    
    
    #nu^2
    nu.plus <- nu + fd.prec2
    nu.plus <- ifelse(nu.plus>=(1-precision), (1-precision), nu.plus)
    f2.plus.nu <- dDEL(i2, sigma = sigma,  mu = mu, nu = nu.plus)
    f2.plus.nu <- ifelse(f2.plus.nu>precision, f2.plus.nu, precision)
    S_prime_nu <- (exp(bigS(i2, mu, sigma, nu.plus+fd.prec))-exp(bigS(i2, mu, sigma, nu.plus)))/fd.prec
    S.plus.nu <- exp(bigS(i2, mu, sigma, nu.plus))
    df2.nu.plus <- f2.plus.nu*(-mu+((-1/sigma)*(1+mu*sigma*(1-nu.plus))^(-1)*(-mu*sigma)+S_prime_nu/S.plus.nu))
    d2f2nu2 <- (df2.nu.plus - df2.nu) / fd.prec2
    
    
    #delta2nu
    S_prime_mu.dnu <- (exp(bigS(i2, mu+fd.prec, sigma, nu.plus))-exp(bigS(i2, mu, sigma, nu.plus)))/fd.prec
    S_prime_eta2.dnu <- mu*S_prime_mu.dnu
    S.dnu <- exp(bigS(i2, mu, sigma, nu.plus))
    df2.dnu.plus <- as.vector(f2.plus.nu*(-mu*nu.plus+(-1/sigma)*(1+mu*sigma*(1-nu.plus))^(-1)*mu*sigma*(1-nu.plus)+S_prime_eta2.dnu/S.dnu))
    df2.dnu.plus <- as.vector(df2.dnu.plus)
    d2f2delta2nu <- (df2.dnu.plus - df2) / fd.prec2
    
    
    #nusigma
    f2.plus.nu.sigma <- dDEL(i2, sigma = sigma.plus,  mu = mu, nu = nu)
    f2.plus.nu.sigma <- ifelse(f2.plus.nu.sigma>precision, f2.plus.nu.sigma, precision)
    S_prime_nu.plus <- (exp(bigS(i2, mu, sigma.plus, nu+fd.prec))-exp(bigS(i2, mu, sigma.plus, nu)))/fd.prec
    df2.nu.sigma.plus <- f2.plus.nu.sigma*(-mu+((-1/sigma.plus)*(1+mu*sigma.plus*(1-nu))^(-1)*(-mu*sigma.plus)+S_prime_nu.plus/S.plus.sigma))
    d2f2nusigma <- (df2.nu.sigma.plus - df2.nu) / fd.prec2
    
    
    if(univariate==FALSE) {
    
    
      ly <- length(i2)
      cdf.f.plus <- rep(0, ly)
      nmu <- rep(mu.plus, length = ly) 
      nsigma <- rep(sigma, length = ly)
      nnu <- rep(nu, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        pdfall <- df2.mu.int.D(allval, mm, ms, mn)
        cdf.f.plus[i] <- sum(pdfall)
      }  
    
      dF2.mu.plus <- cdf.f.plus
      d2F2ddelta22 <- (dF2.mu.plus - dF2) / fd.prec2
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    
      ly <- length(i2)
      cdf.f.sigma.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma.plus, length = ly)
      nnu <- rep(nu, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        pdfall <- df2.sigma.int.D(allval, mm, ms, mn)
        cdf.f.sigma.plus[i] <- sum(pdfall)
      }  
    
      dF2.sigma.plus <- cdf.f.sigma.plus
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec2
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
    
      ly <- length(i2)
      cdf.f.dsigma.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma.plus, length = ly)
      nnu <- rep(nu, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        pdfall <- df2.mu.int.D(allval, mm, ms, mn)
        cdf.f.dsigma.plus[i] <- sum(pdfall)
      }   
    
      dF2.dsigma.plus <- cdf.f.dsigma.plus
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec2
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma   
    
      ly <- length(i2)
      cdf.f.nu.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
      nnu <- rep(nu.plus, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        pdfall <- df2.nu.int.D(allval, mm, ms, mn)
        cdf.f.nu.plus[i] <- sum(pdfall)
      } 
    
      dF2.nu.plus <- cdf.f.nu.plus
      d2F2dnu2 <- (dF2.nu.plus - dF2.nu) / fd.prec2
      d2F22dnu2 <- d2F2dnu2 - d2f2nu2
    
    
      ly <- length(i2)
      cdf.f.dnu.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
      nnu <- rep(nu.plus, length = ly)
    
     for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        pdfall <- df2.mu.int.D(allval, mm, ms, mn)
        cdf.f.dnu.plus[i] <- sum(pdfall)
      }  
    
      dF2.dnu.plus <- cdf.f.dnu.plus
      d2F2ddelta2dnu <- (dF2.dnu.plus - dF2) / fd.prec2
      d2F22ddelta2dnu <- d2F2ddelta2dnu - d2f2delta2nu
      
      ly <- length(i2)
      cdf.f.dnu.dsigma.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma.plus, length = ly)
      nnu <- rep(nu, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        pdfall <- df2.nu.int.D(allval, mm, ms, mn)
        cdf.f.dnu.dsigma.plus[i] <- sum(pdfall)
      }  
    
      dF2.dnu.dsigma.plus <- cdf.f.dnu.dsigma.plus
      d2F2dnudsigma <- (dF2.dnu.dsigma.plus - dF2.nu) / fd.prec2
      d2F22dnudsigma <- d2F2dnudsigma - d2f2nusigma
    
    }
    
    
     
  } else if (VC$margins[2]=="PIG") {
    i2 <- dat[,2]

    fd.prec <- 10^(-7)
    fd.prec2 <- 10^-3
    
    if (univariate==TRUE) {
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]    
    }
      
    sigma <- exp(sigma.star) 
    sigma <- ifelse(sigma>precision, sigma, precision)
    sigma <- ifelse(sigma<10^7, sigma, 10^7)
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^7, mu, 10^7)
    f2 <- dPIG(i2, mu=mu, sigma=sigma)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- pPIG(i2, mu=mu, sigma=sigma)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
    
    f2.fd.mu <- dPIG(i2, mu = (mu+fd.prec), sigma = sigma)
    f2.fd.mu <- ifelse(f2.fd.mu>precision, f2.fd.mu, precision)
    df2 <- (f2.fd.mu - f2) / (fd.prec)
    df2 <- as.vector(df2)*as.vector(mu)
    
    f2.fd.sigma <- dPIG(i2, mu=mu, sigma=(sigma+fd.prec))
    f2.fd.sigma <- ifelse(f2.fd.sigma>precision, f2.fd.sigma, precision)
    df2.sigma <- (f2.fd.sigma - f2) / (fd.prec)
    
    
    if(univariate==FALSE) {
      F2.fd.mu <- pPIG(i2, mu = (mu+fd.prec), sigma = sigma)
      F2.fd.mu <- ifelse(F2.fd.mu<(1-precision), F2.fd.mu, 1-precision)
      dF2 <- (F2.fd.mu-F2) / (fd.prec)
      dF2 <- as.vector(dF2)*as.vector(mu) 
      dF22 <- dF2 - df2
    
      F2.fd.sigma <- pPIG(i2, mu = mu, sigma = (sigma+fd.prec))
      F2.fd.sigma <- ifelse(F2.fd.sigma<(1-precision), F2.fd.sigma, 1-precision)
      dF2.sigma <- (F2.fd.sigma - F2) / (fd.prec)
      dF22.sigma <- dF2.sigma - df2.sigma
    }
    
    
    
    # Hessian derivative components
    
    #fd.prec2 <- 10^-3
    #fd.prec2 <- sqrt(eps)
    
    mu.plus <- ifelse(exp(eta2+fd.prec2)>0.0001, exp(eta2+fd.prec2), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    f2.plus <- dPIG(i2, mu = mu.plus, sigma = sigma)
    f2.plus <- ifelse(f2.plus>precision, f2.plus, precision)
    f2.fd.mu.plus <- dPIG(i2, mu = (mu.plus+fd.prec), sigma = sigma)
    f2.fd.mu.plus <- ifelse(f2.fd.mu.plus>precision, f2.fd.mu.plus, precision)
    df2.plus <- (f2.fd.mu.plus - f2.plus) / (fd.prec)
    df2.plus <- as.vector(df2.plus)*as.vector(mu.plus)
    df2.plus <- as.vector(df2.plus)
    d2f2delta22 <- (df2.plus - df2) / fd.prec2
    
    
    sigma.plus <- sigma + fd.prec2
    f2.plus.sigma <- dPIG(i2, mu = mu, sigma = sigma.plus)
    f2.plus.sigma <- ifelse(f2.plus.sigma>precision, f2.plus.sigma, precision)
    f2.fd.sigma.plus <- dPIG(i2, mu = mu, sigma = (sigma.plus+fd.prec))
    f2.fd.sigma.plus <- ifelse(f2.fd.sigma.plus>precision, f2.fd.sigma.plus, precision)
    df2.sigma.plus <- (f2.fd.sigma.plus - f2.plus.sigma) / (fd.prec)
    d2f2sigma2 <- (df2.sigma.plus - df2.sigma) / fd.prec2
    
    
    f2.dsigma.plus <- dPIG(i2, mu = mu, sigma = sigma.plus)
    f2.dsigma.plus <- ifelse(f2.dsigma.plus>precision, f2.dsigma.plus, precision)
    f2.fd.mu.sigma.plus <- dPIG(i2, mu = (mu+fd.prec), sigma = sigma.plus)
    f2.fd.mu.sigma.plus <- ifelse(f2.fd.mu.sigma.plus>precision, f2.fd.mu.sigma.plus, precision)
    df2.dsigma.plus <- (f2.fd.mu.sigma.plus - f2.dsigma.plus) / (fd.prec)
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)*as.vector(mu)
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec2
    
    if(univariate==FALSE) {
    
      F2.plus <- pPIG(i2, mu = mu.plus, sigma = sigma)
      F2.fd.mu.plus <- pPIG(i2, mu = (mu.plus+fd.prec), sigma = sigma)
      F2.fd.mu.plus <- ifelse(F2.fd.mu.plus<(1-precision), F2.fd.mu.plus, 1-precision)
      dF2.plus <- (F2.fd.mu.plus - F2.plus) / (fd.prec)
      dF2.plus <- as.vector(dF2.plus)*as.vector(mu.plus) 
      d2F2ddelta22 <- (dF2.plus - dF2) / fd.prec2
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
      F2.sigma.plus <- pPIG(i2, mu = mu, sigma = sigma.plus)
      F2.fd.sigma.plus <- pPIG(i2, mu = mu, sigma = (sigma.plus+fd.prec))
      F2.fd.sigma.plus <- ifelse(F2.fd.sigma.plus<(1-precision), F2.fd.sigma.plus, 1-precision)
      dF2.sigma.plus <- (F2.fd.sigma.plus - F2.sigma.plus) / (fd.prec)
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) /fd.prec2
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
      F2.fd.dsigma.plus <- pPIG(i2, mu = (mu+fd.prec), sigma = sigma.plus)
      F2.fd.dsigma.plus <- ifelse(F2.fd.dsigma.plus<(1-precision), F2.fd.dsigma.plus, 1-precision)
      dF2.dsigma.plus <- (F2.fd.dsigma.plus - F2.sigma.plus) / (fd.prec)
      dF2.dsigma.plus <- as.vector(dF2.dsigma.plus)*as.vector(mu) 
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec2
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma   
    
    }
    
    
  } else if (VC$margins[2]=="S") {
    i2 <- dat[,2]
    
    fd.prec <- 10^(-7)
    fd.prec2 <- 10^-3
    
    if(univariate==TRUE){
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
      if (is.null(VC$X3)) nu <- nu.star <- params[(VC$X2.d2 + 2)]
      if (!is.null(VC$X3)) nu <- nu.star <- etanu <- X4 %*% params[(VC$X2.d2 + VC$X3.d2 + 1):(VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]        
      if(is.null(VC$X3))  nu <- nu.star <- params[(VC$X1.d2+VC$X2.d2+2)]
      if(!is.null(VC$X3)) nu <- nu.star <- etanu <- X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)]     
    }
      
    sigma <- exp(sigma.star)
    sigma <- ifelse(sigma>precision, sigma, precision)
    sigma <- ifelse(sigma<10^7, sigma, 10^7)
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^7, mu, 10^7) 
    nu <- ifelse(nu<150, nu, 150)
    nu <- ifelse(nu>-150, nu, -150)    
    f2 <- dSICHEL(i2, mu=mu, sigma=sigma, nu=nu)
    f2 <- ifelse(f2>precision, f2, precision)
  
    if(univariate==FALSE) {
      F2 <- pSICHEL(i2, mu=mu, sigma=sigma, nu=nu)
      F2 <- ifelse(F2>precision, F2, precision)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
    
    
    f2.fd.mu <- dSICHEL(i2, mu = (mu +fd.prec), sigma = sigma, nu = nu)
    f2.fd.mu <- ifelse(f2.fd.mu>precision, f2.fd.mu, precision)
    df2 <- (f2.fd.mu - f2) / (fd.prec)
    df2 <- as.vector(df2)*as.vector(mu)
    
    f2.fd.sigma <- dSICHEL(i2, mu = mu, sigma = (sigma+fd.prec), nu = nu)
    f2.fd.sigma <- ifelse(f2.fd.sigma>precision, f2.fd.sigma, precision)
    df2.sigma <- (f2.fd.sigma - f2) / (fd.prec)
    
    f2.fd.nu <- dSICHEL(i2, mu = mu, sigma = sigma, nu =(nu+fd.prec))
    f2.fd.nu <- ifelse(f2.fd.nu>precision, f2.fd.nu, precision)
    df2.nu <- (f2.fd.nu - f2) / (fd.prec)
    
    if(univariate==FALSE) {
      F2.fd.mu <- pSICHEL(i2, mu = (mu+fd.prec), sigma = sigma, nu = nu)
      F2.fd.mu <- ifelse(F2.fd.mu<(1-precision), F2.fd.mu, 1-precision)
      dF2 <- (F2.fd.mu - F2) / (fd.prec)
      dF2 <- as.vector(dF2)*as.vector(mu) 
      dF22 <- dF2 - df2
    
      F2.fd.sigma <- pSICHEL(i2, mu = mu, sigma = (sigma+fd.prec), nu = nu)
      F2.fd.sigma <- ifelse(F2.fd.sigma<(1-precision), F2.fd.sigma, 1-precision)
      dF2.sigma <- (F2.fd.sigma - F2) / (fd.prec)
      dF22.sigma <- dF2.sigma - df2.sigma
    
      F2.fd.nu <- pSICHEL(i2, mu = mu, sigma = sigma, nu = (nu+fd.prec))
      F2.fd.nu <- ifelse(F2.fd.nu<(1-precision), F2.fd.nu, 1-precision)
      dF2.nu <- (F2.fd.nu - F2) / (fd.prec)
      dF22.nu <- dF2.nu - df2.nu
    }
    
    
    # Hessian derivative components
    
    #second derivative fd precision
    #fd.prec2 <- 10^-3
    #fd.prec2 <- sqrt(eps)
    
    #pmfs
    
    #delta2^2
    mu.plus <- ifelse(exp(eta2+fd.prec2)>0.0001, exp(eta2+fd.prec2), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    f2.plus <- dSICHEL(i2, mu = mu.plus, sigma = sigma, nu = nu)
    f2.plus <- ifelse(f2.plus>precision, f2.plus, precision)
    f2.fd.mu.plus <- dSICHEL(i2, mu = (mu.plus+fd.prec), sigma = sigma, nu = nu)
    f2.fd.mu.plus <- ifelse(f2.fd.mu.plus>precision, f2.fd.mu.plus, precision)
    df2.plus <- (f2.fd.mu.plus - f2.plus) / (fd.prec)
    df2.plus <- as.vector(df2.plus)*as.vector(mu.plus)
    df2.plus <- as.vector(df2.plus)
    d2f2delta22 <- (df2.plus - df2) / fd.prec2
    
    #sigma^2
    sigma.plus <- sigma+fd.prec2
    f2.plus.sigma <- dSICHEL(i2, mu = mu, sigma = sigma.plus, nu = nu)
    f2.plus.sigma <- ifelse(f2.plus.sigma>precision, f2.plus.sigma, precision)
    f2.fd.sigma.plus <- dSICHEL(i2, mu = mu, sigma = (sigma.plus+fd.prec), nu = nu)
    f2.fd.sigma.plus <- ifelse(f2.fd.sigma.plus>precision, f2.fd.sigma.plus, precision)
    df2.sigma.plus <- (f2.fd.sigma.plus - f2.plus.sigma) / (fd.prec)
    d2f2sigma2 <- (df2.sigma.plus - df2.sigma) / fd.prec2
    
    #delta2sigma
    f2.dsigma.plus <- dSICHEL(i2, mu = mu, sigma = sigma.plus, nu = nu)
    f2.dsigma.plus <- ifelse(f2.dsigma.plus>precision, f2.dsigma.plus, precision)
    f2.fd.mu.sigma.plus <- dSICHEL(i2, mu = (mu+fd.prec), sigma = sigma.plus, nu = nu)
    f2.fd.mu.sigma.plus <- ifelse(f2.fd.mu.sigma.plus>precision, f2.fd.mu.sigma.plus, precision)
    df2.dsigma.plus <- (f2.fd.mu.sigma.plus - f2.dsigma.plus) / (fd.prec)
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)*as.vector(mu)
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec2
    
    #nu^2
    nu.plus <- nu+fd.prec2
    f2.plus.nu <- dSICHEL(i2, mu = mu, sigma = sigma, nu = nu.plus)
    f2.plus.nu <- ifelse(f2.plus.nu>precision, f2.plus.nu, precision)
    f2.fd.nu.plus <- dSICHEL(i2, mu = mu, sigma = sigma, nu = nu.plus+fd.prec)
    f2.fd.nu.plus <- ifelse(f2.fd.nu.plus>precision, f2.fd.nu.plus, precision)
    df2.nu.plus <- (f2.fd.nu.plus - f2.plus.nu) / (fd.prec)
    d2f2nu2 <- (df2.nu.plus - df2.nu) / fd.prec2
    
    #delta2nu
    f2.dnu.plus <- dSICHEL(i2, mu = mu, sigma = sigma, nu = nu.plus)
    f2.dnu.plus <- ifelse(f2.dnu.plus>precision, f2.dnu.plus, precision)
    f2.fd.mu.nu.plus <- dSICHEL(i2, mu = (mu+fd.prec), sigma = sigma, nu = nu.plus)
    f2.fd.mu.nu.plus <- ifelse(f2.fd.mu.nu.plus>precision, f2.fd.mu.nu.plus, precision)
    df2.dnu.plus <- (f2.fd.mu.nu.plus - f2.dnu.plus) / (fd.prec)
    df2.dnu.plus <- as.vector(df2.dnu.plus)*as.vector(mu)
    d2f2delta2nu <-(df2.dnu.plus - df2) / fd.prec2
    
    #nusigma
    f2.fd.sigma.nu.plus <- dSICHEL(i2, mu = mu, sigma = (sigma+fd.prec), nu = nu.plus)
    f2.fd.sigma.nu.plus <- ifelse(f2.fd.sigma.nu.plus>precision, f2.fd.sigma.nu.plus, precision)
    df2.dnu.sigma.plus <- (f2.fd.sigma.nu.plus - f2.plus.nu) / (fd.prec)
    df2.dnu.sigma.plus <- as.vector(df2.dnu.sigma.plus)
    d2f2sigmanu <- (df2.dnu.sigma.plus - df2.sigma) / fd.prec2
    
    
    if(univariate==FALSE) {
    
      #cdfs
    
      #delta2^2
      F2.plus <- pSICHEL(i2, mu = mu.plus, sigma = sigma, nu = nu)
      F2.fd.mu.plus <- pSICHEL(i2, mu = (mu.plus+fd.prec), sigma = sigma, nu = nu)
      F2.fd.mu.plus <- ifelse(F2.fd.mu.plus<(1-precision), F2.fd.mu.plus, 1-precision)
      dF2.plus <- (F2.fd.mu.plus - F2.plus) / (fd.prec)
      dF2.plus <- as.vector(dF2.plus)*as.vector(mu.plus) 
      d2F2ddelta22 <-(dF2.plus - dF2) / fd.prec2
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
      #sigma^2
      F2.sigma.plus <- pSICHEL(i2, mu = mu, sigma = sigma.plus, nu = nu)
      F2.fd.sigma.plus <- pSICHEL(i2, mu = mu, sigma = (sigma.plus+fd.prec), nu=nu)
      F2.fd.sigma.plus <- ifelse(F2.fd.sigma.plus<(1-precision), F2.fd.sigma.plus, 1-precision)
      dF2.sigma.plus <- (F2.fd.sigma.plus - F2.sigma.plus) / (fd.prec)
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec2
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
      #delta2sigma
      F2.fd.dsigma.plus <- pSICHEL(i2, mu = (mu+fd.prec), sigma = sigma.plus, nu)
      F2.fd.dsigma.plus <- ifelse(F2.fd.dsigma.plus<(1-precision), F2.fd.dsigma.plus, 1-precision)
      dF2.dsigma.plus <- (F2.fd.dsigma.plus - F2.sigma.plus) / (fd.prec)
      dF2.dsigma.plus <- as.vector(dF2.dsigma.plus)*as.vector(mu) 
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec2
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma
    
      #nu^2
      F2.nu.plus<- pSICHEL(i2, mu = mu, sigma = sigma, nu = nu.plus)
      F2.fd.nu.plus <- pSICHEL(i2, mu = mu, sigma = sigma, nu = nu.plus+fd.prec)
      F2.fd.nu.plus <- ifelse(F2.fd.nu.plus<(1-precision), F2.fd.nu.plus, 1-precision)
      dF2.nu.plus <- (F2.fd.nu.plus - F2.nu.plus) / (fd.prec)
      d2F2dnu2 <- (dF2.nu.plus - dF2.nu) / fd.prec2
      d2F22dnu2 <- d2F2dnu2 - d2f2nu2
    
      #delta2nu
      F2.fd.dnu.plus <- pSICHEL(i2, mu = (mu+fd.prec), sigma = sigma, nu.plus)
      F2.fd.dnu.plus <- ifelse(F2.fd.dnu.plus<(1-precision), F2.fd.dnu.plus, 1-precision)
      dF2.dnu.plus <- (F2.fd.dnu.plus - F2.nu.plus) / (fd.prec)
      dF2.dnu.plus <- as.vector(dF2.dnu.plus)*as.vector(mu) 
      d2F2ddelta2dnu <- (dF2.dnu.plus - dF2) / fd.prec2
      d2F22ddelta2dnu <- d2F2ddelta2dnu - d2f2delta2nu
    
      #sigmanu
      F2.fd.dnu.sigma.plus <- pSICHEL(i2, mu = mu, sigma = (sigma+fd.prec), nu.plus)
      F2.fd.dnu.sigma.plus <- ifelse(F2.fd.dnu.sigma.plus<(1-precision), F2.fd.dnu.sigma.plus, 1-precision)
      dF2.dnu.sigma.plus <- (F2.fd.dnu.sigma.plus - F2.nu.plus) / (fd.prec)
      dF2.dnu.sigma.plus <- as.vector(dF2.dnu.sigma.plus)
      d2F2dsigmadnu <- (dF2.dnu.sigma.plus - dF2.sigma) / fd.prec2
      d2F22dsigmadnu <- d2F2dsigmadnu - d2f2sigmanu
    
    }
    
    
  } else if(VC$margins[2]=="BI") {
  
    i2 <- dat[,2]
    mu <- plogis(eta2)
    mu <- ifelse(mu<(1-precision), mu, 1-mu)
    mu <- ifelse(mu>precision, mu, precision) 
    f2 <- dBI(i2, mu, bd = bd)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- pBI(i2, mu, bd = bd)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
    
    
    # df2, dF2, dF22  
    df2 <- -(((1 - mu)^(-1 + bd - i2)*mu^(-1 + i2)*(-i2 + bd*mu)*gamma(1 + bd))/(gamma(1 + bd - i2)*gamma(1 + i2)))*(1-mu)*mu
    
    
   if(univariate==FALSE) {
  
      df2.int.BI <- function(y, mu, n) { 
        -(((1 - mu)^(-1 + n - y)*mu^(-1 + y)*(-y + n*mu)*gamma(1 + n))/(gamma(1 + n - y)*gamma(1 + y)))*(1-mu)*mu
      }
    
      ly <- length(i2)
      cdf.f <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        allval <- seq(0, y.y)
        pdfall <- df2.int.BI(allval, mu = mm, n=bd)
        cdf.f[i] <- sum(pdfall)
      }  
    
      dF2 <- cdf.f
      dF22 <- dF2-df2
    
   }
    
    # Hessian derivative components
    
    mu.plus <- plogis(eta2+fd.prec)
    df2.plus <- -(((1 - mu.plus)^(-1 + bd - i2)*mu.plus^(-1 + i2)*(-i2 + bd*mu.plus)*gamma(1 + bd))/(gamma(1 + bd - i2)*gamma(1 + i2)))*(1-mu.plus)*mu.plus
    d2f2delta22 <- (df2.plus - df2)/fd.prec
    
    if(univariate==FALSE) {
   
      ly <- length(i2)
      cdf.f.plus <- rep(0, ly)
      nmu.plus <- rep(mu.plus, length = ly) 
    
      for (i in 1:ly) {
      
        y.y <- i2[i]
        mm <- nmu.plus[i]
        allval <- seq(0, y.y)
        pdfall <- df2.int.BI(allval, mu = mm, n=bd)
        cdf.f.plus[i] <- sum(pdfall)
      
      }  
    
      dF2.plus <- cdf.f.plus
      d2F2ddelta22 <- (dF2.plus - dF2) / fd.prec
    
    
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    }
    
    
    
    
  } else if (VC$margins[2]=="BB") {
    
    i2 <- dat[,2]
    
    if(univariate==TRUE) {
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
    }
      
    sigma <- exp(sigma.star)
    sigma <- ifelse(sigma>10^-8, sigma, 10^-8)
    sigma <- ifelse(sigma<10^7, sigma, 10^7)
    mu <- plogis(eta2)
    mu <- ifelse(mu<(1-precision), mu, 1-mu)
    mu <- ifelse(mu>precision, mu, precision)
    f2 <- dBB(i2, sigma = sigma,  mu = mu, bd=bd)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- pBB(i2, sigma = sigma, mu = mu, bd=bd)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
    
    
    df2 <- ( (gamma(1 + bd) * gamma(i2 + mu/sigma) * gamma(1/sigma) * gamma(
            bd - (-1 + mu + i2 *  sigma)/sigma) * (digamma(i2 + mu/sigma) + 
            digamma((1 - mu)/sigma) - digamma(mu/sigma) - 
            digamma(bd - (-1 + mu + i2 * sigma)/sigma)))/(sigma * gamma(
            1 + bd - i2) * gamma(1 + i2) * gamma(bd + 1/sigma) * gamma((1 - mu)/
            sigma) * gamma(mu/sigma)) )*(1-mu)*mu
    df2 <- as.vector(df2)
    
    df2.sigma <- (gamma(1 + bd) * gamma(i2 + mu/sigma) * gamma(1/sigma) * gamma(
                  bd - (-1 + mu + i2 * sigma)/sigma) * (digamma(bd + 1/sigma) - 
                  mu * digamma(i2 + mu/sigma) - digamma(1/sigma) + 
                  digamma((1 - mu)/sigma) - mu * digamma((1 - mu)/sigma) + 
                  mu * digamma(mu/sigma) - 
                  digamma(bd - (-1 + mu + i2 * sigma)/sigma) + 
                  mu * digamma(bd - (-1 + mu + i2 * sigma)/sigma)))/(sigma^2 * gamma(
                  1 + bd - i2) * gamma(1 + i2) * gamma(bd + 1/sigma) * gamma((1 - mu)/
                  sigma) * gamma(mu/sigma))
    
    
    if(univariate==FALSE) {
    
      df2.mu.int.BB<-function(y, mu, sigma, nu, n=bd) {
        
        df2<- ( (gamma(1 + n) * gamma(y + mu/sigma) * gamma(1/sigma) * gamma(
                n - (-1 + mu + y *  sigma)/sigma) * (digamma(y + mu/sigma) + 
                digamma((1 - mu)/sigma) - digamma(mu/sigma) - 
                digamma(n - (-1 + mu + y * sigma)/sigma)))/(sigma * gamma(
                1 + n - y) * gamma(1 + y) * gamma(n + 1/sigma) * gamma((1 - mu)/
                sigma) * gamma(mu/sigma)) )*(1-mu)*mu
       df2
      }
    
      df2.sigma.int.BB<-function(y, mu, sigma, nu, n=bd) {
        
        df2 <- (gamma(1 + n) * gamma(y + mu/sigma) * gamma(1/sigma) * gamma(
                  n - (-1 + mu + y * sigma)/sigma) * (digamma(n + 1/sigma) - 
                  mu * digamma(y + mu/sigma) - digamma(1/sigma) + 
                  digamma((1 - mu)/sigma) - mu * digamma((1 - mu)/sigma) + 
                  mu * digamma(mu/sigma) - 
                  digamma(n - (-1 + mu + y * sigma)/sigma) + 
                  mu * digamma(n - (-1 + mu + y * sigma)/sigma)))/(sigma^2 * gamma(
                  1 + n - y) * gamma(1 + y) * gamma(n + 1/sigma) * gamma((1 - mu)/
                  sigma) * gamma(mu/sigma))
        df2
      }
    
      ly <- length(i2)
      cdf.f <- rep(0, ly)
      cdf.f.sigma <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f[i] <- sum(df2.mu.int.BB(allval, mm, ms, bd))
        cdf.f.sigma[i] <- sum(df2.sigma.int.BB(allval, mm, ms, bd))
      }  
    
      dF2<-cdf.f
      dF22<-dF2-df2
    
      dF2.sigma <- cdf.f.sigma
      dF22.sigma <- dF2.sigma-df2.sigma
    
    }
    
    # Hessian derivative components
    
    mu.plus <- ifelse(plogis(eta2+fd.prec)>0.0001, plogis(eta2+fd.prec), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    df2.plus <-( (gamma(1 + bd) * gamma(i2 + mu.plus/sigma) * gamma(1/sigma) * gamma(
                  bd - (-1 + mu.plus + i2 *  sigma)/sigma) * (digamma(i2 + mu.plus/sigma) + 
                  digamma((1 - mu.plus)/sigma) - digamma(mu.plus/sigma) - 
                  digamma(bd - (-1 + mu.plus + i2 * sigma)/sigma)))/(sigma * gamma(
                  1 + bd - i2) * gamma(1 + i2) * gamma(bd + 1/sigma) * gamma((1 - mu.plus)/
                  sigma) * gamma(mu.plus/sigma)) )*(1-mu.plus)*mu.plus
    df2.plus <- as.vector(df2.plus)
    d2f2delta22 <- (df2.plus - df2) / fd.prec
    
    
    sigma.plus <- sigma+fd.prec
    df2.sigma.plus <- (gamma(1 + bd) * gamma(i2 + mu/sigma.plus) * gamma(1/sigma.plus) * gamma(
                        bd - (-1 + mu + i2 * sigma.plus)/sigma.plus) * (digamma(bd + 1/sigma.plus) - 
                        mu * digamma(i2 + mu/sigma.plus) - digamma(1/sigma.plus) + 
                        digamma((1 - mu)/sigma.plus) - mu * digamma((1 - mu)/sigma.plus) + 
                        mu * digamma(mu/sigma.plus) - 
                        digamma(bd - (-1 + mu + i2 * sigma.plus)/sigma.plus) + 
                        mu * digamma(bd - (-1 + mu + i2 * sigma.plus)/sigma.plus)))/(sigma.plus^2 * gamma(
                        1 + bd - i2) * gamma(1 + i2) * gamma(bd + 1/sigma.plus) * gamma((1 - mu)/
                        sigma.plus) * gamma(mu/sigma.plus))
    d2f2sigma2 <- (df2.sigma.plus - df2.sigma) / fd.prec
    
    
    df2.dsigma.plus <-(  (gamma(1 + bd) * gamma(i2 + mu/sigma.plus) * gamma(1/sigma.plus) * gamma(
                          bd - (-1 + mu + i2 *  sigma.plus)/sigma.plus) * (digamma(i2 + mu/sigma.plus) + 
                          digamma((1 - mu)/sigma.plus) - digamma(mu/sigma.plus) - 
                          digamma(bd - (-1 + mu + i2 * sigma.plus)/sigma.plus)))/(sigma.plus * gamma(
                          1 + bd - i2) * gamma(1 + i2) * gamma(bd + 1/sigma.plus) * gamma((1 - mu)/
                          sigma.plus) * gamma(mu/sigma.plus)) )*(1-mu)*mu
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec
    
    if(univariate==FALSE) {
    
      cdf.f.plus <- rep(0, ly)
      nmu <- rep(mu.plus, length = ly) 
      nsigma <- rep(sigma, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        pdfall <- df2.mu.int.BB(allval, mm, ms, bd)
        cdf.f.plus[i] <- sum(pdfall)
      }  
    
      dF2.mu.plus <- cdf.f.plus
      d2F2ddelta22 <- (dF2.mu.plus - dF2) / fd.prec
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    
      cdf.f.sigma.plus <- rep(0, ly)
      cdf.f.dsigma.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma.plus, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f.sigma.plus[i] <- sum(df2.sigma.int.BB(allval, mm, ms, mn, bd))
        cdf.f.dsigma.plus[i] <- sum(df2.mu.int.BB(allval, mm, ms, mn, bd))
      }  
    
      dF2.sigma.plus <- cdf.f.sigma.plus
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
      dF2.dsigma.plus <- cdf.f.dsigma.plus
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma  
    
    }
    
    
  } else if(VC$margins[2]=="GEOM") {
    i2 <- dat[,2]
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^7, mu, 10^7)    
    f2 <- dGEOM(i2, mu)
    f2 <- ifelse(f2>precision, f2, precision)
    
    
    if(univariate==FALSE) {
      F2 <- pGEOM(i2, mu)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
    
    # df2, dF2, dF22
    p <- 1/(mu+1)
    df2 <- ((1 - p)^i2 - (1 - p)^(i2 - 1) * i2 * p)*(-(1/(mu + 1)^2))*mu
    
    
   if(univariate==FALSE) {
  
      df2.int.GEOM <- function(y, mu) {
        p <- 1/(mu+1)
        ((1 - p)^y - (1 - p)^(y - 1) * y * p)*(-(1/(mu + 1)^2))*mu
      }
    
      ly <- length(i2)
      cdf.f <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        allval <- seq(0, y.y)
        pdfall <- df2.int.GEOM(allval, mu = mm)
        cdf.f[i] <- sum(pdfall)
      }  
    
      dF2 <- cdf.f
      dF22 <- dF2-df2
    
   }
    
    # Hessian derivative components
    
    mu.plus <- exp(eta2+fd.prec)
    p.plus <- 1/(mu.plus+1)
    df2.plus <- ((1 - p.plus)^i2 - (1 - p.plus)^(i2 - 1) * i2 * p.plus)*(-(1/(mu.plus + 1)^2))*mu.plus
    d2f2delta22 <- (df2.plus - df2)/fd.prec
    
    if(univariate==FALSE) {
   
      ly <- length(i2)
      cdf.f.plus <- rep(0, ly)
      nmu.plus <- rep(mu.plus, length = ly) 
    
      for (i in 1:ly) {
      
        y.y <- i2[i]
        mm <- nmu.plus[i]
        allval <- seq(0, y.y)
        pdfall <- df2.int.GEOM(allval, mu = mm)
        cdf.f.plus[i] <- sum(pdfall)
      
      }  
    
      dF2.plus <- cdf.f.plus
      d2F2ddelta22 <- (dF2.plus - dF2) / fd.prec
    
    
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    }
    
    
  } else if(VC$margins[2]=="LG") {
  
    i2 <- dat[,2]
    mu <- plogis(eta2)
    mu <- ifelse(mu<(1-precision), mu, 1-mu)
    mu <- ifelse(mu>precision, mu, precision)   
    f2 <- dLG(i2, mu)
    f2 <- ifelse(f2>precision, f2, precision)
    
    
    if(univariate==FALSE) {
      F2 <- pLG(i2, mu)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision) 
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
    
    
    # df2, dF2, dF22
    df2 <- ( ((log(1 - mu))^((-1) - 1) * ((-1) * (1/(1 - mu))) * mu^i2 + (-(log(1 - mu))^(-1)) * (mu^(i2 - 1) * i2))/i2 )*(1-mu)*mu
    
    
   if(univariate==FALSE) {
  
      df2.int.LG <- function(y, mu) { 
        ( ((log(1 - mu))^((-1) - 1) * ((-1) * (1/(1 - mu))) * mu^y + (-(log(1 - mu))^(-1)) * (mu^(y - 1) * y))/y )*(1-mu)*mu
      }
    
      ly <- length(i2)
      cdf.f <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        allval <- seq(1, y.y)
        pdfall <- df2.int.LG(allval, mu = mm)
        cdf.f[i] <- sum(pdfall)
      }  
    
      dF2 <- cdf.f
      dF22 <- dF2-df2
    
   }
    
    # Hessian derivative components
    
    mu.plus <- plogis(eta2+fd.prec)
    df2.plus <- ( ((log(1 - mu.plus))^((-1) - 1) * ((-1) * (1/(1 - mu.plus))) * mu.plus^i2 + (-(log(1 - mu.plus))^(-1)) * (mu.plus^(i2 - 1) * i2))/i2 )*(1-mu.plus)*mu.plus
    d2f2delta22 <- (df2.plus - df2)/fd.prec
    
    if(univariate==FALSE) {
   
      ly <- length(i2)
      cdf.f.plus <- rep(0, ly)
      nmu.plus <- rep(mu.plus, length = ly) 
    
      for (i in 1:ly) {
      
        y.y <- i2[i]
        mm <- nmu.plus[i]
        allval <- seq(1, y.y)
        pdfall <- df2.int.LG(allval, mu = mm)
        cdf.f.plus[i] <- sum(pdfall)
      
      }  
    
      dF2.plus <- cdf.f.plus
      d2F2ddelta22 <- (dF2.plus - dF2) / fd.prec
    
    
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    }
    
  
    
  } else if (VC$margins[2]=="NBII") {
    i2 <- dat[,2]

    if (univariate==TRUE) {
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]    
    }
      
    sigma <- exp(sigma.star) 
    sigma <- ifelse(sigma>precision, sigma, precision)
    sigma <- ifelse(sigma<10^7, sigma, 10^7)
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^7, mu, 10^7)
    f2 <- dNBII(i2, mu=mu, sigma=sigma)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- pNBII(i2, mu=mu, sigma=sigma)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
    
    f2.fd.mu <- dNBII(i2, mu = (mu+fd.prec), sigma = sigma)
    f2.fd.mu <- ifelse(f2.fd.mu>precision, f2.fd.mu, precision)
    df2 <- (f2.fd.mu - f2) / (fd.prec)
    df2 <- as.vector(df2)*as.vector(mu)
    
    f2.fd.sigma <- dNBII(i2, mu=mu, sigma=(sigma+fd.prec))
    f2.fd.sigma <- ifelse(f2.fd.sigma>precision, f2.fd.sigma, precision)
    df2.sigma <- (f2.fd.sigma - f2) / (fd.prec)
    
    
    if(univariate==FALSE) {
      F2.fd.mu <- pNBII(i2, mu = (mu+fd.prec), sigma = sigma)
      F2.fd.mu <- ifelse(F2.fd.mu<(1-precision), F2.fd.mu, 1-precision)
      dF2 <- (F2.fd.mu-F2) / (fd.prec)
      dF2 <- as.vector(dF2)*as.vector(mu) 
      dF22 <- dF2 - df2
    
      F2.fd.sigma <- pNBII(i2, mu = mu, sigma = (sigma+fd.prec))
      F2.fd.sigma <- ifelse(F2.fd.sigma<(1-precision), F2.fd.sigma, 1-precision)
      dF2.sigma <- (F2.fd.sigma - F2) / (fd.prec)
      dF22.sigma <- dF2.sigma - df2.sigma
    }
    
    
    
    # Hessian derivative components
    
    #fd.prec2 <- 10^-3
    fd.prec2 <- sqrt(eps)
    
    mu.plus <- ifelse(exp(eta2+fd.prec2)>0.0001, exp(eta2+fd.prec2), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    f2.plus <- dNBII(i2, mu = mu.plus, sigma = sigma)
    f2.plus <- ifelse(f2.plus>precision, f2.plus, precision)
    f2.fd.mu.plus <- dNBII(i2, mu = (mu.plus+fd.prec), sigma = sigma)
    f2.fd.mu.plus <- ifelse(f2.fd.mu.plus>precision, f2.fd.mu.plus, precision)
    df2.plus <- (f2.fd.mu.plus - f2.plus) / (fd.prec)
    df2.plus <- as.vector(df2.plus)*as.vector(mu.plus)
    df2.plus <- as.vector(df2.plus)
    d2f2delta22 <- (df2.plus - df2) / fd.prec2
    
    
    sigma.plus <- sigma + fd.prec2
    f2.plus.sigma <- dNBII(i2, mu = mu, sigma = sigma.plus)
    f2.plus.sigma <- ifelse(f2.plus.sigma>precision, f2.plus.sigma, precision)
    f2.fd.sigma.plus <- dNBII(i2, mu = mu, sigma = (sigma.plus+fd.prec))
    f2.fd.sigma.plus <- ifelse(f2.fd.sigma.plus>precision, f2.fd.sigma.plus, precision)
    df2.sigma.plus <- (f2.fd.sigma.plus - f2.plus.sigma) / (fd.prec)
    d2f2sigma2 <- (df2.sigma.plus - df2.sigma) / fd.prec2
    
    
    f2.dsigma.plus <- dNBII(i2, mu = mu, sigma = sigma.plus)
    f2.dsigma.plus <- ifelse(f2.dsigma.plus>precision, f2.dsigma.plus, precision)
    f2.fd.mu.sigma.plus <- dNBII(i2, mu = (mu+fd.prec), sigma = sigma.plus)
    f2.fd.mu.sigma.plus <- ifelse(f2.fd.mu.sigma.plus>precision, f2.fd.mu.sigma.plus, precision)
    df2.dsigma.plus <- (f2.fd.mu.sigma.plus - f2.dsigma.plus) / (fd.prec)
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)*as.vector(mu)
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec2
    
    if(univariate==FALSE) {
    
      F2.plus <- pNBII(i2, mu = mu.plus, sigma = sigma)
      F2.fd.mu.plus <- pNBII(i2, mu = (mu.plus+fd.prec), sigma = sigma)
      F2.fd.mu.plus <- ifelse(F2.fd.mu.plus<(1-precision), F2.fd.mu.plus, 1-precision)
      dF2.plus <- (F2.fd.mu.plus - F2.plus) / (fd.prec)
      dF2.plus <- as.vector(dF2.plus)*as.vector(mu.plus) 
      d2F2ddelta22 <- (dF2.plus - dF2) / fd.prec2
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
      F2.sigma.plus <- pNBII(i2, mu = mu, sigma = sigma.plus)
      F2.fd.sigma.plus <- pNBII(i2, mu = mu, sigma = (sigma.plus+fd.prec))
      F2.fd.sigma.plus <- ifelse(F2.fd.sigma.plus<(1-precision), F2.fd.sigma.plus, 1-precision)
      dF2.sigma.plus <- (F2.fd.sigma.plus - F2.sigma.plus) / (fd.prec)
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) /fd.prec2
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
      F2.fd.dsigma.plus <- pNBII(i2, mu = (mu+fd.prec), sigma = sigma.plus)
      F2.fd.dsigma.plus <- ifelse(F2.fd.dsigma.plus<(1-precision), F2.fd.dsigma.plus, 1-precision)
      dF2.dsigma.plus <- (F2.fd.dsigma.plus - F2.sigma.plus) / (fd.prec)
      dF2.dsigma.plus <- as.vector(dF2.dsigma.plus)*as.vector(mu) 
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec2
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma   
    
    }
    
    
    
  } else if (VC$margins[2]=="WARING") {
    i2 <- dat[,2]

    if (univariate==TRUE) {
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]    
    }
      
    sigma <- exp(sigma.star) 
    sigma <- ifelse(sigma>precision, sigma, precision)
    sigma <- ifelse(sigma<10^7, sigma, 10^7)
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^7, mu, 10^7)
    f2 <- dWARING(i2, mu=mu, sigma=sigma)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- pWARING(i2, mu=mu, sigma=sigma)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
    
    f2.fd.mu <- dWARING(i2, mu = (mu+fd.prec), sigma = sigma)
    f2.fd.mu <- ifelse(f2.fd.mu>precision, f2.fd.mu, precision)
    df2 <- (f2.fd.mu - f2) / (fd.prec)
    df2 <- as.vector(df2)*as.vector(mu)
    
    f2.fd.sigma <- dWARING(i2, mu=mu, sigma=(sigma+fd.prec))
    f2.fd.sigma <- ifelse(f2.fd.sigma>precision, f2.fd.sigma, precision)
    df2.sigma <- (f2.fd.sigma - f2) / (fd.prec)
    
    
    if(univariate==FALSE) {
      F2.fd.mu <- pWARING(i2, mu = (mu+fd.prec), sigma = sigma)
      F2.fd.mu <- ifelse(F2.fd.mu<(1-precision), F2.fd.mu, 1-precision)
      dF2 <- (F2.fd.mu-F2) / (fd.prec)
      dF2 <- as.vector(dF2)*as.vector(mu) 
      dF22 <- dF2 - df2
    
      F2.fd.sigma <- pWARING(i2, mu = mu, sigma = (sigma+fd.prec))
      F2.fd.sigma <- ifelse(F2.fd.sigma<(1-precision), F2.fd.sigma, 1-precision)
      dF2.sigma <- (F2.fd.sigma - F2) / (fd.prec)
      dF22.sigma <- dF2.sigma - df2.sigma
    }
    
    
    
    # Hessian derivative components
    
    #fd.prec2 <- 10^-3
    fd.prec2 <- sqrt(eps)
    
    mu.plus <- ifelse(exp(eta2+fd.prec2)>0.0001, exp(eta2+fd.prec2), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    f2.plus <- dWARING(i2, mu = mu.plus, sigma = sigma)
    f2.plus <- ifelse(f2.plus>precision, f2.plus, precision)
    f2.fd.mu.plus <- dWARING(i2, mu = (mu.plus+fd.prec), sigma = sigma)
    f2.fd.mu.plus <- ifelse(f2.fd.mu.plus>precision, f2.fd.mu.plus, precision)
    df2.plus <- (f2.fd.mu.plus - f2.plus) / (fd.prec)
    df2.plus <- as.vector(df2.plus)*as.vector(mu.plus)
    df2.plus <- as.vector(df2.plus)
    d2f2delta22 <- (df2.plus - df2) / fd.prec2
    
    
    sigma.plus <- sigma + fd.prec2
    f2.plus.sigma <- dWARING(i2, mu = mu, sigma = sigma.plus)
    f2.plus.sigma <- ifelse(f2.plus.sigma>precision, f2.plus.sigma, precision)
    f2.fd.sigma.plus <- dWARING(i2, mu = mu, sigma = (sigma.plus+fd.prec))
    f2.fd.sigma.plus <- ifelse(f2.fd.sigma.plus>precision, f2.fd.sigma.plus, precision)
    df2.sigma.plus <- (f2.fd.sigma.plus - f2.plus.sigma) / (fd.prec)
    d2f2sigma2 <- (df2.sigma.plus - df2.sigma) / fd.prec2
    
    
    f2.dsigma.plus <- dWARING(i2, mu = mu, sigma = sigma.plus)
    f2.dsigma.plus <- ifelse(f2.dsigma.plus>precision, f2.dsigma.plus, precision)
    f2.fd.mu.sigma.plus <- dWARING(i2, mu = (mu+fd.prec), sigma = sigma.plus)
    f2.fd.mu.sigma.plus <- ifelse(f2.fd.mu.sigma.plus>precision, f2.fd.mu.sigma.plus, precision)
    df2.dsigma.plus <- (f2.fd.mu.sigma.plus - f2.dsigma.plus) / (fd.prec)
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)*as.vector(mu)
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec2
    
    if(univariate==FALSE) {
    
      F2.plus <- pWARING(i2, mu = mu.plus, sigma = sigma)
      F2.fd.mu.plus <- pWARING(i2, mu = (mu.plus+fd.prec), sigma = sigma)
      F2.fd.mu.plus <- ifelse(F2.fd.mu.plus<(1-precision), F2.fd.mu.plus, 1-precision)
      dF2.plus <- (F2.fd.mu.plus - F2.plus) / (fd.prec)
      dF2.plus <- as.vector(dF2.plus)*as.vector(mu.plus) 
      d2F2ddelta22 <- (dF2.plus - dF2) / fd.prec2
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
      F2.sigma.plus <- pWARING(i2, mu = mu, sigma = sigma.plus)
      F2.fd.sigma.plus <- pWARING(i2, mu = mu, sigma = (sigma.plus+fd.prec))
      F2.fd.sigma.plus <- ifelse(F2.fd.sigma.plus<(1-precision), F2.fd.sigma.plus, 1-precision)
      dF2.sigma.plus <- (F2.fd.sigma.plus - F2.sigma.plus) / (fd.prec)
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) /fd.prec2
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
      F2.fd.dsigma.plus <- pWARING(i2, mu = (mu+fd.prec), sigma = sigma.plus)
      F2.fd.dsigma.plus <- ifelse(F2.fd.dsigma.plus<(1-precision), F2.fd.dsigma.plus, 1-precision)
      dF2.dsigma.plus <- (F2.fd.dsigma.plus - F2.sigma.plus) / (fd.prec)
      dF2.dsigma.plus <- as.vector(dF2.dsigma.plus)*as.vector(mu) 
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec2
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma   
    
    }
    
    
  } else if (VC$margins[2]=="YULE") {
  
    i2 <- dat[,2]
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^7, mu, 10^7)    
    f2 <- dYULE(i2, mu)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- pYULE(i2, mu)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
    
    
    f2.fd.mu <- dYULE(i2, mu = (mu+fd.prec))
    f2.fd.mu <- ifelse(f2.fd.mu>precision, f2.fd.mu, precision)
    df2 <- (f2.fd.mu - f2) / (fd.prec)
    df2 <- as.vector(df2)*as.vector(mu)
    
    
    
    if(univariate==FALSE) {
      F2.fd.mu <- pYULE(i2, mu = (mu+fd.prec))
      F2.fd.mu <- ifelse(F2.fd.mu<(1-precision), F2.fd.mu, 1-precision)
      dF2 <- (F2.fd.mu-F2) / (fd.prec)
      dF2 <- as.vector(dF2)*as.vector(mu) 
      dF22 <- dF2 - df2
  
    }
    
    
    
    # Hessian derivative components
    
    #fd.prec2 <- 10^-3
    fd.prec2 <- sqrt(eps)
    
    mu.plus <- ifelse(exp(eta2+fd.prec2)>0.0001, exp(eta2+fd.prec2), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    f2.plus <- dYULE(i2, mu = mu.plus)
    f2.plus <- ifelse(f2.plus>precision, f2.plus, precision)
    f2.fd.mu.plus <- dYULE(i2, mu = (mu.plus+fd.prec))
    f2.fd.mu.plus <- ifelse(f2.fd.mu.plus>precision, f2.fd.mu.plus, precision)
    df2.plus <- (f2.fd.mu.plus - f2.plus) / (fd.prec)
    df2.plus <- as.vector(df2.plus)*as.vector(mu.plus)
    df2.plus <- as.vector(df2.plus)
    d2f2delta22 <- (df2.plus - df2) / fd.prec2
    
    
    if(univariate==FALSE) {
    
      F2.plus <- pYULE(i2, mu = mu.plus)
      F2.fd.mu.plus <- pYULE(i2, mu = (mu.plus+fd.prec))
      F2.fd.mu.plus <- ifelse(F2.fd.mu.plus<(1-precision), F2.fd.mu.plus, 1-precision)
      dF2.plus <- (F2.fd.mu.plus - F2.plus) / (fd.prec)
      dF2.plus <- as.vector(dF2.plus)*as.vector(mu.plus) 
      d2F2ddelta22 <- (dF2.plus - dF2) / fd.prec2
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    
    }
    
    
  } else if (VC$margins[2]=="ZABB") {
  
    i2 <- dat[,2]
    
    if(univariate==TRUE) {
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
      if (is.null(VC$X3)) nu.star <- params[(VC$X2.d2 + 2)]
      if (!is.null(VC$X3)) nu.star <- etanu <- X4 %*% params[(VC$X2.d2 + VC$X3.d2 + 1):(VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]        
      if(is.null(VC$X3))  nu.star <- params[(VC$X1.d2+VC$X2.d2+2)]
      if(!is.null(VC$X3)) nu.star <- etanu <- X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)] 
    }
      
    nu <- plogis(nu.star)
    nu <- ifelse(nu<precision, precision, nu)
    nu <- ifelse(nu>(1-precision), 1-precision, nu)
    sigma <- exp(sigma.star)
    sigma <- ifelse(sigma>precision, sigma, precision)
    sigma <- ifelse(sigma<10^7, sigma, 10^7)
    mu <- plogis(eta2)
    mu <- ifelse(mu<precision, precision, mu)
    mu <- ifelse(mu>(1-precision), 1-precision, mu) 
    f2 <- dZABB(i2, bd=bd, mu=mu, sigma=sigma, nu=nu)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- pZABB(i2, bd=bd, mu=mu, sigma=sigma, nu=nu)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
    
   
    f2.fd.mu <- dZABB(i2, mu = (mu +fd.prec), sigma = sigma, nu = nu, bd=bd)
    f2.fd.mu <- ifelse(f2.fd.mu>precision, f2.fd.mu, precision)
    df2 <- (f2.fd.mu - f2) / (fd.prec)
    df2 <- as.vector(df2)*as.vector(mu * (1-mu))
    
    f2.fd.sigma <- dZABB(i2, mu = mu, sigma = (sigma+fd.prec), nu = nu, bd=bd)
    f2.fd.sigma <- ifelse(f2.fd.sigma>precision, f2.fd.sigma, precision)
    df2.sigma <- (f2.fd.sigma - f2) / (fd.prec)
    
    f2.fd.nu <- dZABB(i2, mu = mu, sigma = sigma, nu =(nu+fd.prec), bd=bd)
    f2.fd.nu <- ifelse(f2.fd.nu>precision, f2.fd.nu, precision)
    df2.nu <- (f2.fd.nu - f2) / (fd.prec)
    
    if(univariate==FALSE) {
      F2.fd.mu <- pZABB(i2, mu = (mu+fd.prec), sigma = sigma, nu = nu, bd=bd)
      F2.fd.mu <- ifelse(F2.fd.mu<(1-precision), F2.fd.mu, 1-precision)
      dF2 <- (F2.fd.mu - F2) / (fd.prec)
      dF2 <- as.vector(dF2)*as.vector(mu * (1-mu)) 
      dF22 <- dF2 - df2
    
      F2.fd.sigma <- pZABB(i2, mu = mu, sigma = (sigma+fd.prec), nu = nu, bd=bd)
      F2.fd.sigma <- ifelse(F2.fd.sigma<(1-precision), F2.fd.sigma, 1-precision)
      dF2.sigma <- (F2.fd.sigma - F2) / (fd.prec)
      dF22.sigma <- dF2.sigma - df2.sigma
    
      F2.fd.nu <- pZABB(i2, mu = mu, sigma = sigma, nu = (nu+fd.prec), bd=bd)
      F2.fd.nu <- ifelse(F2.fd.nu<(1-precision), F2.fd.nu, 1-precision)
      dF2.nu <- (F2.fd.nu - F2) / (fd.prec)
      dF22.nu <- dF2.nu - df2.nu
    }
    
    
    # Hessian derivative components
    
    #second derivative fd precision
    #fd.prec2 <- 10^-3
    fd.prec2 <- sqrt(eps)
    
    #pmfs
    
    #delta2^2
    mu.plus <- ifelse(plogis(eta2+fd.prec2)>0.0001, plogis(eta2+fd.prec2), 0.0001)
    mu.plus <- ifelse(mu.plus>(1-precision), 1-precision, mu.plus)
    f2.plus <- dZABB(i2, mu = mu.plus, sigma = sigma, nu = nu, bd=bd)
    f2.plus <- ifelse(f2.plus>precision, f2.plus, precision)
    f2.fd.mu.plus <- dZABB(i2, mu = (mu.plus+fd.prec), sigma = sigma, nu = nu, bd=bd)
    f2.fd.mu.plus <- ifelse(f2.fd.mu.plus>precision, f2.fd.mu.plus, precision)
    df2.plus <- (f2.fd.mu.plus - f2.plus) / (fd.prec)
    df2.plus <- as.vector(df2.plus)*as.vector(mu.plus * (1-mu.plus))
    df2.plus <- as.vector(df2.plus)
    d2f2delta22 <- (df2.plus - df2) / fd.prec2
    
    #sigma^2
    sigma.plus <- sigma+fd.prec2
    f2.plus.sigma <- dZABB(i2, mu = mu, sigma = sigma.plus, nu = nu, bd=bd)
    f2.plus.sigma <- ifelse(f2.plus.sigma>precision, f2.plus.sigma, precision)
    f2.fd.sigma.plus <- dZABB(i2, mu = mu, sigma = (sigma.plus+fd.prec), nu = nu, bd=bd)
    f2.fd.sigma.plus <- ifelse(f2.fd.sigma.plus>precision, f2.fd.sigma.plus, precision)
    df2.sigma.plus <- (f2.fd.sigma.plus - f2.plus.sigma) / (fd.prec)
    d2f2sigma2 <- (df2.sigma.plus - df2.sigma) / fd.prec2
    
    #delta2sigma
    f2.dsigma.plus <- dZABB(i2, mu = mu, sigma = sigma.plus, nu = nu, bd=bd)
    f2.dsigma.plus <- ifelse(f2.dsigma.plus>precision, f2.dsigma.plus, precision)
    f2.fd.mu.sigma.plus <- dZABB(i2, mu = (mu+fd.prec), sigma = sigma.plus, nu = nu, bd=bd)
    f2.fd.mu.sigma.plus <- ifelse(f2.fd.mu.sigma.plus>precision, f2.fd.mu.sigma.plus, precision)
    df2.dsigma.plus <- (f2.fd.mu.sigma.plus - f2.dsigma.plus) / (fd.prec)
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)*as.vector(mu * (1-mu))
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec2
    
    #nu^2
    nu.plus <- nu + fd.prec2
    nu.plus <- ifelse(nu.plus>=(1-precision), (1-precision), nu.plus)
    f2.plus.nu <- dZABB(i2, mu = mu, sigma = sigma, nu = nu.plus, bd=bd)
    f2.plus.nu <- ifelse(f2.plus.nu>precision, f2.plus.nu, precision)
    f2.fd.nu.plus <- dZABB(i2, mu = mu, sigma = sigma, nu = nu.plus+fd.prec, bd=bd)
    f2.fd.nu.plus <- ifelse(f2.fd.nu.plus>precision, f2.fd.nu.plus, precision)
    df2.nu.plus <- (f2.fd.nu.plus - f2.plus.nu) / (fd.prec)
    d2f2nu2 <- (df2.nu.plus - df2.nu) / fd.prec2
    
    #delta2nu
    f2.dnu.plus <- dZABB(i2, mu = mu, sigma = sigma, nu = nu.plus, bd=bd)
    f2.dnu.plus <- ifelse(f2.dnu.plus>precision, f2.dnu.plus, precision)
    f2.fd.mu.nu.plus <- dZABB(i2, mu = (mu+fd.prec), sigma = sigma, nu = nu.plus, bd=bd)
    f2.fd.mu.nu.plus <- ifelse(f2.fd.mu.nu.plus>precision, f2.fd.mu.nu.plus, precision)
    df2.dnu.plus <- (f2.fd.mu.nu.plus - f2.dnu.plus) / (fd.prec)
    df2.dnu.plus <- as.vector(df2.dnu.plus)*as.vector(mu * (1-mu))
    d2f2delta2nu <-(df2.dnu.plus - df2) / fd.prec2
    
    #nusigma
    f2.fd.sigma.nu.plus <- dZABB(i2, mu = mu, sigma = (sigma+fd.prec), nu = nu.plus, bd=bd)
    f2.fd.sigma.nu.plus <- ifelse(f2.fd.sigma.nu.plus>precision, f2.fd.sigma.nu.plus, precision)
    df2.dnu.sigma.plus <- (f2.fd.sigma.nu.plus - f2.plus.nu) / (fd.prec)
    df2.dnu.sigma.plus <- as.vector(df2.dnu.sigma.plus)
    d2f2sigmanu <- (df2.dnu.sigma.plus - df2.sigma) / fd.prec2
    
    
    if(univariate==FALSE) {
    
      #cdfs
    
      #delta2^2
      F2.plus <- pZABB(i2, mu = mu.plus, sigma = sigma, nu = nu, bd=bd)
      F2.fd.mu.plus <- pZABB(i2, mu = (mu.plus+fd.prec), sigma = sigma, nu = nu, bd=bd)
      F2.fd.mu.plus <- ifelse(F2.fd.mu.plus<(1-precision), F2.fd.mu.plus, 1-precision)
      dF2.plus <- (F2.fd.mu.plus - F2.plus) / (fd.prec)
      dF2.plus <- as.vector(dF2.plus)*as.vector(mu.plus * (1-mu.plus)) 
      d2F2ddelta22 <-(dF2.plus - dF2) / fd.prec2
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
      #sigma^2
      F2.sigma.plus <- pZABB(i2, mu = mu, sigma = sigma.plus, nu = nu, bd=bd)
      F2.fd.sigma.plus <- pZABB(i2, mu = mu, sigma = (sigma.plus+fd.prec), nu=nu, bd=bd)
      F2.fd.sigma.plus <- ifelse(F2.fd.sigma.plus<(1-precision), F2.fd.sigma.plus, 1-precision)
      dF2.sigma.plus <- (F2.fd.sigma.plus - F2.sigma.plus) / (fd.prec)
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec2
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
      #delta2sigma
      F2.fd.dsigma.plus <- pZABB(i2, mu = (mu+fd.prec), sigma = sigma.plus, nu, bd=bd)
      F2.fd.dsigma.plus <- ifelse(F2.fd.dsigma.plus<(1-precision), F2.fd.dsigma.plus, 1-precision)
      dF2.dsigma.plus <- (F2.fd.dsigma.plus - F2.sigma.plus) / (fd.prec)
      dF2.dsigma.plus <- as.vector(dF2.dsigma.plus)*as.vector(mu * (1-mu)) 
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec2
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma
    
      #nu^2
      F2.nu.plus<- pZABB(i2, mu = mu, sigma = sigma, nu = nu.plus, bd=bd)
      F2.fd.nu.plus <- pZABB(i2, mu = mu, sigma = sigma, nu = nu.plus+fd.prec, bd=bd)
      F2.fd.nu.plus <- ifelse(F2.fd.nu.plus<(1-precision), F2.fd.nu.plus, 1-precision)
      dF2.nu.plus <- (F2.fd.nu.plus - F2.nu.plus) / (fd.prec)
      d2F2dnu2 <- (dF2.nu.plus - dF2.nu) / fd.prec2
      d2F22dnu2 <- d2F2dnu2 - d2f2nu2
    
      #delta2nu
      F2.fd.dnu.plus <- pZABB(i2, mu = (mu+fd.prec), sigma = sigma, nu.plus, bd=bd)
      F2.fd.dnu.plus <- ifelse(F2.fd.dnu.plus<(1-precision), F2.fd.dnu.plus, 1-precision)
      dF2.dnu.plus <- (F2.fd.dnu.plus - F2.nu.plus) / (fd.prec)
      dF2.dnu.plus <- as.vector(dF2.dnu.plus)*as.vector(mu * (1-mu)) 
      d2F2ddelta2dnu <- (dF2.dnu.plus - dF2) / fd.prec2
      d2F22ddelta2dnu <- d2F2ddelta2dnu - d2f2delta2nu
    
      #sigmanu
      F2.fd.dnu.sigma.plus <- pZABB(i2, mu = mu, sigma = (sigma+fd.prec), nu.plus, bd=bd)
      F2.fd.dnu.sigma.plus <- ifelse(F2.fd.dnu.sigma.plus<(1-precision), F2.fd.dnu.sigma.plus, 1-precision)
      dF2.dnu.sigma.plus <- (F2.fd.dnu.sigma.plus - F2.nu.plus) / (fd.prec)
      dF2.dnu.sigma.plus <- as.vector(dF2.dnu.sigma.plus)
      d2F2dsigmadnu <- (dF2.dnu.sigma.plus - dF2.sigma) / fd.prec2
      d2F22dsigmadnu <- d2F2dsigmadnu - d2f2sigmanu
     
    }
 
  } else if (VC$margins[2]=="ZABI") {
  
    i2 <- dat[,2]
    
    if(univariate==TRUE) {
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
    }
      
    sigma <- plogis(sigma.star)
    sigma <- ifelse(sigma<(1-precision), sigma, 1-sigma)
    sigma <- ifelse(sigma>precision, sigma, precision)
    mu <- plogis(eta2)
    mu <- ifelse(mu<(1-precision), mu, 1-mu)
    mu <- ifelse(mu>precision, mu, precision)
    f2 <- dZABI(i2, sigma = sigma,  mu = mu, bd=bd)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- pZABI(i2, sigma = sigma, mu = mu, bd=bd)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision) 
    }
    
    df2<- ifelse(i2>0, ((1 - mu)^(-1 + bd - i2) * mu^(-1 + i2) * (i2 * (-1 + (1 - mu)^bd) + bd * mu) * 
                        (-1 + sigma) * gamma(1 + bd))/((-1 + (1 - mu)^bd)^2 * gamma(1 + bd - i2) * gamma(1 + i2)), 0)
    df2 <- as.vector(df2)*(1-mu)*mu
    
    df2.sigma <- ifelse(i2>0, -(gamma(bd + 1) * mu^(i2) * (1 - mu)^(bd - i2)/((1 - (1 - mu)^(bd)) * gamma(i2 + 1) * gamma(bd - i2 + 1))), 1)
    
    
    if(univariate==FALSE) {
    
      df2.mu.int.ZABI<-function(y, mu, sigma, n) {
        
        df2<- ifelse(y>0, ((1 - mu)^(-1 + n - y) * mu^(-1 + y) * (y * (-1 + (1 - mu)^n) + n * mu) * 
                    (-1 + sigma) * gamma(1 + n))/((-1 + (1 - mu)^n)^2 * 
                    gamma(1 + n - y) * gamma(1 + y)), 0)
        df2 <- as.vector(df2)*(1-mu)*mu
        df2
      }
    
      df2.sigma.int.ZABI<-function(y, mu, sigma, n) {
        
        df2<- ifelse(y>0, -(gamma(n + 1) * mu^(y) * (1 - mu)^(n - y)/((1 - (1 - mu)^(n)) * gamma(y + 1) * gamma(n - y + 1))), 1)
        df2
      }
    
      
      ly <- length(i2)
      cdf.f <- rep(0, ly)
      cdf.f.sigma <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f[i] <- sum(df2.mu.int.ZABI(allval, mm, ms, bd))
        cdf.f.sigma[i] <- sum(df2.sigma.int.ZABI(allval, mm, ms, bd))
      }  
    
      dF2<-cdf.f
      dF22<-dF2-df2
    
      dF2.sigma <- cdf.f.sigma
      dF22.sigma <- dF2.sigma-df2.sigma
    
    
    }
    
    # Hessian derivative components
    
    mu.plus <- ifelse(plogis(eta2+fd.prec)>0.0001, plogis(eta2+fd.prec), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    df2.plus <- ifelse(i2>0, ((1 - mu.plus)^(-1 + bd - i2) * mu.plus^(-1 + i2) * (i2 * (-1 + 
                    (1 - mu.plus)^bd) + bd * mu.plus) * (-1 + sigma) * gamma(1 + bd))/((-1 + (1 - mu.plus)^bd)^2 * 
                    gamma(1 + bd - i2) * gamma(1 + i2)), 0)
    df2.plus <- as.vector(df2.plus)*(1-mu.plus)*mu.plus
    d2f2delta22 <- (df2.plus - df2) / fd.prec
    
    sigma.plus <- sigma + fd.prec
    d2f2sigma2 <- 0
    
    df2.dsigma.plus <- ifelse(i2>0, ((1 - mu)^(-1 + bd - i2) * mu^(-1 + i2) * (i2 * (-1 + (1 - mu)^bd) + bd * mu) * 
                    (-1 + sigma.plus) * gamma(1 + bd))/((-1 + (1 - mu)^bd)^2 * 
                    gamma(1 + bd - i2) * gamma(1 + i2)), 0)
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)*(1-mu)*mu
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec
    
    
    if(univariate==FALSE) {
    
      cdf.f.plus <- rep(0, ly)
      nmu <- rep(mu.plus, length = ly) 
      nsigma <- rep(sigma, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f.plus[i] <- sum(df2.mu.int.ZABI(allval, mm, ms, bd))
      }  
    
      dF2.mu.plus <- cdf.f.plus
      d2F2ddelta22 <- (dF2.mu.plus - dF2) / fd.prec
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    
      cdf.f.sigma.plus <- rep(0, ly)
      cdf.f.dsigma.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma.plus, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f.sigma.plus[i] <- sum(df2.sigma.int.ZABI(allval, mm, ms, bd))
        cdf.f.dsigma.plus[i] <- sum(df2.mu.int.ZABI(allval, mm, ms, bd))
      }  
    
      dF2.sigma.plus <- cdf.f.sigma.plus
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
      dF2.dsigma.plus <- cdf.f.dsigma.plus
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma   
    
    
    }
    
  
  } else if (VC$margins[2]=="ZALG") {
    
    i2 <- dat[,2]
    
    if(univariate==TRUE) {
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
    }
      
    sigma <- plogis(sigma.star)
    sigma <- ifelse(sigma<(1-precision), sigma, 1-sigma)
    sigma <- ifelse(sigma>precision, sigma, precision)
    mu <- plogis(eta2)
    mu <- ifelse(mu<(1-precision), mu, 1-mu)
    mu <- ifelse(mu>precision, mu, precision)
    f2 <- dZALG(i2, sigma = sigma,  mu = mu)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- pZALG(i2, sigma = sigma, mu = mu)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
    
    df2<- ifelse(i2>0, (mu^(-1 + i2) * (-1 + sigma) * (-mu + i2 * (-1 + mu) * log(1 - mu)))/(i2 * (-1 + mu) * log(1 - mu)^2), 0)
    df2 <- as.vector(df2)*(1-mu)*mu
    
    df2.sigma <- ifelse(i2>0, -((-(log(1 - mu))^(-1)) * mu^(i2)/i2), 1)
    
    
    if(univariate==FALSE) {
    
      df2.mu.int.ZALG<-function(y, mu, sigma) {
        
        df2<- ifelse(y>0, (mu^(-1 + y) * (-1 + sigma) * (-mu + y * (-1 + mu) * log(1 - mu)))/(y * (-1 + mu) * log(1 - mu)^2), 0)
        df2 <- as.vector(df2)*(1-mu)*mu
        df2
      }
    
      df2.sigma.int.ZALG<-function(y, mu, sigma) {
        
        df2<- ifelse(y>0, -((-(log(1 - mu))^(-1)) * mu^(y)/y), 1)
        df2
      }
    
      ly <- length(i2)
      cdf.f <- rep(0, ly)
      cdf.f.sigma <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f[i] <- sum(df2.mu.int.ZALG(allval, mm, ms))
        cdf.f.sigma[i] <- sum(df2.sigma.int.ZALG(allval, mm, ms))
      }  
    
      dF2<-cdf.f
      dF22<-dF2-df2
    
      dF2.sigma <- cdf.f.sigma
      dF22.sigma <- dF2.sigma-df2.sigma
  
    
    }
    
    # Hessian derivative components
    
    mu.plus <- ifelse(plogis(eta2+fd.prec)>0.0001, plogis(eta2+fd.prec), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    df2.plus <- ifelse(i2>0, (mu.plus^(-1 + i2) * (-1 + sigma) * (-mu.plus + i2 * (-1 + mu.plus) * log(1 - mu.plus)))/(i2 * 
                      (-1 + mu.plus) * log(1 - mu.plus)^2), 0)
    df2.plus <- as.vector(df2.plus)*(1-mu.plus)*mu.plus
    d2f2delta22 <- (df2.plus - df2) / fd.prec
    
    sigma.plus <- sigma + fd.prec
    d2f2sigma2 <- 0
    
    df2.dsigma.plus <- ifelse(i2>0, (mu^(-1 + i2) * (-1 + sigma.plus) * (-mu + i2 * (-1 + mu) * log(1 - mu)))/(i2 * 
                              (-1 + mu) * log(1 - mu)^2), 0)
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)*(1-mu)*mu
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec
    
    if(univariate==FALSE) {
    
      cdf.f.plus <- rep(0, ly)
      nmu <- rep(mu.plus, length = ly) 
      nsigma <- rep(sigma, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f.plus[i] <- sum(df2.mu.int.ZALG(allval, mm, ms))
      }  
    
      dF2.mu.plus <- cdf.f.plus
      d2F2ddelta22 <- (dF2.mu.plus - dF2) / fd.prec
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    
      cdf.f.sigma.plus <- rep(0, ly)
      cdf.f.dsigma.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma.plus, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f.sigma.plus[i] <- sum(df2.sigma.int.ZALG(allval, mm, ms))
        cdf.f.dsigma.plus[i] <- sum(df2.mu.int.ZALG(allval, mm, ms))
      }  
    
      dF2.sigma.plus <- cdf.f.sigma.plus
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
      dF2.dsigma.plus <- cdf.f.dsigma.plus
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma   
      
    
    }
    
    
    
  } else if (VC$margins[2]=="ZANBI") {

    i2 <- dat[,2]
    
    if(univariate==TRUE) {
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
      if (is.null(VC$X3)) nu.star <- params[(VC$X2.d2 + 2)]
      if (!is.null(VC$X3)) nu.star <- etanu <- X4 %*% params[(VC$X2.d2 + VC$X3.d2 + 1):(VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]        
      if(is.null(VC$X3))  nu.star <- params[(VC$X1.d2+VC$X2.d2+2)]
      if(!is.null(VC$X3)) nu.star <- etanu <- X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)] 
    }
      
    nu <- plogis(nu.star)
    nu <- ifelse(nu<precision, precision, nu)
    nu <- ifelse(nu>(1-precision), 1-precision, nu)
    sigma <- exp(sigma.star)
    sigma <- ifelse(sigma>precision, sigma, precision)
    sigma <- ifelse(sigma<10^7, sigma, 10^7)
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^7, mu, 10^7)  
    f2 <- dZANBI(i2, mu=mu, sigma=sigma, nu=nu)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- pZANBI(i2, mu=mu, sigma=sigma, nu=nu)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
    
  
    df2 <- ifelse(i2>0, ((-1 + nu) * (1/(1 + mu * sigma))^(1 + 1/sigma) * ((mu * sigma)/(1 + mu * sigma))^i2 * (mu + 
                  i2 * (-1 + (1/(1 + mu * sigma))^(1/sigma))) * gamma(i2 + 1/sigma))/(mu *
                  (-1 + (1/(1 + mu * sigma))^(1/sigma))^2 * gamma(1 + i2) * gamma(1/sigma)), 0)
    df2 <- as.vector(df2) * mu
    
    
    df2.sigma <- ifelse(i2>0, -(((-1 + nu) * (1/(1 + mu * sigma))^(1 + 1/sigma) * ((mu * sigma)/(1 + mu * sigma))^i2 * 
                        gamma(i2 + 1/sigma) * (i2 * sigma - mu * sigma - 
                        i2 * sigma * (1/(1 + mu * sigma))^(1/sigma) - log(1/(1 + mu * sigma)) -
                        mu * sigma * log(1/(1 + mu * sigma)) + (1 + mu * sigma) * (-1 + (1/(1 + mu * sigma))^(1/
                        sigma)) * digamma(i2 + 1/sigma) - (1 + mu * sigma) * (-1 + (1/(1 + mu * sigma))^(1/
                        sigma)) * digamma(1/sigma)))/(sigma^2 * (-1 + (1/(1 + mu * sigma))^(1/sigma))^2 * gamma(1 + i2) * 
                        gamma(1/sigma))), 0)
    
    
    df2.nu <- ifelse(i2>0, (1/(1 + mu * sigma))^(1/sigma) * (mu * sigma/(1 + mu * sigma))^i2 * 
                    gamma(1/sigma + i2)/((-1 + (1/(1 + sigma * mu))^(1/sigma)) * 
                    gamma(1/sigma) * gamma(1 + i2)), 1)
    
    
    if(univariate==FALSE) {
    
      df2.mu.int.ZANBI <- function(y, mu, sigma, nu) {
        
        df2 <- ifelse(y>0, ((-1 + nu) * (1/(1 + mu * sigma))^(1 + 1/sigma) * ((mu * sigma)/(1 + mu * sigma))^y * (mu + 
                          y * (-1 + (1/(1 + mu * sigma))^(1/sigma))) * gamma(y + 1/sigma))/(mu *
                          (-1 + (1/(1 + mu * sigma))^(1/sigma))^2 * gamma(1 + y) * gamma(1/sigma)), 0)
        df2 <- as.vector(df2) * mu
      
      }
    
      df2.sigma.int.ZANBI <-function(y, mu, sigma, nu) {
     
        df2 <- ifelse(y>0, -(((-1 + nu) * (1/(1 + mu * sigma))^(1 + 1/sigma) * ((mu * sigma)/(1 + mu * sigma))^y * 
                        gamma(y + 1/sigma) * (y * sigma - mu * sigma - 
                        y * sigma * (1/(1 + mu * sigma))^(1/sigma) - log(1/(1 + mu * sigma)) -
                        mu * sigma * log(1/(1 + mu * sigma)) + (1 + mu * sigma) * (-1 + (1/(1 + mu * sigma))^(1/
                        sigma)) * digamma(y + 1/sigma) - (1 + mu * sigma) * (-1 + (1/(1 + mu * sigma))^(1/
                        sigma)) * digamma(1/sigma)))/(sigma^2 * (-1 + (1/(1 + mu * sigma))^(1/sigma))^2 * gamma(1 + y) * 
                        gamma(1/sigma))), 0)
        df2
      
      }
    
      df2.nu.int.ZANBI <- function(y, mu, sigma, nu) {
    
        df2  <- ifelse(y>0, (1/(1 + mu * sigma))^(1/sigma) * (mu * sigma/(1 + mu * sigma))^y * 
                      gamma(1/sigma + y)/((-1 + (1/(1 + sigma * mu))^(1/sigma)) * 
                      gamma(1/sigma) * gamma(1 + y)), 1) 
        df2
      }
    
    
      ly <- length(i2)
      cdf.f <- rep(0, ly)
      cdf.f.sigma <- rep(0, ly)
      cdf.f.nu <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
      nnu <- rep(nu, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        cdf.f[i] <- sum(df2.mu.int.ZANBI(allval, mm, ms, mn))
        cdf.f.sigma[i] <- sum(df2.sigma.int.ZANBI(allval, mm, ms, mn))
        cdf.f.nu[i] <- sum(df2.nu.int.ZANBI(allval, mm, ms, mn))
      }  
    
      dF2<-cdf.f
      dF22<-dF2-df2
    
      dF2.sigma <- cdf.f.sigma
      dF22.sigma <- dF2.sigma-df2.sigma
    
      dF2.nu <- cdf.f.nu
      dF22.nu <- dF2.nu-df2.nu
    
    }
    
    
    # Hessian derivative components
    
    #fd.prec2 <- 10^-3
    fd.prec2 <- sqrt(eps)
    
    #delta2delta2
    mu.plus <- ifelse(exp(eta2+fd.prec)>0.0001, exp(eta2+fd.prec), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    
    df2.plus <- ifelse(i2>0, ((-1 + nu) * (1/(1 + mu.plus * sigma))^(1 + 1/sigma) * ((mu.plus * sigma)/(1 + mu.plus * sigma))^i2 * 
                      (mu.plus + i2 * (-1 + (1/(1 + mu.plus * sigma))^(1/sigma))) * gamma(i2 + 1/sigma))/(mu.plus *
                      (-1 + (1/(1 + mu.plus * sigma))^(1/sigma))^2 * gamma(1 + i2) * gamma(1/sigma)), 0)
    df2.plus <- as.vector(df2.plus) * mu.plus
    d2f2delta22 <- (df2.plus - df2) / fd.prec
    
    #sigma^2
    sigma.plus <- sigma + fd.prec
    df2.sigma.plus <- ifelse(i2>0, -(((-1 + nu) * (1/(1 + mu * sigma.plus))^(1 + 1/sigma.plus) * ((mu * sigma.plus)/(1 + mu * 
                        sigma.plus))^i2 * gamma(i2 + 1/sigma.plus) * (i2 * sigma.plus - mu * sigma.plus - 
                        i2 * sigma.plus * (1/(1 + mu * sigma.plus))^(1/sigma.plus) - log(1/(1 + mu * sigma.plus)) -
                        mu * sigma.plus * log(1/(1 + mu * sigma.plus)) + (1 + mu * sigma.plus) * (-1 + (1/(1 + mu * sigma.plus))^(1/
                        sigma.plus)) * digamma(i2 + 1/sigma.plus) - (1 + mu * sigma.plus) * (-1 + (1/(1 + mu * sigma.plus))^(1/
                        sigma.plus)) * digamma(1/sigma.plus)))/(sigma.plus^2 * (-1 + (1/(1 + mu * sigma.plus))^(1/sigma.plus))^2 * 
                        gamma(1 + i2) * gamma(1/sigma.plus))), 0)
    d2f2sigma2 <- (df2.sigma.plus - df2.sigma) / fd.prec
    
    #delta2sigma
    
    df2.dsigma.plus <- ifelse(i2>0, ((-1 + nu) * (1/(1 + mu * sigma.plus))^(1 + 1/sigma.plus) * 
                              ((mu * sigma.plus)/(1 + mu * sigma.plus))^i2 * (mu + 
                              i2 * (-1 + (1/(1 + mu * sigma.plus))^(1/sigma.plus))) * gamma(i2 + 1/sigma.plus))/(mu *
                              (-1 + (1/(1 + mu * sigma.plus))^(1/sigma.plus))^2 * gamma(1 + i2) * gamma(1/sigma.plus)), 0)
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)*mu
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec
    
    
    #nu^2
    nu.plus <- nu + fd.prec
    nu.plus <- ifelse(nu.plus>=(1-precision), (1-precision), nu.plus)
    d2f2nu2 <- 0
    
    
    #delta2nu
    df2.dnu.plus <- ifelse(i2>0, ((-1 + nu.plus) * (1/(1 + mu * sigma))^(1 + 1/sigma) * ((mu * sigma)/(1 + mu * sigma))^i2 *
                          (mu + i2 * (-1 + (1/(1 + mu * sigma))^(1/sigma))) * gamma(i2 + 1/sigma))/(mu *
                          (-1 + (1/(1 + mu * sigma))^(1/sigma))^2 * gamma(1 + i2) * gamma(1/sigma)), 0)
    df2.dnu.plus <- as.vector(df2.dnu.plus)*mu
    d2f2delta2nu <- (df2.dnu.plus - df2) / fd.prec
    
    
    #nusigma
    df2.nu.sigma.plus <- ifelse(i2>0, (1/(1 + mu * sigma.plus))^(1/sigma.plus) * (mu * sigma.plus/(1 + mu * sigma.plus))^i2 * 
                              gamma(1/sigma.plus + i2)/((-1 + (1/(1 + sigma.plus * mu))^(1/sigma.plus)) * 
                              gamma(1/sigma.plus) * gamma(1 + i2)), 1)
    d2f2sigmanu <- (df2.nu.sigma.plus - df2.nu) / fd.prec
    
    
    if(univariate==FALSE) {
    
      cdf.f.plus <- rep(0, ly)
      nmu <- rep(mu.plus, length = ly) 
      nsigma <- rep(sigma, length = ly)
      nnu <- rep(nu, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        pdfall <- df2.mu.int.ZANBI(allval, mm, ms, mn)
        cdf.f.plus[i] <- sum(pdfall)
      }  
    
      dF2.mu.plus <- cdf.f.plus
      d2F2ddelta22 <- (dF2.mu.plus - dF2) / fd.prec
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    
      cdf.f.sigma.plus <- rep(0, ly)
      cdf.f.dsigma.plus <- rep(0, ly)
      cdf.f.dnu.dsigma.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma.plus, length = ly)
      nnu <- rep(nu, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        cdf.f.sigma.plus[i] <- sum(df2.sigma.int.ZANBI(allval, mm, ms, mn))
        cdf.f.dsigma.plus[i] <- sum(df2.mu.int.ZANBI(allval, mm, ms, mn))
        cdf.f.dnu.dsigma.plus[i] <- sum(df2.nu.int.ZANBI(allval, mm, ms, mn))
      }  
    
      dF2.sigma.plus <- cdf.f.sigma.plus
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
      dF2.dsigma.plus <- cdf.f.dsigma.plus
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma   
    
      dF2.dnu.dsigma.plus <- cdf.f.dnu.dsigma.plus
      d2F2dsigmadnu <- (dF2.dnu.dsigma.plus - dF2.nu) / fd.prec
      d2F22dsigmadnu <- d2F2dsigmadnu - d2f2sigmanu
      
      

      cdf.f.nu.plus <- rep(0, ly)
      cdf.f.dnu.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
      nnu <- rep(nu.plus, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        cdf.f.nu.plus[i] <- sum(df2.nu.int.ZANBI(allval, mm, ms, mn))
        cdf.f.dnu.plus[i] <- sum(df2.mu.int.ZANBI(allval, mm, ms, mn))
      } 
    
      dF2.nu.plus <- cdf.f.nu.plus
      d2F2dnu2 <- (dF2.nu.plus - dF2.nu) / fd.prec
      d2F22dnu2 <- d2F2dnu2 - d2f2nu2
    
      dF2.dnu.plus <- cdf.f.dnu.plus
      d2F2ddelta2dnu <- (dF2.dnu.plus - dF2) / fd.prec
      d2F22ddelta2dnu <- d2F2ddelta2dnu - d2f2delta2nu
    
    }
    
    

  } else if (VC$margins[2]=="ZAP") {
    
    i2 <- dat[,2]
    
    if(univariate==TRUE) {
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
    }
      
    sigma <- plogis(sigma.star)
    sigma <- ifelse(sigma<(1-precision), sigma, 1-sigma)
    sigma <- ifelse(sigma>precision, sigma, precision)
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^7, mu, 10^7)
    f2 <- dZAP(i2, sigma = sigma,  mu = mu)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- pZAP(i2, sigma = sigma, mu = mu)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision) 
    }
    
    
    df2 <- ifelse(i2>0, (mu^(-1 + i2) * (i2 - exp(mu) * i2 + exp(mu)* mu) * (-1 + sigma))/((-1 + exp(mu))^2 * gamma(1 + i2)), 0)
    df2 <- as.vector(df2)*mu
    
    df2.sigma <- ifelse(i2>0, -(exp(-mu) * mu^i2/(gamma(i2 + 1) * (1 - exp(-mu)))), 1)
    
    
    if(univariate==FALSE) {
    
      df2.mu.int.ZAP<-function(y, mu, sigma) {
        
        df2 <- ifelse(y>0, (mu^(-1 + y) * (y - exp(mu) * y + exp(mu)* mu) * (-1 + sigma))/((-1 + exp(mu))^2 * gamma(1 + y)), 0)
        df2 <- as.vector(df2)*mu
        df2
      }
    
      df2.sigma.int.ZAP<-function(y, mu, sigma) {
        df2 <- ifelse(y>0, -(exp(-mu) * mu^y/(gamma(y + 1) * (1 - exp(-mu)))), 1) 
        df2
      }
    
    
      ly <- length(i2)
      cdf.f <- rep(0, ly)
      cdf.f.sigma <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f[i] <- sum(df2.mu.int.ZAP(allval, mm, ms))
        cdf.f.sigma[i] <- sum(df2.sigma.int.ZAP(allval, mm, ms))
      }  
    
      dF2<-cdf.f
      dF22<-dF2-df2
    
      dF2.sigma <- cdf.f.sigma
      dF22.sigma <- dF2.sigma-df2.sigma
    
    }
    
    # Hessian derivative components
    
    
    mu.plus <- ifelse(exp(eta2+fd.prec)>0.0001, exp(eta2+fd.prec), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    
    df2.plus <- ifelse(i2>0, (mu.plus^(-1 + i2) * (i2 - exp(mu.plus) * i2 + exp(mu.plus)* mu.plus) * 
                        (-1 + sigma))/((-1 + exp(mu.plus))^2 * gamma(1 + i2)), 0)
    df2.plus <- as.vector(df2.plus)*mu.plus
    d2f2delta22 <- (df2.plus - df2) / fd.prec
    
    
    sigma.plus <- sigma+fd.prec
    d2f2sigma2 <- 0
    
    
    df2.dsigma.plus <- ifelse(i2>0, (mu^(-1 + i2) * (i2 - exp(mu) * i2 + exp(mu)* mu) * 
                              (-1 + sigma.plus))/((-1 + exp(mu))^2 * gamma(1 + i2)), 0)
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)*mu
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec
    
    
    if(univariate==FALSE) {
    
      cdf.f.plus <- rep(0, ly)
      nmu <- rep(mu.plus, length = ly) 
      nsigma <- rep(sigma, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f.plus[i] <- sum(df2.mu.int.ZAP(allval, mm, ms))
      }  
    
      dF2.mu.plus <- cdf.f.plus
      d2F2ddelta22 <- (dF2.mu.plus - dF2) / fd.prec
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    
      cdf.f.sigma.plus <- rep(0, ly)
      cdf.f.dsigma.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma.plus, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f.sigma.plus[i] <- sum(df2.sigma.int.ZAP(allval, mm, ms))
        cdf.f.dsigma.plus[i] <- sum(df2.mu.int.ZAP(allval, mm, ms))
      }  
    
      dF2.sigma.plus <- cdf.f.sigma.plus
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
      dF2.dsigma.plus <- cdf.f.dsigma.plus
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma   
         
       
    
    }

  } else if (VC$margins[2]=="ZIBB") {
    
    i2 <- dat[,2]
    
    if(univariate==TRUE) {
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
      if (is.null(VC$X3)) nu.star <- params[(VC$X2.d2 + 2)]
      if (!is.null(VC$X3)) nu.star <- etanu <- X4 %*% params[(VC$X2.d2 + VC$X3.d2 + 1):(VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]        
      if(is.null(VC$X3))  nu.star <- params[(VC$X1.d2+VC$X2.d2+2)]
      if(!is.null(VC$X3)) nu.star <- etanu <- X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)] 
    }
      
    nu <- plogis(nu.star)
    nu <- ifelse(nu<precision, precision, nu)
    nu <- ifelse(nu>(1-precision), 1-precision, nu)
    sigma <- exp(sigma.star)
    sigma <- ifelse(sigma>precision, sigma, precision)
    sigma <- ifelse(sigma<10^7, sigma, 10^7)
    mu <- plogis(eta2)
    mu <- ifelse(mu<precision, precision, mu)
    mu <- ifelse(mu>(1-precision), 1-precision, mu)
    f2 <- dZIBB(i2, bd=bd, mu=mu, sigma=sigma, nu=nu)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- pZIBB(i2, bd=bd, mu=mu, sigma=sigma, nu=nu)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
   
    df2 <- ifelse(i2>0, ((-1 + nu)*gamma(1 + bd)*gamma(i2 + mu/sigma)*gamma(1/sigma)*
                  gamma(bd - (-1 + mu + i2*sigma)/sigma)*(-digamma(i2 + mu/sigma) - 
                  digamma((1 - mu)/sigma) + digamma(mu/sigma) + digamma(bd - (-1 + mu + i2*sigma)/sigma)))/
                  (sigma*gamma(1 + bd - i2)*gamma(1 + i2)*gamma(bd + 1/sigma)*gamma((1 - mu)/sigma)*gamma(mu/sigma)), 
                  ((-1 + nu)*gamma(bd + (1 - mu)/sigma)*gamma(1/sigma)*(digamma(bd + (1 - mu)/sigma) - 
                  digamma((1 - mu)/sigma)))/(sigma*gamma(bd + 1/sigma)*gamma((1 - mu)/sigma)))
    df2 <- as.vector(df2) * mu * (1-mu)
    
    
    df2.sigma <- ifelse(i2>0, -(((-1 + nu)*gamma(1 + bd)*gamma(i2 + mu/sigma)*gamma(1/sigma)*
                      gamma(bd - (-1 + mu + i2*sigma)/sigma)*(digamma(bd + 1/sigma) - 
                      mu*digamma(i2 + mu/sigma) - digamma(1/sigma) + 
                      digamma((1 - mu)/sigma) - mu*digamma((1 - mu)/sigma) + 
                      mu*digamma(mu/sigma) - digamma(bd - (-1 + mu + i2*sigma)/sigma) + 
                      mu*digamma(bd - (-1 + mu + i2*sigma)/sigma)))/
                      (sigma^2*gamma(1 + bd - i2)*gamma(1 + i2)*gamma(bd + 1/sigma)*
                      gamma((1 - mu)/sigma)*gamma(mu/sigma))), 
                      -(((-1 + nu)*gamma(bd + (1 - mu)/sigma)*gamma(1/sigma)*
                      (digamma(bd + 1/sigma) + (-1 + mu)*digamma(bd + (1 - mu)/sigma) - 
                      digamma(1/sigma) + digamma((1 - mu)/sigma) - 
                      mu*digamma((1 - mu)/sigma)))/(sigma^2*gamma(bd + 1/sigma)*
                      gamma((1 - mu)/sigma))))
    
    
    df2.nu <- ifelse(i2>0, -(gamma(bd + 1)/(gamma(i2 + 1) * gamma(bd - i2 + 1)) * (gamma(1/sigma) * 
                  gamma(i2 + mu/sigma) * gamma(bd + (1 - mu)/sigma - i2))/(gamma(bd + 
                  1/sigma) * gamma(mu/sigma) * gamma((1 - mu)/sigma))), 
                  1 - (gamma(bd + 1)/(gamma(1) * gamma(bd + 1)) * (gamma(1/sigma) * 
                  gamma(mu/sigma) * gamma(bd + (1 - mu)/sigma))/(gamma(bd + 1/sigma) * 
                  gamma(mu/sigma) * gamma((1 - mu)/sigma))))
    
    
    if(univariate==FALSE) {
    
      df2.mu.int.ZIBB <- function(y, mu, sigma, nu, n) {
        
        df2 <- ifelse(y>0, ((-1 + nu)*gamma(1 + n)*gamma(y + mu/sigma)*gamma(1/sigma)*
                  gamma(n - (-1 + mu + y*sigma)/sigma)*(-digamma(y + mu/sigma) - 
                  digamma((1 - mu)/sigma) + digamma(mu/sigma) + digamma(n - (-1 + mu + y*sigma)/sigma)))/
                  (sigma*gamma(1 + n - y)*gamma(1 + y)*gamma(n + 1/sigma)*gamma((1 - mu)/sigma)*gamma(mu/sigma)), 
                ((-1 + nu)*gamma(n + (1 - mu)/sigma)*gamma(1/sigma)*(digamma(n + (1 - mu)/sigma) - 
                digamma((1 - mu)/sigma)))/(sigma*gamma(n + 1/sigma)*gamma((1 - mu)/sigma)))
        df2 <- as.vector(df2) * mu*(1-mu)
      
      }
    
      df2.sigma.int.ZIBB<-function(y, mu, sigma, nu, n) {
     
        df2 <- ifelse(y>0, -(((-1 + nu)*gamma(1 + n)*gamma(y + mu/sigma)*gamma(1/sigma)*
                  gamma(n - (-1 + mu + y*sigma)/sigma)*(digamma(n + 1/sigma) - 
                  mu*digamma(y + mu/sigma) - digamma(1/sigma) + 
                  digamma((1 - mu)/sigma) - mu*digamma((1 - mu)/sigma) + 
                  mu*digamma(mu/sigma) - digamma(n - (-1 + mu + y*sigma)/sigma) + 
                  mu*digamma(n - (-1 + mu + y*sigma)/sigma)))/
                  (sigma^2*gamma(1 + n - y)*gamma(1 + y)*gamma(n + 1/sigma)*
                  gamma((1 - mu)/sigma)*gamma(mu/sigma))), 
                  -(((-1 + nu)*gamma(n + (1 - mu)/sigma)*gamma(1/sigma)*
                  (digamma(n + 1/sigma) + (-1 + mu)*digamma(n + (1 - mu)/sigma) - 
                  digamma(1/sigma) + digamma((1 - mu)/sigma) - 
                  mu*digamma((1 - mu)/sigma)))/(sigma^2*gamma(n + 1/sigma)*
                  gamma((1 - mu)/sigma))))
        df2
      
      }
    
      df2.nu.int.ZIBB <- function(y, mu, sigma, nu, n) {
    
        df2  <- ifelse(y>0, -(gamma(n + 1)/(gamma(y + 1) * gamma(n - y + 1)) * (gamma(1/sigma) * 
                gamma(y + mu/sigma) * gamma(n + (1 - mu)/sigma - y))/(gamma(n + 
                1/sigma) * gamma(mu/sigma) * gamma((1 - mu)/sigma))), 1 - (gamma(n + 1)/(gamma(1) * gamma(n + 1)) * (gamma(1/sigma) * 
                gamma(mu/sigma) * gamma(n + (1 - mu)/sigma))/(gamma(n + 1/sigma) * 
                gamma(mu/sigma) * gamma((1 - mu)/sigma)))) 
        df2
      }
    
    
      ly <- length(i2)
      cdf.f <- rep(0, ly)
      cdf.f.sigma <- rep(0, ly)
      cdf.f.nu <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
      nnu <- rep(nu, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        cdf.f[i] <- sum(df2.mu.int.ZIBB(allval, mm, ms, mn, bd))
        cdf.f.sigma[i] <- sum(df2.sigma.int.ZIBB(allval, mm, ms, mn, bd))
        cdf.f.nu[i] <- sum(df2.nu.int.ZIBB(allval, mm, ms, mn, bd))
      }  
    
      dF2<-cdf.f
      dF22<-dF2-df2
    
      dF2.sigma <- cdf.f.sigma
      dF22.sigma <- dF2.sigma-df2.sigma
    
      dF2.nu <- cdf.f.nu
      dF22.nu <- dF2.nu-df2.nu
    
    }
    
    
    # Hessian derivative components
    
    
    #delta2delta2
    mu.plus <- plogis(eta2+fd.prec)
    mu.plus <- ifelse(mu.plus<precision, precision, mu.plus)
    mu.plus <- ifelse(mu.plus>(1-precision), 1-precision, mu.plus)
    
    df2.plus <- ifelse(i2>0, ((-1 + nu)*gamma(1 + bd)*gamma(i2 + mu.plus/sigma)*gamma(1/sigma)*
                      gamma(bd - (-1 + mu.plus + i2*sigma)/sigma)*(-digamma(i2 + mu.plus/sigma) - 
                      digamma((1 - mu.plus)/sigma) + digamma(mu.plus/sigma) + 
                      digamma(bd - (-1 + mu.plus + i2*sigma)/sigma)))/(sigma*gamma(1 + bd - i2)*gamma(1 + i2)*gamma(bd + 1/sigma)*
                      gamma((1 - mu.plus)/sigma)*gamma(mu.plus/sigma)), 
                      ((-1 + nu)*gamma(bd + (1 - mu.plus)/sigma)*gamma(1/sigma)*(digamma(bd + (1 - mu.plus)/sigma) - 
                      digamma((1 - mu.plus)/sigma)))/(sigma*gamma(bd + 1/sigma)*gamma((1 - mu.plus)/sigma)))
    df2.plus <- as.vector(df2.plus) * mu.plus*(1-mu.plus)
    d2f2delta22 <- (df2.plus - df2) / fd.prec
    
    #sigma^2
    sigma.plus <- sigma + fd.prec
    df2.sigma.plus <- ifelse(i2>0, -(((-1 + nu)*gamma(1 + bd)*gamma(i2 + mu/sigma.plus)*gamma(1/sigma.plus)*
                        gamma(bd - (-1 + mu + i2*sigma.plus)/sigma.plus)*(digamma(bd + 1/sigma.plus) - 
                        mu*digamma(i2 + mu/sigma.plus) - digamma(1/sigma.plus) + 
                        digamma((1 - mu)/sigma.plus) - mu*digamma((1 - mu)/sigma.plus) + 
                        mu*digamma(mu/sigma.plus) - digamma(bd - (-1 + mu + i2*sigma.plus)/sigma.plus) + 
                        mu*digamma(bd - (-1 + mu + i2*sigma.plus)/sigma.plus)))/
                        (sigma.plus^2*gamma(1 + bd - i2)*gamma(1 + i2)*gamma(bd + 1/sigma.plus)*
                        gamma((1 - mu)/sigma.plus)*gamma(mu/sigma.plus))), 
                        -(((-1 + nu)*gamma(bd + (1 - mu)/sigma.plus)*gamma(1/sigma.plus)*
                        (digamma(bd + 1/sigma.plus) + (-1 + mu)*digamma(bd + (1 - mu)/sigma.plus) - 
                        digamma(1/sigma.plus) + digamma((1 - mu)/sigma.plus) - 
                        mu*digamma((1 - mu)/sigma.plus)))/(sigma.plus^2*gamma(bd + 1/sigma.plus)*
                        gamma((1 - mu)/sigma.plus))))
    
    d2f2sigma2 <- (df2.sigma.plus - df2.sigma) / fd.prec
    
    #delta2sigma
    
    df2.dsigma.plus <- ifelse(i2>0, ((-1 + nu)*gamma(1 + bd)*gamma(i2 + mu/sigma.plus)*gamma(1/sigma.plus)*
                            gamma(bd - (-1 + mu + i2*sigma.plus)/sigma.plus)*(-digamma(i2 + mu/sigma.plus) - 
                            digamma((1 - mu)/sigma.plus) + digamma(mu/sigma.plus) + 
                            digamma(bd - (-1 + mu + i2*sigma.plus)/sigma.plus)))/
                            (sigma.plus*gamma(1 + bd - i2)*gamma(1 + i2)*gamma(bd + 1/sigma.plus)*
                            gamma((1 - mu)/sigma.plus)*gamma(mu/sigma.plus)), 
                            ((-1 + nu)*gamma(bd + (1 - mu)/sigma.plus)*gamma(1/sigma.plus)*(digamma(bd + (1 - mu)/sigma.plus) - 
                            digamma((1 - mu)/sigma.plus)))/(sigma.plus*gamma(bd + 1/sigma.plus)*gamma((1 - mu)/sigma.plus)))
    
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)*(mu*(1-mu))
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec
    
    
    #nu^2
    nu.plus <- nu + fd.prec
    nu.plus <- ifelse(nu.plus>=(1-precision), (1-precision), nu.plus)
    df2.nu.plus <- ifelse(i2>0, -(gamma(bd + 1)/(gamma(i2 + 1) * gamma(bd - i2 + 1)) * (gamma(1/sigma) * 
                gamma(i2 + mu/sigma) * gamma(bd + (1 - mu)/sigma - i2))/(gamma(bd + 
                1/sigma) * gamma(mu/sigma) * gamma((1 - mu)/sigma))), 
                1 - (gamma(bd + 1)/(gamma(1) * gamma(bd + 1)) * (gamma(1/sigma) * 
                gamma(mu/sigma) * gamma(bd + (1 - mu)/sigma))/(gamma(bd + 1/sigma) * 
                gamma(mu/sigma) * gamma((1 - mu)/sigma))))
    
    d2f2nu2 <- (df2.nu.plus - df2.nu) / fd.prec
    
    
    #delta2nu
    df2.dnu.plus <- ifelse(i2>0, ((-1 + nu.plus)*gamma(1 + bd)*gamma(i2 + mu/sigma)*gamma(1/sigma)*
              gamma(bd - (-1 + mu + i2*sigma)/sigma)*(-digamma(i2 + mu/sigma) - 
              digamma((1 - mu)/sigma) + digamma(mu/sigma) + 
              digamma(bd - (-1 + mu + i2*sigma)/sigma)))/
              (sigma*gamma(1 + bd - i2)*gamma(1 + i2)*gamma(bd + 1/sigma)*
              gamma((1 - mu)/sigma)*gamma(mu/sigma)), 
              ((-1 + nu.plus)*gamma(bd + (1 - mu)/sigma)*gamma(1/sigma)*(digamma(bd + (1 - mu)/sigma) - 
              digamma((1 - mu)/sigma)))/(sigma*gamma(bd + 1/sigma)*gamma((1 - mu)/sigma)))
    df2.dnu.plus <- as.vector(df2.dnu.plus)*(mu*(1-mu))
    d2f2delta2nu <- (df2.dnu.plus - df2) / fd.prec
    
    
    #nusigma
    df2.nu.sigma.plus <- ifelse(i2>0, -(gamma(bd + 1)/(gamma(i2 + 1) * gamma(bd - i2 + 1)) * (gamma(1/sigma.plus) * 
                              gamma(i2 + mu/sigma.plus) * gamma(bd + (1 - mu)/sigma.plus - i2))/(gamma(bd + 
                              1/sigma.plus) * gamma(mu/sigma.plus) * gamma((1 - mu)/sigma.plus))) ,
                              1 - (gamma(bd + 1)/(gamma(1) * gamma(bd + 1)) * (gamma(1/sigma.plus) * 
                              gamma(mu/sigma.plus) * gamma(bd + (1 - mu)/sigma.plus))/(gamma(bd + 1/sigma.plus) * 
                              gamma(mu/sigma.plus) * gamma((1 - mu)/sigma.plus))) )
    d2f2sigmanu <- (df2.nu.sigma.plus - df2.nu) / fd.prec
    
    
    if(univariate==FALSE) {
    
    
      cdf.f.plus <- rep(0, ly)
      nmu <- rep(mu.plus, length = ly) 
      nsigma <- rep(sigma, length = ly)
      nnu <- rep(nu, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        pdfall <- df2.mu.int.ZIBB(allval, mm, ms, mn, bd)
        cdf.f.plus[i] <- sum(pdfall)
      }  
    
      dF2.mu.plus <- cdf.f.plus
      d2F2ddelta22 <- (dF2.mu.plus - dF2) / fd.prec
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    
      cdf.f.sigma.plus <- rep(0, ly)
      cdf.f.dsigma.plus <- rep(0, ly)
      cdf.f.dnu.dsigma.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma.plus, length = ly)
      nnu <- rep(nu, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        cdf.f.sigma.plus[i] <- sum(df2.sigma.int.ZIBB(allval, mm, ms, mn, bd))
        cdf.f.dsigma.plus[i] <- sum(df2.mu.int.ZIBB(allval, mm, ms, mn, bd))
        cdf.f.dnu.dsigma.plus[i] <- sum(df2.nu.int.ZIBB(allval, mm, ms, mn, bd))
      }  
    
      dF2.sigma.plus <- cdf.f.sigma.plus
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
      dF2.dsigma.plus <- cdf.f.dsigma.plus
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma   
    
      dF2.dnu.dsigma.plus <- cdf.f.dnu.dsigma.plus
      d2F2dsigmadnu <- (dF2.dnu.dsigma.plus - dF2.nu) / fd.prec
      d2F22dsigmadnu <- d2F2dsigmadnu - d2f2sigmanu
      
      
      
      cdf.f.nu.plus <- rep(0, ly)
      cdf.f.dnu.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
      nnu <- rep(nu.plus, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        cdf.f.nu.plus[i] <- sum(df2.nu.int.ZIBB(allval, mm, ms, mn, bd))
        cdf.f.dnu.plus[i] <- sum(df2.mu.int.ZIBB(allval, mm, ms, mn, bd))
      } 
    
      dF2.nu.plus <- cdf.f.nu.plus
      d2F2dnu2 <- (dF2.nu.plus - dF2.nu) / fd.prec
      d2F22dnu2 <- d2F2dnu2 - d2f2nu2
    
      dF2.dnu.plus <- cdf.f.dnu.plus
      d2F2ddelta2dnu <- (dF2.dnu.plus - dF2) / fd.prec
      d2F22ddelta2dnu <- d2F2ddelta2dnu - d2f2delta2nu
    
    }
    
    
    
    
    
  } else if (VC$margins[2]=="ZIBI") {
    

    i2 <- dat[,2]
    
    if(univariate==TRUE) {
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
    }
      
    sigma <- plogis(sigma.star)
    sigma <- ifelse(sigma<(1-precision), sigma, 1-sigma)
    sigma <- ifelse(sigma>precision, sigma, precision)
    mu <- plogis(eta2)
    mu <- ifelse(mu<(1-precision), mu, 1-mu)
    mu <- ifelse(mu>precision, mu, precision)
    f2 <- dZIBI(i2, sigma = sigma,  mu = mu, bd=bd)
    f2 <- ifelse(f2>precision, f2, precision)

    
    if(univariate==FALSE) {
      F2 <- pZIBI(i2, sigma = sigma, mu = mu, bd=bd)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision) 
    }
    
    
    df2<- ifelse(i2>0, ((1 - mu)^(-1 + bd - i2) * mu^(-1 + i2) * (-i2 + bd * mu) * (-1 + sigma) * gamma(1 + bd))/(gamma(1 + bd - i2)
                        * gamma(1 + i2)), 
                    -((1 - sigma) * ((1 - mu)^(bd - 1) * bd)))
    df2 <- as.vector(df2)*(1-mu)*mu
    
    df2.sigma <- ifelse(i2>0,  -(gamma(bd + 1) * mu^i2 * (1 - mu)^(bd - i2)/(gamma(i2 + 1) * gamma(bd - i2 + 1))), 1 - (1 - mu)^bd)
    
    
    if(univariate==FALSE) {
    
      df2.mu.int.ZIBI<-function(y, mu, sigma, n) {
        
        df2<- ifelse(y>0, ((1 - mu)^(-1 + n - y) * mu^(-1 + y) * (-y + n * mu) * (-1 + sigma) * gamma(1 + n))/(gamma(1 + n - y) * 
                    gamma(1 + y)), 
                    -((1 - sigma) * ((1 - mu)^(n - 1) * n)))
        df2 <- as.vector(df2)*(1-mu)*mu
        df2
      }
    
      df2.sigma.int.ZIBI<-function(y, mu, sigma, n) {
        
        df2<- ifelse(y>0,  -(gamma(n + 1) * mu^y * (1 - mu)^(n - y)/(gamma(y + 1) * gamma(n - 
                    y + 1))), 1 - (1 - mu)^n)
        df2
      }
    
      ly <- length(i2)
      cdf.f <- rep(0, ly)
      cdf.f.sigma <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f[i] <- sum(df2.mu.int.ZIBI(allval, mm, ms, bd))
        cdf.f.sigma[i] <- sum(df2.sigma.int.ZIBI(allval, mm, ms, bd))
      }  
    
      dF2<-cdf.f
      dF22<-dF2-df2
    
      dF2.sigma <- cdf.f.sigma
      dF22.sigma <- dF2.sigma-df2.sigma
    
    
    }
    
    # Hessian derivative components
    
    mu.plus <- ifelse(plogis(eta2+fd.prec)>0.0001, plogis(eta2+fd.prec), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    df2.plus <- ifelse(i2>0, ((1 - mu.plus)^(-1 + bd - i2) * mu.plus^(-1 + i2) * (-i2 + bd * mu.plus) * (-1 + sigma) * 
                        gamma(1 + bd))/(gamma(1 + bd - i2) * gamma(1 + i2)), 
                      -((1 - sigma) * ((1 - mu.plus)^(bd - 1) * bd)))
    df2.plus <- as.vector(df2.plus)*(1-mu.plus)*mu.plus
    d2f2delta22 <- (df2.plus - df2) / fd.prec
    
    sigma.plus <- sigma+fd.prec
    d2f2sigma2 <- 0
    
    df2.dsigma.plus <- ifelse(i2>0, ((1 - mu)^(-1 + bd - i2) * mu^(-1 + i2) * (-i2 + bd * mu) * (-1 + sigma.plus) * 
                              gamma(1 + bd))/(gamma(1 + bd - i2) * gamma(1 + i2)), 
                              -((1 - sigma.plus) * ((1 - mu)^(bd - 1) * bd)))
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)*(1-mu)*mu
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec
  
    
    if(univariate==FALSE) {
    
      cdf.f.plus <- rep(0, ly)
      nmu <- rep(mu.plus, length = ly) 
      nsigma <- rep(sigma, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f.plus[i] <- sum(df2.mu.int.ZIBI(allval, mm, ms, bd))
      }  
    
      dF2.mu.plus <- cdf.f.plus
      d2F2ddelta22 <- (dF2.mu.plus - dF2) / fd.prec
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    
      cdf.f.sigma.plus <- rep(0, ly)
      cdf.f.dsigma.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma.plus, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f.sigma.plus[i] <- sum(df2.sigma.int.ZIBI(allval, mm, ms, bd))
        cdf.f.dsigma.plus[i] <- sum(df2.mu.int.ZIBI(allval, mm, ms, bd))
      }  
    
      dF2.sigma.plus <- cdf.f.sigma.plus
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
      dF2.dsigma.plus <- cdf.f.dsigma.plus
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma     
    
    }
    
    

  } else if (VC$margins[2]=="ZINBI") {
    
    i2 <- dat[,2]
    
    if(univariate==TRUE) {
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
      if (is.null(VC$X3)) nu.star <- params[(VC$X2.d2 + 2)]
      if (!is.null(VC$X3)) nu.star <- etanu <- X4 %*% params[(VC$X2.d2 + VC$X3.d2 + 1):(VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]        
      if(is.null(VC$X3))  nu.star <- params[(VC$X1.d2+VC$X2.d2+2)]
      if(!is.null(VC$X3)) nu.star <- etanu <- X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)] 
    }
      
    nu <- plogis(nu.star)
    nu <- ifelse(nu<precision, precision, nu)
    nu <- ifelse(nu>(1-precision), 1-precision, nu)
    sigma <- exp(sigma.star)
    sigma <- ifelse(sigma>precision, sigma, precision)
    sigma <- ifelse(sigma<10^7, sigma, 10^7)
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^7, mu, 10^7)
    f2 <- dZINBI(i2, mu=mu, sigma=sigma, nu=nu)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- pZINBI(i2, mu=mu, sigma=sigma, nu=nu)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
    
    df2 <- ifelse(i2>0, ((-i2 + mu) * (-1 + nu) * (1/(1 + mu * sigma))^(1 + 1/sigma) * ((mu * sigma)/(1 + mu * sigma))^i2 * 
                      gamma(i2 + 1/sigma))/(mu * gamma(1 + i2) * gamma(1/sigma)), 
                      (-1 + nu) * ((1/(1 + mu * sigma))^((1/sigma) - 1) * ((1/sigma) * 
                      (sigma/(1 + mu * sigma)^2))) )
    df2 <- as.vector(df2) * mu
    
    
    df2.sigma <- ifelse(i2>0, ((-1 + nu) * (1/(1 + mu * sigma))^(1 + 1/sigma) * ((mu * sigma)/(1 + mu * sigma))^
                    i2 * gamma(i2 + 1/sigma) * (-i2 * sigma + mu * sigma + log(1/(1 + mu * sigma)) + 
                    mu * sigma * log(1/(1 + mu * sigma)) + (1 + mu * sigma) * digamma(
                    i2 + 1/sigma) - (1 + mu * sigma) * digamma(1/
                    sigma)))/(sigma^2 * gamma(1 + i2) * gamma(1/sigma)), 
                    ((-1 + nu) * (1/(1 + mu * sigma))^(1 + 1/sigma) * (mu * sigma + (1 + mu * sigma) * log(1/(1 + mu * sigma))))/sigma^2)
    
    
    df2.nu <- ifelse(i2>0, -((1/(1 + sigma * mu))^(1/sigma) * (mu * sigma/(1 + mu * sigma))^i2 * 
                      gamma(1/sigma + i2)/(gamma(1/sigma) * gamma(1 + i2))), 
                      1 - (1/(1 + mu * sigma))^(1/sigma))
    
    
    if(univariate==FALSE) {
    
      df2.mu.int.ZINBI <- function(y, mu, sigma, nu) {
        
        df2 <- ifelse(y>0, -(((-1 + nu) * (1/(1 + sigma * mu))^(1/sigma) * ((mu * sigma/(1 + 
                      mu * sigma))^(y - 1) * (y * (sigma/(1 + mu * sigma) - mu * 
                      sigma * sigma/(1 + mu * sigma)^2))) - (-1 + nu) * ((1/(1 + 
                      sigma * mu))^((1/sigma) - 1) * ((1/sigma) * (sigma/(1 + sigma * 
                      mu)^2))) * (mu * sigma/(1 + mu * sigma))^y) * gamma(1/sigma + 
                      y)/(gamma(1/sigma) * gamma(1 + y))), 
                      (-1 + nu) * ((1/(1 + mu * sigma))^((1/sigma) - 1) * ((1/sigma) * 
                      (sigma/(1 + mu * sigma)^2))) )
        df2 <- as.vector(df2) * mu
      
      }
    
      df2.sigma.int.ZINBI<-function(y, mu, sigma, nu) {
     
        df2 <- ifelse(y>0, ((-1 + nu) * (1/(1 + mu * sigma))^(1 + 1/sigma) * ((mu * sigma)/(1 + mu * sigma))^
                    y * gamma(y + 1/sigma) * (-y * sigma + mu * sigma + log(1/(1 + mu * sigma)) + 
                    mu * sigma * log(1/(1 + mu * sigma)) + (1 + mu * sigma) * digamma(
                    y + 1/sigma) - (1 + mu * sigma) * digamma(1/
                    sigma)))/(sigma^2 * gamma(1 + y) * gamma(1/sigma)), 
                    ((-1 + nu) * (1/(1 + mu * sigma))^(1 + 1/sigma) * (mu * sigma + (1 + mu * sigma) * 
                    log(1/(1 + mu * sigma))))/sigma^2)
        df2
      
      }
    
      df2.nu.int.ZINBI <- function(y, mu, sigma, nu) {
    
        df2  <- ifelse(y>0, -((1/(1 + sigma * mu))^(1/sigma) * (mu * sigma/(1 + mu * sigma))^y * 
                      gamma(1/sigma + y)/(gamma(1/sigma) * gamma(1 + y))), 
                      1 - (1/(1 + mu * sigma))^(1/sigma)) 
        df2
      }
    
    
      ly <- length(i2)
      cdf.f <- rep(0, ly)
      cdf.f.sigma <- rep(0, ly)
      cdf.f.nu <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
      nnu <- rep(nu, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        cdf.f[i] <- sum(df2.mu.int.ZINBI(allval, mm, ms, mn))
        cdf.f.sigma[i] <- sum(df2.sigma.int.ZINBI(allval, mm, ms, mn))
        cdf.f.nu[i] <- sum(df2.nu.int.ZINBI(allval, mm, ms, mn))
      }  
    
      dF2<-cdf.f
      dF22<-dF2-df2
    
      dF2.sigma <- cdf.f.sigma
      dF22.sigma <- dF2.sigma-df2.sigma
    
      dF2.nu <- cdf.f.nu
      dF22.nu <- dF2.nu-df2.nu
    
    }
    
    
    # Hessian derivative components
    
    #fd.prec2 <- 10^-3
    fd.prec2 <- sqrt(eps)
    
    #delta2delta2
    mu.plus <- ifelse(exp(eta2+fd.prec)>0.0001, exp(eta2+fd.prec), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    df2.plus <- ifelse(i2>0, ((-i2 + mu.plus) * (-1 + nu) * (1/(1 + mu.plus * sigma))^(1 + 1/sigma) * 
                      ((mu.plus * sigma)/(1 + mu.plus * sigma))^i2 * gamma(i2 + 1/sigma))/(mu.plus * gamma(1 + i2) * gamma(1/sigma)), 
                      (-1 + nu) * ((1/(1 + mu.plus * sigma))^((1/sigma) - 1) * ((1/sigma) * 
                      (sigma/(1 + mu.plus * sigma)^2))) )
    df2.plus <- as.vector(df2.plus) * mu.plus
    d2f2delta22 <- (df2.plus - df2) / fd.prec
    
    #sigma^2
    sigma.plus <- sigma + fd.prec
    df2.sigma.plus <- ifelse(i2>0, ((-1 + nu) * (1/(1 + mu * sigma.plus))^(1 + 1/sigma.plus) * ((mu * sigma.plus)/(1 + mu * 
                    sigma.plus))^i2 * gamma(i2 + 1/sigma.plus) * (-i2 * sigma.plus + mu * sigma.plus + log(1/(1 + mu * sigma.plus)) + 
                    mu * sigma.plus * log(1/(1 + mu * sigma.plus)) + (1 + mu * sigma.plus) * digamma(i2 + 1/sigma.plus) - 
                    (1 + mu * sigma.plus) * digamma(1/sigma.plus)))/(sigma.plus^2 * gamma(1 + i2) * gamma(1/sigma.plus)), 
                    ((-1 + nu) * (1/(1 + mu * sigma.plus))^(1 + 1/sigma.plus) * (mu * sigma.plus + (1 + mu * sigma.plus) * 
                    log(1/(1 + mu * sigma.plus))))/sigma.plus^2)
    d2f2sigma2 <- (df2.sigma.plus - df2.sigma) / fd.prec
    
    #delta2sigma
    
    df2.dsigma.plus <- ifelse(i2>0, ((-i2 + mu) * (-1 + nu) * (1/(1 + mu * sigma.plus))^(1 + 1/sigma.plus) * 
                      ((mu * sigma.plus)/(1 + mu * sigma.plus))^i2 * gamma(i2 + 1/sigma.plus))/(mu * gamma(1 + i2) * gamma(1/sigma.plus)), 
                      (-1 + nu) * ((1/(1 + mu * sigma.plus))^((1/sigma.plus) - 1) * ((1/sigma.plus) * 
                      (sigma.plus/(1 + mu * sigma.plus)^2))) )
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)*mu
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec
    
    
    #nu^2
    nu.plus <- nu + fd.prec
    nu.plus <- ifelse(nu.plus>=(1-precision), (1-precision), nu.plus)
    d2f2nu2 <- 0
    
    
    #delta2nu
    df2.dnu.plus <- ifelse(i2>0, ((-i2 + mu) * (-1 + nu.plus) * (1/(1 + mu * sigma))^(1 + 1/sigma) * 
                      ((mu * sigma)/(1 + mu * sigma))^i2 * gamma(i2 + 1/sigma))/(mu * gamma(1 + i2) * gamma(1/sigma)), 
                      (-1 + nu.plus) * ((1/(1 + mu * sigma))^((1/sigma) - 1) * ((1/sigma) * 
                      (sigma/(1 + mu * sigma)^2))) )
    df2.dnu.plus <- as.vector(df2.dnu.plus)*mu
    d2f2delta2nu <- (df2.dnu.plus - df2) / fd.prec
    
    
    #nusigma
    df2.nu.sigma.plus <- ifelse(i2>0, -((1/(1 + sigma.plus * mu))^(1/sigma.plus) * (mu * sigma.plus/(1 + mu * sigma.plus))^i2 * 
                              gamma(1/sigma.plus + i2)/(gamma(1/sigma.plus) * gamma(1 + i2))), 
                              1 - (1/(1 + mu * sigma.plus))^(1/sigma.plus))
    d2f2sigmanu <- (df2.nu.sigma.plus - df2.nu) / fd.prec
    
    
    if(univariate==FALSE) {
    
      cdf.f.plus <- rep(0, ly)
      nmu <- rep(mu.plus, length = ly) 
      nsigma <- rep(sigma, length = ly)
      nnu <- rep(nu, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        pdfall <- df2.mu.int.ZINBI(allval, mm, ms, mn)
        cdf.f.plus[i] <- sum(pdfall)
      }  
    
      dF2.mu.plus <- cdf.f.plus
      d2F2ddelta22 <- (dF2.mu.plus - dF2) / fd.prec
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    
      cdf.f.sigma.plus <- rep(0, ly)
      cdf.f.dsigma.plus <- rep(0, ly)
      cdf.f.dnu.dsigma.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma.plus, length = ly)
      nnu <- rep(nu, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        cdf.f.sigma.plus[i] <- sum(df2.sigma.int.ZINBI(allval, mm, ms, mn))
        cdf.f.dsigma.plus[i] <- sum(df2.mu.int.ZINBI(allval, mm, ms, mn))
        cdf.f.dnu.dsigma.plus[i] <- sum(df2.nu.int.ZINBI(allval, mm, ms, mn))
      }  
    
      dF2.sigma.plus <- cdf.f.sigma.plus
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
      dF2.dsigma.plus <- cdf.f.dsigma.plus
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma   
    
      dF2.dnu.dsigma.plus <- cdf.f.dnu.dsigma.plus
      d2F2dsigmadnu <- (dF2.dnu.dsigma.plus - dF2.nu) / fd.prec
      d2F22dsigmadnu <- d2F2dsigmadnu - d2f2sigmanu
      
      
      
      cdf.f.nu.plus <- rep(0, ly)
      cdf.f.dnu.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
      nnu <- rep(nu.plus, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        mn<- nnu[i]
        allval <- seq(0, y.y)
        cdf.f.nu.plus[i] <- sum(df2.nu.int.ZINBI(allval, mm, ms, mn))
        cdf.f.dnu.plus[i] <- sum(df2.mu.int.ZINBI(allval, mm, ms, mn))
      } 
    
      dF2.nu.plus <- cdf.f.nu.plus
      d2F2dnu2 <- (dF2.nu.plus - dF2.nu) / fd.prec
      d2F22dnu2 <- d2F2dnu2 - d2f2nu2
    
      dF2.dnu.plus <- cdf.f.dnu.plus
      d2F2ddelta2dnu <- (dF2.dnu.plus - dF2) / fd.prec
      d2F22ddelta2dnu <- d2F2ddelta2dnu - d2f2delta2nu
    
    }
    
    

  } else if (VC$margins[2]=="ZIP") {
    
    i2 <- dat[,2]
    
    if(univariate==TRUE) {
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
    }
      
    sigma <- plogis(sigma.star)
    sigma <- ifelse(sigma<(1-precision), sigma, 1-sigma)
    sigma <- ifelse(sigma>precision, sigma, precision)
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^7, mu, 10^7)
    f2 <- dZIP(i2, sigma = sigma,  mu = mu)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- pZIP(i2, sigma = sigma, mu = mu)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
    
    df2 <- ifelse(i2>0, (1 - sigma) * (mu^((i2) - 1) * (i2))/gamma(i2 + 1) * exp(-mu) - 
                    (1 - sigma) * mu^(i2)/gamma(i2 + 1) * exp(-mu), -((1 - sigma) * exp(-mu)))
    df2 <- as.vector(df2)*mu
    
    df2.sigma <- ifelse(i2>0, -(mu^(i2)/gamma(i2 + 1) * exp(-mu)), 1 - exp(-mu)) 
    
    
    if(univariate==FALSE) {
    
      df2.mu.int.ZIP<-function(y, mu, sigma) {
        
        df2 <- ifelse(y>0, (1 - sigma) * (mu^((y) - 1) * (y))/gamma(y + 1) * exp(-mu) - 
                      (1 - sigma) * mu^(y)/gamma(y + 1) * exp(-mu), -((1 - sigma) * exp(-mu)))
        df2 <- as.vector(df2)*mu
        df2
      }
    
      df2.sigma.int.ZIP<-function(y, mu, sigma) {
        df2 <- ifelse(y>0, -(mu^(y)/gamma(y + 1) * exp(-mu)), 1 - exp(-mu)) 
        df2
      }
    
    
      ly <- length(i2)
      cdf.f <- rep(0, ly)
      cdf.f.sigma <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f[i] <- sum(df2.mu.int.ZIP(allval, mm, ms))
        cdf.f.sigma[i] <- sum(df2.sigma.int.ZIP(allval, mm, ms))
      }  
    
      dF2<-cdf.f
      dF22<-dF2-df2
    
      dF2.sigma <- cdf.f.sigma
      dF22.sigma <- dF2.sigma-df2.sigma
    
    }
    
    # Hessian derivative components
    
    
    mu.plus <- ifelse(exp(eta2+fd.prec)>0.0001, exp(eta2+fd.prec), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    
    df2.plus <- ifelse(i2>0, (1 - sigma) * (mu.plus^((i2) - 1) * (i2))/gamma(i2 + 1) * exp(-mu.plus) - 
                (1 - sigma) * mu.plus^(i2)/gamma(i2 + 1) * exp(-mu.plus), -((1 - sigma) * exp(-mu.plus)))
    df2.plus <- as.vector(df2.plus)*mu.plus
    d2f2delta22 <- (df2.plus - df2) / fd.prec
    
    
    sigma.plus <- sigma + fd.prec
    d2f2sigma2 <- 0
    
    df2.dsigma.plus <- ifelse(i2>0, (1 - sigma.plus) * (mu^((i2) - 1) * (i2))/gamma(i2 + 1) * exp(-mu) - 
                              (1 - sigma.plus) * mu^(i2)/gamma(i2 + 1) * exp(-mu), -((1 - sigma.plus) * exp(-mu)))
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)*mu
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec
    
    if(univariate==FALSE) {
    
      cdf.f.plus <- rep(0, ly)
      nmu <- rep(mu.plus, length = ly) 
      nsigma <- rep(sigma, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f.plus[i] <- sum(df2.mu.int.ZIP(allval, mm, ms))
      }  
    
      dF2.mu.plus <- cdf.f.plus
      d2F2ddelta22 <- (dF2.mu.plus - dF2) / fd.prec
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    
      cdf.f.sigma.plus <- rep(0, ly)
      cdf.f.dsigma.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma.plus, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f.sigma.plus[i] <- sum(df2.sigma.int.ZIP(allval, mm, ms))
        cdf.f.dsigma.plus[i] <- sum(df2.mu.int.ZIP(allval, mm, ms))
      }  
    
      dF2.sigma.plus <- cdf.f.sigma.plus
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
      dF2.dsigma.plus <- cdf.f.dsigma.plus
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma   
      
    
    }


  } else if (VC$margins[2]=="ZIP2") {
    
    i2 <- dat[,2]
    
    if(univariate==TRUE) {
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
    }
      
    sigma <- plogis(sigma.star)
    sigma <- ifelse(sigma<(1-precision), sigma, 1-sigma)
    sigma <- ifelse(sigma>precision, sigma, precision)
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^7, mu, 10^7)
    f2 <- dZIP2(i2, sigma = sigma,  mu = mu)
    f2 <- ifelse(f2>precision, f2, precision)
    
    if(univariate==FALSE) {
      F2 <- pZIP2(i2, sigma = sigma, mu = mu)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F2 <- ifelse(F2>precision, F2, precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision)
    }
    
    df2 <- ifelse(i2>0, -(exp(mu/(-1 + sigma))*mu^(-1 + i2)*(mu + i2*(-1 + sigma)))/((1 - sigma)^i2*gamma(1 + i2)), 
                        -((1 - sigma) * (exp(-(mu/(1 - sigma))) * (1/(1 - sigma)))) )
    df2 <- as.vector(df2)*mu
    
    df2.sigma <- ifelse(i2>0, -(exp(mu/(-1 + sigma))*mu^i2*(1 - sigma)^(-1 - i2)*(1 + mu + i2*(-1 + sigma) - sigma))/gamma(1 + i2),
                      1 - ((1 - sigma) * (exp(-(mu/(1 - sigma))) * (mu/(1 - sigma)^2)) + 
                      exp(-(mu/(1 - sigma))))) 
    
    
    if(univariate==FALSE) {
    
      df2.mu.int.ZIP2<-function(y, mu, sigma) {
        
        df2 <- ifelse(y>0, -(exp(mu/(-1 + sigma))*mu^(-1 + y)*(mu + y*(-1 + sigma)))/((1 - sigma)^y*gamma(1 + y)), 
                      -((1 - sigma) * (exp(-(mu/(1 - sigma))) * (1/(1 - sigma)))) )
        df2 <- as.vector(df2)*mu
        df2
      }
    
      df2.sigma.int.ZIP2<-function(y, mu, sigma) {
        df2 <- ifelse(y>0, -(exp(mu/(-1 + sigma))*mu^y*(1 - sigma)^(-1 - y)*(1 + mu + y*(-1 + sigma) - sigma))/gamma(1 + y),  
                      1 - ((1 - sigma) * (exp(-(mu/(1 - sigma))) * (mu/(1 - sigma)^2)) + 
                      exp(-(mu/(1 - sigma))))) 
        df2
      }
    
    
      ly <- length(i2)
      cdf.f <- rep(0, ly)
      cdf.f.sigma <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f[i] <- sum(df2.mu.int.ZIP2(allval, mm, ms))
        cdf.f.sigma[i] <- sum(df2.sigma.int.ZIP2(allval, mm, ms))
      }  
    
      dF2<-cdf.f
      dF22<-dF2-df2
    
      dF2.sigma <- cdf.f.sigma
      dF22.sigma <- dF2.sigma-df2.sigma
    
    }
    
    # Hessian derivative components
    
    
    mu.plus <- ifelse(exp(eta2+fd.prec)>0.0001, exp(eta2+fd.prec), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    
    df2.plus <- ifelse(i2>0, -(exp(mu.plus/(-1 + sigma))*mu.plus^(-1 + i2)*(mu.plus + i2*(-1 + sigma)))/((1 - sigma)^i2*gamma(1 + i2)),
                      -((1 - sigma) * (exp(-(mu.plus/(1 - sigma))) * (1/(1 - sigma)))) )
    df2.plus <- as.vector(df2.plus)*mu.plus
    d2f2delta22 <- (df2.plus - df2) / fd.prec
    
    
    sigma.plus <- sigma + fd.prec
    df2.sigma.plus <- ifelse(i2>0, -(exp(mu/(-1 + sigma.plus))*mu^i2*(1 - sigma.plus)^(-1 - i2)*(1 + mu + i2*(-1 + sigma.plus) - sigma.plus))/gamma(1 + i2),
                      1 - ((1 - sigma.plus) * (exp(-(mu/(1 - sigma.plus))) * (mu/(1 - sigma.plus)^2)) + exp(-(mu/(1 - sigma.plus)))))
    d2f2sigma2 <- (df2.sigma.plus - df2.sigma) / fd.prec
    
    df2.dsigma.plus <- ifelse(i2>0, -(exp(mu/(-1 + sigma.plus))*mu^(-1 + i2)*(mu + i2*(-1 + sigma.plus)))/((1 - sigma.plus)^i2*gamma(1 + i2)),
                              -((1 - sigma.plus) * (exp(-(mu/(1 - sigma.plus))) * (1/(1 - sigma.plus)))) )
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)*mu
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec
    
    if(univariate==FALSE) {
    
      cdf.f.plus <- rep(0, ly)
      nmu <- rep(mu.plus, length = ly) 
      nsigma <- rep(sigma, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f.plus[i] <- sum(df2.mu.int.ZIP2(allval, mm, ms))
      }  
    
      dF2.mu.plus <- cdf.f.plus
      d2F2ddelta22 <- (dF2.mu.plus - dF2) / fd.prec
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    
      cdf.f.sigma.plus <- rep(0, ly)
      cdf.f.dsigma.plus <- rep(0, ly)
      nmu <- rep(mu, length = ly) 
      nsigma <- rep(sigma.plus, length = ly)
    
      for (i in 1:ly) {
        y.y <- i2[i]
        mm <- nmu[i]
        ms <- nsigma[i]
        allval <- seq(0, y.y)
        cdf.f.sigma.plus[i] <- sum(df2.sigma.int.ZIP2(allval, mm, ms))
        cdf.f.dsigma.plus[i] <- sum(df2.mu.int.ZIP2(allval, mm, ms))
      }  
    
      dF2.sigma.plus <- cdf.f.sigma.plus
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
      dF2.dsigma.plus <- cdf.f.dsigma.plus
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma   
        
       
    
    }

  } else if (VC$margins[2]=="ZIPIG") {
    i2 <- dat[,2]
    
    if(univariate==TRUE){
      if (is.null(VC$X3)) sigma.star <- params[(VC$X2.d2 + 1)]
      if (!is.null(VC$X3)) sigma.star <- etasqv <- X3 %*% params[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X3.d2)]
      if (is.null(VC$X3)) nu.star <- params[(VC$X2.d2 + 2)]
      if (!is.null(VC$X3)) nu.star <- etanu <- X4 %*% params[(VC$X2.d2 + VC$X3.d2 + 1):(VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
    } else {
      if(is.null(VC$X3))  sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
      if(!is.null(VC$X3)) sigma.star <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]        
      if(is.null(VC$X3))  nu.star <- params[(VC$X1.d2+VC$X2.d2+2)]
      if(!is.null(VC$X3)) nu.star <- etanu <- X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)]     
    }
      
    sigma <- exp(sigma.star)
    sigma <- ifelse(sigma>precision, sigma, precision)
    sigma <- ifelse(sigma<10^7, sigma, 10^7)
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^7, mu, 10^7) 
    nu <- plogis(nu.star)
    nu <- ifelse(nu<precision, precision, nu)
    nu <- ifelse(nu>(1-precision), 1-precision, nu)   
    f2 <- dZIPIG(i2, mu=mu, sigma=sigma, nu=nu)
    f2 <- ifelse(f2>precision, f2, precision)
  
    if(univariate==FALSE) {
      F2 <- pZIPIG(i2, mu=mu, sigma=sigma, nu=nu)
      F2 <- ifelse(F2>precision, F2, precision)
      F2 <- ifelse(F2<(1-precision), F2, 1-precision)
      F22 <- F2 - f2
      F22 <- ifelse(F22<(1-precision), F22, 1-precision)
      F22 <- ifelse(F22>precision, F22, precision) 
    }
    
    f2.fd.mu <- dZIPIG(i2, mu = (mu +fd.prec), sigma = sigma, nu = nu)
    f2.fd.mu <- ifelse(f2.fd.mu>precision, f2.fd.mu, precision)
    df2 <- (f2.fd.mu - f2) / (fd.prec)
    df2 <- as.vector(df2)*as.vector(mu)
    
    f2.fd.sigma <- dZIPIG(i2, mu = mu, sigma = (sigma+fd.prec), nu = nu)
    f2.fd.sigma <- ifelse(f2.fd.sigma>precision, f2.fd.sigma, precision)
    df2.sigma <- (f2.fd.sigma - f2) / (fd.prec)
    
    f2.fd.nu <- dZIPIG(i2, mu = mu, sigma = sigma, nu =(nu+fd.prec))
    f2.fd.nu <- ifelse(f2.fd.nu>precision, f2.fd.nu, precision)
    df2.nu <- (f2.fd.nu - f2) / (fd.prec)
    
    if(univariate==FALSE) {
      F2.fd.mu <- pZIPIG(i2, mu = (mu+fd.prec), sigma = sigma, nu = nu)
      F2.fd.mu <- ifelse(F2.fd.mu<(1-precision), F2.fd.mu, 1-precision)
      dF2 <- (F2.fd.mu - F2) / (fd.prec)
      dF2 <- as.vector(dF2)*as.vector(mu) 
      dF22 <- dF2 - df2
    
      F2.fd.sigma <- pZIPIG(i2, mu = mu, sigma = (sigma+fd.prec), nu = nu)
      F2.fd.sigma <- ifelse(F2.fd.sigma<(1-precision), F2.fd.sigma, 1-precision)
      dF2.sigma <- (F2.fd.sigma - F2) / (fd.prec)
      dF22.sigma <- dF2.sigma - df2.sigma
    
      F2.fd.nu <- pZIPIG(i2, mu = mu, sigma = sigma, nu = (nu+fd.prec))
      F2.fd.nu <- ifelse(F2.fd.nu<(1-precision), F2.fd.nu, 1-precision)
      dF2.nu <- (F2.fd.nu - F2) / (fd.prec)
      dF22.nu <- dF2.nu - df2.nu
    }
    
    
    # Hessian derivative components
    
    #second derivative fd precision
    #fd.prec2 <- 10^-3
    fd.prec2 <- sqrt(eps)
    
    #pmfs
    
    #delta2^2
    mu.plus <- ifelse(exp(eta2+fd.prec2)>0.0001, exp(eta2+fd.prec2), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    f2.plus <- dZIPIG(i2, mu = mu.plus, sigma = sigma, nu = nu)
    f2.plus <- ifelse(f2.plus>precision, f2.plus, precision)
    f2.fd.mu.plus <- dZIPIG(i2, mu = (mu.plus+fd.prec), sigma = sigma, nu = nu)
    f2.fd.mu.plus <- ifelse(f2.fd.mu.plus>precision, f2.fd.mu.plus, precision)
    df2.plus <- (f2.fd.mu.plus - f2.plus) / (fd.prec)
    df2.plus <- as.vector(df2.plus)*as.vector(mu.plus)
    df2.plus <- as.vector(df2.plus)
    d2f2delta22 <- (df2.plus - df2) / fd.prec2
    
    #sigma^2
    sigma.plus <- sigma+fd.prec2
    f2.plus.sigma <- dZIPIG(i2, mu = mu, sigma = sigma.plus, nu = nu)
    f2.plus.sigma <- ifelse(f2.plus.sigma>precision, f2.plus.sigma, precision)
    f2.fd.sigma.plus <- dZIPIG(i2, mu = mu, sigma = (sigma.plus+fd.prec), nu = nu)
    f2.fd.sigma.plus <- ifelse(f2.fd.sigma.plus>precision, f2.fd.sigma.plus, precision)
    df2.sigma.plus <- (f2.fd.sigma.plus - f2.plus.sigma) / (fd.prec)
    d2f2sigma2 <- (df2.sigma.plus - df2.sigma) / fd.prec2
    
    #delta2sigma
    f2.dsigma.plus <- dZIPIG(i2, mu = mu, sigma = sigma.plus, nu = nu)
    f2.dsigma.plus <- ifelse(f2.dsigma.plus>precision, f2.dsigma.plus, precision)
    f2.fd.mu.sigma.plus <- dZIPIG(i2, mu = (mu+fd.prec), sigma = sigma.plus, nu = nu)
    f2.fd.mu.sigma.plus <- ifelse(f2.fd.mu.sigma.plus>precision, f2.fd.mu.sigma.plus, precision)
    df2.dsigma.plus <- (f2.fd.mu.sigma.plus - f2.dsigma.plus) / (fd.prec)
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)*as.vector(mu)
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec2
    
    #nu^2
    nu.plus <- nu + fd.prec2
    nu.plus <- ifelse(nu.plus>=(1-precision), (1-precision), nu.plus)
    f2.plus.nu <- dZIPIG(i2, mu = mu, sigma = sigma, nu = nu.plus)
    f2.plus.nu <- ifelse(f2.plus.nu>precision, f2.plus.nu, precision)
    nu.plus2 <- nu.plus + fd.prec
    nu.plus2 <- ifelse(nu.plus2>(1-precision), 1-precision, nu.plus2)
    f2.fd.nu.plus <- dZIPIG(i2, mu = mu, sigma = sigma, nu = nu.plus2)
    f2.fd.nu.plus <- ifelse(f2.fd.nu.plus>precision, f2.fd.nu.plus, precision)
    df2.nu.plus <- (f2.fd.nu.plus - f2.plus.nu) / (fd.prec)
    d2f2nu2 <- (df2.nu.plus - df2.nu) / fd.prec2
    
    #delta2nu
    f2.dnu.plus <- dZIPIG(i2, mu = mu, sigma = sigma, nu = nu.plus)
    f2.dnu.plus <- ifelse(f2.dnu.plus>precision, f2.dnu.plus, precision)
    f2.fd.mu.nu.plus <- dZIPIG(i2, mu = (mu+fd.prec), sigma = sigma, nu = nu.plus)
    f2.fd.mu.nu.plus <- ifelse(f2.fd.mu.nu.plus>precision, f2.fd.mu.nu.plus, precision)
    df2.dnu.plus <- (f2.fd.mu.nu.plus - f2.dnu.plus) / (fd.prec)
    df2.dnu.plus <- as.vector(df2.dnu.plus)*as.vector(mu)
    d2f2delta2nu <-(df2.dnu.plus - df2) / fd.prec2
    
    #nusigma
    f2.fd.sigma.nu.plus <- dZIPIG(i2, mu = mu, sigma = (sigma+fd.prec), nu = nu.plus)
    f2.fd.sigma.nu.plus <- ifelse(f2.fd.sigma.nu.plus>precision, f2.fd.sigma.nu.plus, precision)
    df2.dnu.sigma.plus <- (f2.fd.sigma.nu.plus - f2.plus.nu) / (fd.prec)
    df2.dnu.sigma.plus <- as.vector(df2.dnu.sigma.plus)
    d2f2sigmanu <- (df2.dnu.sigma.plus - df2.sigma) / fd.prec2
    
    
    if(univariate==FALSE) {
    
      #cdfs
    
      #delta2^2
      F2.plus <- pZIPIG(i2, mu = mu.plus, sigma = sigma, nu = nu)
      F2.fd.mu.plus <- pZIPIG(i2, mu = (mu.plus+fd.prec), sigma = sigma, nu = nu)
      F2.fd.mu.plus <- ifelse(F2.fd.mu.plus<(1-precision), F2.fd.mu.plus, 1-precision)
      dF2.plus <- (F2.fd.mu.plus - F2.plus) / (fd.prec)
      dF2.plus <- as.vector(dF2.plus)*as.vector(mu.plus) 
      d2F2ddelta22 <-(dF2.plus - dF2) / fd.prec2
      d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
      #sigma^2
      F2.sigma.plus <- pZIPIG(i2, mu = mu, sigma = sigma.plus, nu = nu)
      F2.fd.sigma.plus <- pZIPIG(i2, mu = mu, sigma = (sigma.plus+fd.prec), nu=nu)
      F2.fd.sigma.plus <- ifelse(F2.fd.sigma.plus<(1-precision), F2.fd.sigma.plus, 1-precision)
      dF2.sigma.plus <- (F2.fd.sigma.plus - F2.sigma.plus) / (fd.prec)
      d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec2
      d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
      #delta2sigma
      F2.fd.dsigma.plus <- pZIPIG(i2, mu = (mu+fd.prec), sigma = sigma.plus, nu)
      F2.fd.dsigma.plus <- ifelse(F2.fd.dsigma.plus<(1-precision), F2.fd.dsigma.plus, 1-precision)
      dF2.dsigma.plus <- (F2.fd.dsigma.plus - F2.sigma.plus) / (fd.prec)
      dF2.dsigma.plus <- as.vector(dF2.dsigma.plus)*as.vector(mu) 
      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec2
      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma
    
      #nu^2
      F2.nu.plus<- pZIPIG(i2, mu = mu, sigma = sigma, nu = nu.plus)
      F2.fd.nu.plus <- pZIPIG(i2, mu = mu, sigma = sigma, nu = nu.plus2)
      F2.fd.nu.plus <- ifelse(F2.fd.nu.plus<(1-precision), F2.fd.nu.plus, 1-precision)
      dF2.nu.plus <- (F2.fd.nu.plus - F2.nu.plus) / (fd.prec)
      d2F2dnu2 <- (dF2.nu.plus - dF2.nu) / fd.prec2
      d2F22dnu2 <- d2F2dnu2 - d2f2nu2
    
      #delta2nu
      F2.fd.dnu.plus <- pZIPIG(i2, mu = (mu+fd.prec), sigma = sigma, nu.plus)
      F2.fd.dnu.plus <- ifelse(F2.fd.dnu.plus<(1-precision), F2.fd.dnu.plus, 1-precision)
      dF2.dnu.plus <- (F2.fd.dnu.plus - F2.nu.plus) / (fd.prec)
      dF2.dnu.plus <- as.vector(dF2.dnu.plus)*as.vector(mu) 
      d2F2ddelta2dnu <- (dF2.dnu.plus - dF2) / fd.prec2
      d2F22ddelta2dnu <- d2F2ddelta2dnu - d2f2delta2nu
    
      #sigmanu
      F2.fd.dnu.sigma.plus <- pZIPIG(i2, mu = mu, sigma = (sigma+fd.prec), nu.plus)
      F2.fd.dnu.sigma.plus <- ifelse(F2.fd.dnu.sigma.plus<(1-precision), F2.fd.dnu.sigma.plus, 1-precision)
      dF2.dnu.sigma.plus <- (F2.fd.dnu.sigma.plus - F2.nu.plus) / (fd.prec)
      dF2.dnu.sigma.plus <- as.vector(dF2.dnu.sigma.plus)
      d2F2dsigmadnu <- (dF2.dnu.sigma.plus - dF2.sigma) / fd.prec2
      d2F22dsigmadnu <- d2F2dsigmadnu - d2f2sigmanu
    
    }
    
    
    
  }

   
  
  
  
  list(i0=i0, i1=i1, i2=i2, ind=ind, F1=F1, dF1=dF1, d2F1ddelta1delta1=d2F1ddelta1delta1, mu=mu, sigma=sigma, nu=nu, 
       F2=F2, f2=f2, F22=F22, df2=df2, df2.sigma=df2.sigma, df2.nu=df2.nu, dF2=dF2, dF22=dF22, dF2.sigma=dF2.sigma, 
       dF22.sigma=dF22.sigma, dF2.nu=dF2.nu, dF22.nu=dF22.nu, d2f2delta22=d2f2delta22, d2f2sigma2=d2f2sigma2, 
       d2f2delta2sigma=d2f2delta2sigma, d2f2nu2=d2f2nu2, d2f2delta2nu=d2f2delta2nu, d2f2nusigma=d2f2nusigma, 
       d2f2sigmanu=d2f2sigmanu, d2F2ddelta22=d2F2ddelta22, d2F22ddelta22=d2F22ddelta22, d2F2dsigma2=d2F2dsigma2, 
       d2F22dsigma2=d2F22dsigma2, d2F2ddelta2dsigma=d2F2ddelta2dsigma, d2F22ddelta2dsigma=d2F22ddelta2dsigma, d2F2dnu2=d2F2dnu2, 
       d2F22dnu2=d2F22dnu2, d2F2ddelta2dnu=d2F2ddelta2dnu, d2F22ddelta2dnu=d2F22ddelta2dnu, d2F2dnudsigma=d2F2dnudsigma,
       d2F22dnudsigma=d2F22dnudsigma, d2F2dsigmadnu=d2F2dsigmadnu, d2F22dsigmadnu=d2F22dsigmadnu, etasqv=etasqv, etanu=etanu)
 
  
  
  
  
 }
  
  
  