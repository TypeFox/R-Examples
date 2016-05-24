COMP_maxn <- function(lambda,nu) {
  app_En <- lambda^(1/(nu+1e-6))
  # and using var ~ En/nu
  res <- max(2,1+app_En+6*sqrt(app_En/(nu+1e-6)))  ## real, to allow continuity correction
  opt_maxn <- spaMM.getOption("COMP_maxn")
  if (res>opt_maxn) {
    res <- opt_maxn
    warning(paste("maxn truncated to",res,"for (lambda,nu)=",lambda,nu))
  }
  res
}

COMP_Z <- function(eta,nu,lambda=exp(eta),maxn=COMP_maxn(lambda,nu)){
  return(Rcpp_COMP_Z(moment=0,nu=nu,lambda=lambda,maxn=maxn))
  # R version:
  if (nu==0) {
    return(c(logScaleFac=0,scaled=as.numeric(1/(1-lambda))))
  } 
  floorn <- floor(maxn)
  epsn <- maxn - floorn
  facs <- c(1,lambda/seq(floorn+1L)^nu)
  cumprodfacs <- cumprod(facs)
  cumprodfacs[floorn+2L] <- cumprodfacs[floorn+2L]*epsn
  scaled <- sum(cumprodfacs)
  if ( ! is.finite(scaled)) {
    logfacs <- log(facs)
    cumsumlogfacs <- cumsum(logfacs)
    refi <- which.max(cumsumlogfacs)
    logScaleFac <- cumsumlogfacs[refi]
    scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  } else logScaleFac <- 0
  return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z = exp(logScaleFac)*scaled
}

COMP_Z_n <- function(eta,nu,lambda=exp(eta),maxn=COMP_maxn(lambda,nu)){
  return(Rcpp_COMP_Z(moment=1,nu=nu,lambda=lambda,maxn=maxn))
  # R version:
  if (nu==0) {
    return(c(logScaleFac=0,scaled=as.numeric(lambda/(1-lambda)^2)))
  } 
  floorn <- floor(maxn)
  epsn <- maxn - floorn
  seqn <- seq(2L,floorn+1L)
  facs <- c(lambda,(seqn/(seqn-1L))*lambda/(seqn^nu))
  cumprodfacs <- cumprod(facs)
  cumprodfacs[floorn+1L] <- cumprodfacs[floorn+1L]*epsn
  scaled <- sum(cumprodfacs)
  if ( ! is.finite(scaled)) {
    logfacs <- log(facs)
    cumsumlogfacs <- cumsum(logfacs)
    refi <- which.max(cumsumlogfacs)
    logScaleFac <- cumsumlogfacs[refi]
    scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  } else logScaleFac <- 0
  return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z_n = exp(logScaleFac)*scaled
}

COMP_Z_n2 <- function(eta,nu,lambda=exp(eta),maxn=COMP_maxn(lambda,nu)){
  return(Rcpp_COMP_Z(moment=2,nu=nu,lambda=lambda,maxn=maxn))
  # R version:
  if (nu==0) {
    return(c(logScaleFac=0,scaled=as.numeric(lambda*(1+lambda)/(1-lambda)^3)))
  } 
  floorn <- floor(maxn)
  epsn <- maxn - floorn
  seqn <- seq(2L,floorn+1L)
  facs <- c(lambda,((seqn/(seqn-1L))^2)*lambda/(seqn^nu))
  cumprodfacs <- cumprod(facs)
  cumprodfacs[floorn+1L] <- cumprodfacs[floorn+1L]*epsn
  scaled <- sum(cumprodfacs)
  if ( ! is.finite(scaled)) {
    logfacs <- log(facs)
    cumsumlogfacs <- cumsum(logfacs)
    refi <- which.max(cumsumlogfacs)
    logScaleFac <- cumsumlogfacs[refi]
    scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  } else logScaleFac <- 0
  return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z_n2 = exp(logScaleFac)*scaled
}

COMP_Z_n3 <- function(eta,nu,lambda=exp(eta),maxn=COMP_maxn(lambda,nu)){
  return(Rcpp_COMP_Z(moment=3,nu=nu,lambda=lambda,maxn=maxn))
  # R version:
  if (nu==0) {
    return(c(logScaleFac=0,scaled=as.numeric(lambda*(1+4*lambda+lambda^2)/(1-lambda)^4)))
  } 
  floorn <- floor(maxn)
  epsn <- maxn - floorn
  seqn <- seq(2L,floorn+1L)
  facs <- c(lambda,((seqn/(seqn-1L))^3)*lambda/(seqn^nu))
  cumprodfacs <- cumprod(facs)
  cumprodfacs[floorn+1L] <- cumprodfacs[floorn+1L]*epsn
  scaled <- sum(cumprodfacs)
  if ( ! is.finite(scaled)) {
    logfacs <- log(facs)
    cumsumlogfacs <- cumsum(logfacs)
    refi <- which.max(cumsumlogfacs)
    logScaleFac <- cumsumlogfacs[refi]
    scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  } else logScaleFac <- 0
  return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z_n3 = exp(logScaleFac)*scaled
}

COMP_Z_lfacn <- function(eta,nu,lambda=exp(eta),maxn=COMP_maxn(lambda,nu)){
  floorn <- floor(maxn)
  epsn <- maxn - floorn
  seqn <- seq(2L,floorn+1L)
  facs <- c(lambda,lambda/(seqn^nu)) ## n from 1 to floorn+1
  cumprodfacs <- cumprod(facs) ## sequence of lambda^n / n!^nu 
  cumsumlogn <- cumsum(c(0,log(seqn))) ## sequence of log(n!) for n from 1 to floorn+1
  cumprodfacs <- cumprodfacs*cumsumlogn ## sequence of log(n!) lambda^n / n!^nu 
  cumprodfacs[floorn+1L] <- cumprodfacs[floorn+1L]*epsn
  scaled <- sum(cumprodfacs)
  if ( ! is.finite(scaled)) {
    logfacs <- log(facs)
    cumsumlogfacs <- cumsum(logfacs)
    refi <- which.max(cumsumlogfacs)
    logScaleFac <- cumsumlogfacs[refi]
    scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  } else logScaleFac <- 0
  return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z_n = exp(logScaleFac)*scaled
}


COMP_Z_ratio <- function(Z1,Z2,log=FALSE) {
  if (log) {
    return(Z1[["logScaleFac"]]-Z2[["logScaleFac"]]+log(Z1[["scaled"]])-log(Z2[["scaled"]]))
  } else return(exp(Z1[["logScaleFac"]]-Z2[["logScaleFac"]])*Z1[["scaled"]]/Z2[["scaled"]])
} 

COMP_dnu_objfn <- function(nu,y,eta,lambda=exp(eta),maxn=COMP_maxn(lambda,nu)) {
  summand <- function(yi,lambdai) {
    comp_z_lfacn <- COMP_Z_lfacn(nu=nu,lambda=lambdai,maxn=maxn)
    comp_z <- COMP_Z(nu=nu,lambda=lambdai,maxn=maxn)
    lfactorial(yi)-COMP_Z_ratio(comp_z_lfacn,comp_z,log=TRUE)
  }
  res <- sum(sapply(seq_len(length(y)),function(i) summand(yi=y[i],lambdai=lambda[i])))
  res
}

dCOMP <- function(x, mu, nu,
                  lambda=COMPoisson(nu=nu)$linkfun(mu,log=FALSE),
                  log = FALSE, maxn=COMP_maxn(lambda,nu)) {
  compz <- COMP_Z(lambda=lambda,nu=nu,maxn=maxn)
  logd <- x * log(lambda) - nu* lfactorial(x) - compz[["logScaleFac"]] -log(compz[["scaled"]])
  if (log) { 
    return(logd)
  } else return(exp(logd))
}

COMPoisson <- function (nu = stop("'nu' must be specified"), 
                        link = "loglambda" # eta <-> mu link, not the eta <-> lambda log link
                        ) {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) 
    linktemp <- deparse(linktemp)
  okLinks <- c("loglambda")
  if (linktemp %in% okLinks) {} else {
      stop(gettextf("link \"%s\" not available for COMPoisson family; available links are %s", 
                    linktemp, paste(sQuote(okLinks), collapse = ", ")), 
           domain = NA)
  }
  initialize <- expression({
    if (any(y < 0)) 
      stop("negative values not allowed for the 'COMPoisson' family")
    n <- rep.int(1, nobs)
    mustart <- y + 0.1
  })
  linkinv <- function(eta,lambda=exp(eta)) {
    if (! is.null(mu <- attr(lambda,"mu"))) {
      attributes(lambda) <- NULL
      return(structure(mu,lambda=lambda))
    }
    if (nu==0) {
      mus <- lambda/(1-lambda) 
    } else {
      locfn <- function(lambda) {
        num <- COMP_Z_n(lambda=lambda,nu=nu)
        denum <- COMP_Z(lambda=lambda,nu=nu)
        mu <- COMP_Z_ratio(num,denum)
        if ( ! is.finite(mu) && lambda>10^nu) { ## FR->FR heuristic
          mu <- lambda^(1/nu)-(nu-1)/(2*nu)
        }  
        return(max(mu,1e-8)) ## avoids mu=0 which fails validmu(mu) (as for poisson family)
      }
      mus <- sapply(lambda,locfn)
    }
    return(mus)
  }
  linkfun <- function(mu, ## scalar or vector
                      log=TRUE) {
    if ( is.null(lambdas <- attr(mu,"lambda"))) {
      if (nu==0) {
        lambdas <- mu/(1+mu) 
      } else lambdas <- sapply(mu,function(mu) {
        if (nu==1) {
          lambda <- mu ## pb du code general est qu'objfn n'est alors que l'erreur numérique de linkinv()
        } else if (mu==Inf) {
          warning(paste("Approximating lambda as 1-1e-8 in COMPoisson(nu=",nu,") for mu = Inf"))
          lambda <- 1 -1e-8 ## 1 -> mu.eta = NaN. mu.eta being Inf would not be OK bc Inf GLMweights -> 1/w.resid=0 -> singular Sig or d2hdv2
        } else {
          objfn <- function(lambda) {linkinv(lambda=lambda) -mu}
          app_lambda <-  max(c(0,(mu+max(0,(nu-1)/(2*nu)))))^(nu)
          lambdamin <- max(c(0,app_lambda-1))
          lambdamax <- max(c(0,(app_lambda+1)*c(1.01))) 
          # last one for low lambda,nu values
          fupper <- objfn(lambdamax)  ## normally >0 
          if (is.nan(fupper) || fupper<0) {
            warning(paste("Trying geometric approximation in COMPoisson(nu=",nu,")$linkfun(mu=",mu,")..."))
            lambda <- mu/(1+mu) 
          } else {
            flower <- objfn(lambdamin) ## normally <0 
            interval <- c(lambdamin,lambdamax)
            # linkinv(0)=0 => objfn(0)= -mu
            lambda <- uniroot(objfn,interval=interval,f.lower = flower,f.upper = fupper)$root
          }
        }
        #if(is.nan(lambda)) stop("is.nan(lambda)")
        return(lambda)
      })
    } else attributes(mu) <- NULL ## avoids 'mise en abime'
    if (log) {
      return(log(lambdas)) ## eta, ie standard linkfun value
    } else return(structure(lambdas,mu=mu))
  }
  variance <- function(mu) {
    if (nu==0) {
      return(mu*(1+mu)) 
    } else {
      lambdas <- linkfun(mu,log=FALSE)
      En <- sapply(lambdas,function(lambda) {
        COMP_Z_ratio(COMP_Z_n(lambda=lambda,nu=nu),COMP_Z(lambda=lambda,nu=nu))
      }) 
      En2 <- sapply(lambdas,function(lambda) {
        COMP_Z_ratio(COMP_Z_n2(lambda=lambda,nu=nu),COMP_Z(lambda=lambda,nu=nu)) 
      })
      resu <- pmax(En2-En^2,1e-8) ## pmax otherwise for low mu, Vmu=0, -> ... -> w.resid=0
      ## for geom V(mu) = mu(1+mu) but this cannot be a lower bound for resu.
      return(resu)
    }
  }
  validmu <- function(mu) all(is.finite(mu)) && all(mu > 0) ## from poisson()
  valideta <- function(eta) TRUE ## from poisson()
  mu.eta <- function (eta,lambda=exp(eta)) {
    if (nu==0) {
      resu <- lambda/(1-lambda)^2
    } else resu <- sapply(lambda, function(lambda) {
      compz <- COMP_Z(lambda=lambda,nu=nu)
      compzn <- COMP_Z_n(lambda=lambda,nu=nu)
      compzn2 <- COMP_Z_n2(lambda=lambda,nu=nu)
      rn2 <- COMP_Z_ratio(compzn2,compz)
      rn <- COMP_Z_ratio(compzn,compz) ## mu
      res <- rn2 - rn^2 # (compz*compzn2-compzn^2)/(compz^2)
      # dmu/deta=(dmu/dlam) (dlam/deta) = lam dmu/dlam = lam (compz*compzn2-compzn^2)/(lam compz^2) = resu
      # b/c mu = compzn/compz, d compz/ dlam = compzn/lam, d compzn/ dlam = compzn2/lam
      pmax(res, .Machine$double.eps)
    })
    resu
  }
  dev.resids <- function(y, mu, wt){
    # must accept, among others, vector y and scalar mu.
    lambdas <- linkfun(mu=mu,log=FALSE)
    n <- length(y)
    if (length(mu)==1L) lambdas <- rep(lambdas,n)
    calc_dev <- function(yi,lambdai) {
      Z2 <- COMP_Z(lambda=lambdai,nu=nu)
      if (yi==0) { # lambda = 0,  Z1 = 1
        dev <- 2*(Z2[["logScaleFac"]]+log(Z2[["scaled"]]))
      } else {
        # eval lambda for Z1(lambda...)
        if (nu==0) {
          lambda <- yi/(1+yi)  
        } else {
          objfn <- function(lambda) {linkinv(lambda=lambda) -yi}
          app_lambda <-  max(c(0,(yi+max(0,(nu-1)/(2*nu)))))^(nu)
          lambdamin <- max(c(0,app_lambda-1))
          lambdamax <- max(c(0,(app_lambda+1)*c(1.01))) 
          while (objfn(lambdamax)<0) {
            lambdamax <- lambdamax*10 ## greedy but resolves errors for low nu 
            if (lambdamax>1e100) stop("lambdamax>1e100 in COMPoisson$dev.resids()")
          }
          interval <- c(lambdamin,lambdamax)
          lambda <- uniroot(objfn,interval=interval)$root
        }
        Z1 <- COMP_Z(lambda=lambda,nu=nu)
        #
        dev <- 2*(yi*log(lambda/lambdai)-(Z1[["logScaleFac"]]-Z2[["logScaleFac"]]+log(Z1[["scaled"]]/Z2[["scaled"]])))
      }
      dev
    }
    devs <- numeric(n)
    for(i in seq(n)) { devs[i] <- calc_dev(y[i],lambdai=lambdas[i]) }
    devs <- devs*wt
    devs[devs==Inf] <- .Machine$double.xmax/n ## so that total deviance may be finite 
    return(devs) 
  }
  aic <- function(y, n, mu, wt, dev) {
    aici <- sapply(seq(length(y)),function(i){
      dCOMP(y[i], mu[i], nu=nu, log = TRUE)
    })
    -2 * sum(aici * wt)
  }
  initialize <- expression({
    if (any(y < 0)) stop("negative values not allowed for the 'Poisson' family")
    n <- rep.int(1, nobs)
    mustart <- y + 0.1
  })
  simfun <- function(object, nsim) {
    wts <- object$prior.weights
    if (any(wts != 1)) 
      warning("ignoring prior weights")
    ftd <- fitted(object)
    rpois(nsim * length(ftd), ftd)
  }
  structure(list(family = structure("COMPoisson",
                                    withArgs=quote(paste("COMPoisson(",nu,")",sep=""))), 
                 link = linktemp, linkfun = linkfun, 
                 linkinv = linkinv, variance = variance, dev.resids = dev.resids, 
                 aic = aic, mu.eta = mu.eta, initialize = initialize, 
                 validmu = validmu, valideta = valideta, 
                 simulate = function(...) {
                   message("simulate function missing from COMPoisson")
                   }), 
            class = "family")
}

