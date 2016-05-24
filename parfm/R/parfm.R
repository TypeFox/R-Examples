################################################################################
#  Parametric frailty models fitting                                           #
################################################################################
#                                                                              #
#  This function is the core of the package,                                   #
#   performing estimation and returning results.                               #
#                                                                              #
#  It is also the front-end of the package,                                    #
#   being the only that will be visible to the user.                           #
#                                                                              #
#  Its parameters are                                                          #
#   - formula  : a formula object, with the response                           #
#                on the left of a ~ operator,                                  #
#                and the terms on the right.                                   #
#                The response must be a survival object                        #
#                as returned by the Surv() function;                           #
#   - cluster  : the name of the variable in data containing cluster IDs;      #
#   - strata   : the name of the variable in data containing strata IDs;       #
#   - data     : a data.frame containing all the variables;                    #
#   - inip     : initial values for the baseline hazard & reg parameters;      #
#   - iniFpar  : initial value(s) for the heterogeneity parameter(s);          #
#   - dist     : the baseline hazard;                                          #
#   - frailty  : the frailty distribution;                                     #
#   - method   : the optimisation method (See optim());                        #
#   - maxit    : the maximum number of iterations (See optim());               #
#   - Fparscale: the scaling value for all the frailty parameter(s) in optim() #
#                Optimisation is performed on Fpar/Fparscale                   #
#   - showtime : is the execution time displayed?                              #
#   - correct  : (only for possta) the correction to use in case of many       #
#                events per cluster to get finite likelihood values.           #
#                When correct!=0 the likelihood is divided by                  #
#                10^(#clusters * correct) for computation,                     #
#                but the value of the log-likelihood in the output             #
#                is the re-adjusted value.                                     #
#                                                                              #
#                                                                              #
#   Date:                 December 21, 2011                                    #
#   Last modification on: July 6, 2015                                         #
################################################################################

parfm <- function(formula,
                  cluster=NULL,
                  strata=NULL,
                  data,
                  inip=NULL,
                  iniFpar=NULL,
                  dist="weibull",
                  frailty="none",
                  method="BFGS",
                  maxit=500,
                  Fparscale=1,
                  showtime=TRUE,
                  correct=0){
  if (missing(data)) {
    data <- eval(parse(text=paste("data.frame(", 
                                  paste(all.vars(formula), collapse=", "),
                                  ")")))
  }
  
  #----- Check the baseline hazard and the frailty distribution ---------------#
  if (!(dist %in% 
    c("exponential", "weibull",
      "gompertz", "loglogistic", "lognormal"))) {
    stop("invalid baseline hazard")
  }
  if (!(frailty %in% 
    c("none", "gamma", "ingau", "possta", "lognormal"))) {
    stop("invalid frailty distribution")
  }
  if (frailty == "none" &&  !is.null(cluster)) {
    warning(paste("With frailty='none' the cluster variable '",
                  cluster, "' is not used!", sep=""))
  }
  if (frailty == "none" &&  !is.null(iniFpar)) {
    warning("With frailty='none' the argument 'iniFpar' is not used!")
  }
  
  #----- 'Correct' is useless except for frailty="possta" ---------------------#
  if (frailty == "possta") {  #Do not exaggerate when setting 'correct' !
    if (10^correct == Inf || 10^-correct == 0) {
      stop("'correct' is too large!")
    }
    if (10^correct == 0 || 10^-correct == Inf) {
      stop("'correct' is too small!")
    }
  } else if (correct != 0) {
    warning(paste("'correct' has no effect when 'frailty = ", frailty, "'",
                  sep=""))
  }
  
  #----- Data for Mloglikelihood() --------------------------------------------#
  obsdata <- NULL
  
  #time
  if (length(formula[[2]]) == 3) {          # --> without left truncation
    obsdata$time <- eval(#parse(text=paste("data$", 
      formula[[2]][[2]], 
      envir=data #sep=""))
    )
    obsdata$event <- eval(#parse(text=paste("data$", 
      formula[[2]][[3]],
      envir=data #sep=""))
    )
  } else if (length(formula[[2]]) == 4) {   # --> with left truncation
    obsdata$trunc <- eval(#parse(text=paste("data$", 
      formula[[2]][[2]],
      envir=data #sep=""))
    )    
    obsdata$time <- eval(#parse(text=paste("data$", 
      formula[[2]][[3]] ,
      envir=data #sep=""))
    )
    obsdata$event <- eval(#parse(text=paste("data$", 
      formula[[2]][[4]], 
      envir=data #sep=""))
    )
  }
  if (!all(levels(as.factor(obsdata$event)) %in% 0:1)) {
    stop(paste("The status indicator 'event' in the Surv object",
               "in the left-hand side of the formula object",
               "must be either 0 (no event) or 1 (event)."))
  }
  
  #covariates (an intercept is automatically added)
  obsdata$x <- as.data.frame(model.matrix(formula, data=data))
  
  #cluster
  if (is.null(cluster)) {
    if (frailty != "none") {
      stop(paste("if you specify a frailty distribution,\n",
                 "then you have to specify the cluster variable as well"))
    } else {
      obsdata$cluster <- rep(1, nrow(data))
    }
    #number of clusters
    obsdata$ncl <- 1
    #number of events in each cluster
    obsdata$di <- sum(obsdata$event)
  } else {
    if (! cluster %in% names(data)) {
      stop(paste("object '", cluster, "' not found", sep=""))
    }
    obsdata$cluster <- eval(#parse(text=paste("data$", 
      as.name(cluster),
      envir=data #sep=""))
    )
    #number of clusters
    obsdata$ncl <- length(levels(as.factor(obsdata$cluster)))
    #number of events in each cluster
    obsdata$di <- aggregate(obsdata$event,
                            by=list(obsdata$cluster), 
                            FUN=sum)[,, drop=FALSE]
    cnames <- obsdata$di[,1]
    obsdata$di <- as.vector(obsdata$di[,2])
    names(obsdata$di) <- cnames
  }
  
  #strata
  if (is.null(strata)) {
    obsdata$strata <- rep(1, length(obsdata$time))
    #number of strata
    obsdata$nstr <- 1
    #number of events in each stratum
    obsdata$dq <- sum(obsdata$event)
  } else {
    if (! strata %in% names(data)) {
      stop(paste("object '", strata, "' not found", sep=""))
    }
    obsdata$strata <- eval(#parse(text=paste("data$", 
      as.name(strata), 
      envir=data #sep=""))
    )
    #number of strata
    obsdata$nstr <- length(levels(as.factor(obsdata$strata)))
    #number of events in each stratum
    obsdata$dq <- aggregate(obsdata$event, 
                            by=list(obsdata$strata), 
                            FUN=sum)[, ,drop=FALSE]
    snames <- obsdata$dq[,1]
    obsdata$dq <- as.vector(obsdata$dq[,2])
    names(obsdata$dq) <- snames
  }
  
  #cluster+strata
  if (!is.null(cluster) && !is.null(strata)) {
    #number of events in each cluster for each stratum
    obsdata$dqi <- xtabs(x~Group.1+Group.2,
                         data=aggregate(obsdata$event, 
                                        by=list(obsdata$cluster, obsdata$strata), 
                                        FUN=sum))
    dimnames(obsdata$dqi) <- list(cluster=dimnames(obsdata$dqi)[[1]], 
                                  strata =dimnames(obsdata$dqi)[[2]])
  } else if (!is.null(cluster)) {
    obsdata$dqi <- obsdata$di
  } else if (!is.null(strata)) {
    obsdata$dqi <- obsdata$dq
  } else {
    obsdata$dqi <- sum(obsdata$event)
  }
  
  
  
  #----- Dimensions -----------------------------------------------------------#
  #nFpar: number of heterogeneity parameters
  if (frailty == "none") {
    nFpar <- 0
  } else if (frailty %in% c("gamma", "ingau", "possta", "lognormal")) {
    nFpar <- 1
  }
  obsdata$nFpar <- nFpar
  
  #nBpar: number of parameters in the baseline hazard
  if (dist == "exponential") {
    nBpar <- 1
  } else if (dist %in% c("weibull", "gompertz",
                         "lognormal", "loglogistic")) {
    nBpar <- 2
  }
  #   if (!is.null(strata)) {
  #     nBpar <- nBpar #* obsdata$nstr
  #   }
  obsdata$nBpar <- nBpar
  
  #nRpar: number of regression parameters
  nRpar <- ncol(obsdata$x) - 1
  obsdata$nRpar <- nRpar  
  
  #----- Initial parameters ---------------------------------------------------#
  if (!is.null(inip)) {
    #if they are specified, then 
    #(1) check the dimension,
    #(2) check whether they lie in their parameter space, and
    #(3) reparametrise them so that they take on values on the whole real line
    
    if (length(inip) != nBpar * obsdata$nstr + nRpar) {
      stop(paste("number of initial parameters 'inip' must be", 
                 nBpar * obsdata$nstr + nRpar))
    }
    p.init <- inip
    if (dist %in% c("exponential", "weibull", "gompertz")) {
      #1st initial par: log(lambda), log(rho), log(rho), or log(gamma)
      if (any(p.init[1:obsdata$nstr] <= 0)) {
        stop(paste("with that baseline, the 1st parameter has to be > 0"))
      }
      p.init[1:obsdata$nstr] <- log(p.init[1:obsdata$nstr]) 
    }
    if (dist %in% c("weibull", "gompertz", 
                    "lognormal", "loglogistic")) {
      #2nd initial par: log(lambda), log(lambda), log(lambda), 
      #                 log(sigma), or log(kappa)
      if (any(p.init[obsdata$nstr + 1:obsdata$nstr] <= 0)) {
        stop(paste("with that baseline, the 2nd parameter has to be > 0"))
      }
      p.init[obsdata$nstr + 1:obsdata$nstr] <- 
        log(p.init[obsdata$nstr + 1:obsdata$nstr]) 
    }
  } else {    
    coxformula <- formula
    if (!is.null(strata)) {
      if (dist!="exponential") {
        coxformula <- eval(parse(text=
          paste("update(formula, .~.+strata(", strata, "))", sep="")))
      } else {
        coxformula <- eval(parse(text=
          paste("update(formula, .~.+", strata, ")", sep="")))
      }
    }
    
    #if they are not specified, then fit a parametric Cox's model
    #     library(eha)
    
    shape <- 0  #if zero or negative, the shape parameter is estimated
    d <- dist
    if (d == "exponential") {
      d <- "weibull"
      shape <- 1  #if positive, the shape parameter is fixed at that value
    }
    
    coxMod <- phreg(formula=coxformula, data=data,
                    dist=d, shape=shape, 
                    center=FALSE, control=list(maxiter=maxit))
    
    logshape <- as.numeric(coxMod$coef[
      substr(names(coxMod$coef), 5, 9) == "shape"])
    logscale <- as.numeric(coxMod$coef[
      substr(names(coxMod$coef), 5, 9) == "scale"])
    if (!is.null(strata) && dist=="exponential") {
      logscale <- logscale - 
        c(0, coxMod$coef[substr(names(coxMod$coef), 1, nchar(strata)) == strata])
    }
    
    
    if (dist == "exponential") {
      p.init <- - logscale                        #log(lambda)
    } else if (dist == "weibull") {
      p.init <- c(logshape,                       #log(rho)
                  - exp(logshape) * logscale)     #log(lambda)
    } else if (dist == "gompertz") {
      p.init <- c(- logscale,                     #log(gamma)
                  logshape)                       #log(lambda)
    } else if (dist == "lognormal") {
      p.init <- c(logscale,                       #mu
                  - logshape)                     #log(sigma)
    } else if (dist == "loglogistic") {
      p.init <- c(- exp(logshape) * logscale,     #alpha
                  logshape)                       #log(kappa)
    }
    
    if (nRpar > 0) {
      p.init <- c(p.init,
                  #                   as.numeric(coxMod$coef[1:nRpar]))   
                  as.numeric(coxMod$coef[
                    setdiff(names(coxMod$coef), 
                            c("(Intercept)", "log(scale)", "log(shape)"))
                    ]))   
    }
  }
  
  #--- frailty parameters initialisation ---#
  if (frailty == "none") {
    pars <- NULL
  } else if (frailty %in% c("gamma", "ingau", "lognormal")) {
    if (is.null(iniFpar)) {
      iniFpar <- 1
    } else if (iniFpar <= 0) {
      stop("initial heterogeneity parameter (theta) has to be > 0")
    }
    pars <- log(iniFpar)
  } else if (frailty == "possta") {
    if (is.null(iniFpar)) {
      iniFpar <- 0.5
    } else if (iniFpar <= 0 || iniFpar >= 1) {
      stop("initial heterogeneity parameter (nu) must lie in (0, 1)")
    }
    pars <- log(-log(iniFpar))
  }
  
  pars <- c(pars, p.init)
  res <- NULL
  
  #----- Minimise Mloglikelihood() --------------------------------------------#
  if ((frailty == "none") && is.null(inip)) {
    todo <- expression({res <- list(par=pars)})
    if (showtime) {
      extime <- system.time(eval(todo))[1]
    } else {
      eval(todo)
      extime <- NULL
    } 
    it <- NULL
    lL <- coxMod$loglik[2]
  } else {
    todo <- expression({
      res <- optim(par=pars, fn=Mloglikelihood, method=method, 
                   obs=obsdata, dist=dist, frailty=frailty,
                   correct=correct,
                   hessian=TRUE, 
                   control=list(maxit=maxit,
                                parscale=c(rep(Fparscale, nFpar),
                                           rep(1, nBpar  * obsdata$nstr + 
                                             nRpar))))})
    if (showtime) {
      extime <- system.time(eval(todo))[1]
    } else {
      eval(todo)
      extime <- NULL
    }
    
    if(res$convergence > 0) {
      warning("optimisation procedure did not converge,
              conv = ", bquote(.(res$convergence)), ": see optim() for details")
    }
    it <- res$counts[1]   #number of iterations
    lL <- - res$value     #maximum value of the marginal loglikelihood
    if (frailty == "possta") {
      lL <- lL + correct * log(10) * obsdata$ncl
    }
  }
  
  #----- Recover the estimates ------------------------------------------------#
  #heterogeneity parameter
  if (frailty %in% c("gamma", "ingau")) {
    theta <- exp(res$par[1:nFpar])
    sigma2 <- NULL
    nu <- NULL
  } else if (frailty == "lognormal") {
    theta <- NULL
    sigma2 <- exp(res$par[1:nFpar])
    nu <- NULL
  } else if (frailty == "possta") {
    theta <- NULL
    sigma2 <- NULL
    nu <- exp(-exp(res$par[1:nFpar]))
  } else if (frailty == "none"){
    theta <- NULL
    sigma2 <- NULL
    nu <- NULL
  }
  
  #baseline hazard parameter(s)
  if (dist == "exponential") {
    lambda <- exp(res$par[nFpar + 1:obsdata$nstr])
    ESTIMATE <- c(lambda=lambda)
  } else if (dist %in% c("weibull")) {
    rho <- exp(res$par[nFpar + 1:obsdata$nstr])
    lambda <- exp(res$par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(rho=rho, lambda=lambda)
  } else if (dist == "gompertz") {
    gamma <- exp(res$par[nFpar + 1:obsdata$nstr])
    lambda <- exp(res$par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(gamma=gamma, lambda=lambda)
  } else if (dist == "lognormal") {
    mu <- res$par[nFpar + 1:obsdata$nstr]
    sigma <- exp(res$par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(mu=mu, sigma=sigma)
  } else if (dist == "loglogistic") {
    alpha <- res$par[nFpar + 1:obsdata$nstr]
    kappa <- exp(res$par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(alpha=alpha, kappa=kappa)
  }
  
  #regression parameter(s)
  if (nRpar == 0) {
    beta <- NULL
  } else {
    beta <- res$par[-(1:(nFpar + nBpar * obsdata$nstr))]
    names(beta) <- paste("beta", names(obsdata$x), sep=".")[-1]
  }
  
  #all together
  ESTIMATE <- c(theta=theta,
                sigma2=sigma2,
                nu=nu,
                ESTIMATE,
                beta=beta)
  
  #----- Recover the standard errors ------------------------------------------#
  if ((frailty == "none") && is.null(inip)) {
    var <- coxMod$var
    if (nRpar == 0) {
      seBeta <- NULL
    } else {
      seBeta <- sqrt(diag(var)[
        setdiff(names(coxMod$coef), 
                c("(Intercept)", "log(scale)", "log(shape)"))])
      PVAL <- c(rep(NA, nFpar + nBpar * obsdata$nstr), 
                2 * pnorm(q=- abs(beta / seBeta)))
    }
    
    if (obsdata$nstr == 1) {
      if (dist == "exponential") {
        seLambda <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g=~exp(- x1), 
                      mean=logscale[x],
                      cov=var["log(scale)", "log(scale)"],
                      ses=TRUE)
        })
        STDERR <- c(seLambda=seLambda)
      } else if (dist %in% c("weibull")) {
        seRho <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g=~exp(x1), 
                      mean=logshape[x],
                      cov=var["log(shape)", "log(shape)"], 
                      ses=TRUE)
        })
        seLambda <- deltamethod(g=~exp(- exp(x2) * x1),
                                mean=c(logscale, logshape),
                                cov=var[c("log(scale)","log(shape)"),
                                        c("log(scale)","log(shape)")],
                                ses=TRUE)
        STDERR <- c(seRho=seRho, seLambda=seLambda)
      } else if (dist == "gompertz") {
        seGamma <- deltamethod(g=~exp(- x1),
                               mean=logscale,
                               cov=var["log(scale)", "log(scale)"],
                               ses=TRUE)
        seLambda <- deltamethod(g=~exp(x1),
                                mean=logshape,
                                cov=var["log(shape)", "log(shape)"],
                                ses=TRUE)
        STDERR <- c(seGamma=seGamma, seLambda=seLambda)
      } else if (dist == "lognormal") {
        seMu <- sqrt(var["log(scale)", "log(scale)"])
        seSigma <- deltamethod(g=~exp(- x1), 
                               mean=logshape, 
                               cov=var["log(shape)", "log(shape)"], 
                               ses=TRUE)
        STDERR <- c(seMu=seMu, seSigma=seSigma)
      } else if (dist == "loglogistic") {
        seAlpha <- deltamethod(g=~- exp(x2) * x1,
                               mean=c(logscale, logshape),
                               cov=var[c("log(scale)","log(shape)"),
                                       c("log(scale)","log(shape)")],
                               ses=TRUE)
        seKappa <- deltamethod(g=~exp(x1), 
                               mean=logshape, 
                               cov=var["log(shape)", "log(shape)"], 
                               ses=TRUE)
        STDERR <- c(seAlpha=seAlpha, seKappa=seKappa)
      }
    } else {
      if (dist == "exponential") {
        seLambda <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g=~exp(- x1), 
                      mean=logscale[x],
                      cov=var[paste("log(scale)", x, sep=":"), 
                              paste("log(scale)", x, sep=":")],
                      ses=TRUE)
        })
        STDERR <- c(seLambda=seLambda)
      } else if (dist %in% c("weibull")) {
        seRho <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g=~exp(x1), 
                      mean=logshape[x],
                      cov=var[paste("log(shape)", x, sep=":"), 
                              paste("log(shape)", x, sep=":")], 
                      ses=TRUE)
        })
        seLambda <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g=~exp(- exp(x2) * x1),
                      mean=c(logscale[x], logshape[x]),
                      cov=var[paste(c("log(scale)","log(shape)"), x, sep=":"), 
                              paste(c("log(scale)","log(shape)"), x, sep=":")],
                      ses=TRUE)
        })
        STDERR <- c(seRho=seRho, seLambda=seLambda)
      } else if (dist == "gompertz") {
        seGamma <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g=~exp(- x1),
                      mean=logscale[x],
                      cov=var[paste("log(scale)", x, sep=":"), 
                              paste("log(scale)", x, sep=":")],
                      ses=TRUE)
        })
        seLambda <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g=~exp(x1),
                      mean=logshape,
                      cov=var[paste("log(shape)", x, sep=":"), 
                              paste("log(shape)", x, sep=":")],
                      ses=TRUE)
        })
        STDERR <- c(seGamma=seGamma, seLambda=seLambda)
      } else if (dist == "lognormal") {
        seMu <- sqrt(diag(var[substr(rownames(var),5,9) == "scale",
                              substr(rownames(var),5,9) == "scale"]))
        seSigma <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g=~exp(- x1), 
                      mean=logshape, 
                      cov=var[paste("log(shape)", x, sep=":"), 
                              paste("log(shape)", x, sep=":")], 
                      ses=TRUE)
        })
        STDERR <- c(seMu=seMu, seSigma=seSigma)
      } else if (dist == "loglogistic") {
        seAlpha <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g=~- exp(x2) * x1,
                      mean=c(logscale, logshape),
                      cov=var[paste(c("log(scale)","log(shape)"), x, sep=":"), 
                              paste(c("log(scale)","log(shape)"), x, sep=":")],
                      ses=TRUE)
        })
        seKappa <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g=~exp(x1), 
                      mean=logshape, 
                      cov=var[paste("log(shape)", x, sep=":"), 
                              paste("log(shape)", x, sep=":")], 
                      ses=TRUE)
        })
        STDERR <- c(seAlpha=seAlpha, seKappa=seKappa)
      }
    }
    STDERR <- c(STDERR,
                se.beta=seBeta)
  } else {
    var <- try(diag(solve(res$hessian)), silent=TRUE)
    if (class(var) == "try-error") {
      warning(var[1])
      STDERR <- rep(NA, nFpar + nBpar * obsdata$nstr + nRpar)
      PVAL <- rep(NA, nFpar + nBpar * obsdata$nstr + nRpar)
    } else {
      if (any(var <= 0)) {
        warning(paste("negative variances have been replaced by NAs\n",
                      "Please, try other initial values",
                      "or another optimisation method"))
      }
      
      #heterogeneity parameter(s)
      if (frailty %in% c("gamma", "ingau")) {
        seTheta <- sapply(1:nFpar, function(x){
          ifelse(var[x] > 0, sqrt(var[x] * theta[x]^2), NA)
        })
        seSigma2 <- seNu <- NULL
      } else if (frailty == "lognormal") {
        seSigma2 <- sapply(1:nFpar, function(x){
          ifelse(var[x] > 0, sqrt(var[x] * sigma2[x]^2), NA)
        })
        seTheta <- seNu <- NULL
      } else if (frailty == "possta") {
        seNu <- sapply(1:nFpar, function(x){
          ifelse(var[x] > 0, sqrt(var[x] * (nu * log(nu))^2), NA)
        })
        seTheta <- seSigma2 <- NULL
      }
      
      #baseline hazard parameter(s)
      if (dist == "exponential") {
        seLambda <- sapply(1:obsdata$nstr, function(x){
          ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x] * lambda[x]^2), NA)
        })
        STDERR <- c(seLambda=seLambda)
      } else if (dist %in% c("weibull")) {
        seRho <- sapply(1:obsdata$nstr, function(x){
          ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x] * rho[x]^2), NA)
        })
        seLambda <- sapply(1:obsdata$nstr, function(x){
          ifelse(var[nFpar + obsdata$nstr + x] > 0, 
                 sqrt(var[nFpar + obsdata$nstr + x] * lambda[x]^2), NA)
        })
        STDERR <- c(seRho=seRho, seLambda=seLambda)
      } else if (dist == "gompertz") {
        seGamma <- sapply(1:obsdata$nstr, function(x){
          ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x] * gamma[x]^2), NA)
        })
        seLambda <- sapply(1:obsdata$nstr, function(x){
          ifelse(var[nFpar + obsdata$nstr + x] > 0,
                 sqrt(var[nFpar + obsdata$nstr + x] * lambda[x]^2), NA)
        })
        STDERR <- c(seGamma=seGamma, seLambda=seLambda)
      } else if (dist == "lognormal") {
        seMu <- sapply(1:obsdata$nstr, function(x){
          ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), NA)
        })
        seSigma <- sapply(1:obsdata$nstr, function(x){
          ifelse(var[nFpar + obsdata$nstr + x] > 0,
                 sqrt(var[nFpar + obsdata$nstr + x] * sigma[x]^2), NA)
        })
        STDERR <- c(seMu=seMu, seSigma=seSigma)
      } else if (dist == "loglogistic") {
        seAlpha <- sapply(1:obsdata$nstr, function(x){
          ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), NA)
        })
        seKappa <- sapply(1:obsdata$nstr, function(x){
          ifelse(var[nFpar + obsdata$nstr + x] > 0,
                 sqrt(var[nFpar + obsdata$nstr + x] * kappa[x]^2), NA)
        })
        STDERR <- c(seAlpha=seAlpha, seKappa=seKappa)
      }
      
      #regression parameter(s)
      if (nRpar == 0) {
        seBeta <- NULL
      } else {
        seBeta <- numeric(nRpar)
        varBeta <- var[-(1:(nFpar + nBpar * obsdata$nstr))]
        for(i in 1:nRpar) {
          seBeta[i] <- ifelse(varBeta[i] > 0, sqrt(varBeta[i]), NA)
        }
        PVAL <- c(rep(NA, nFpar + nBpar * obsdata$nstr), 
                  2 * pnorm(q= -abs(beta / seBeta)))
      }
      
      #all together
      STDERR <- c(STDERR,
                  se.beta=seBeta)
      if (frailty != "none") {
        STDERR <- c(se.theta=seTheta,
                    se.sigma2=seSigma2,
                    se.nu=seNu,
                    STDERR)
      }
    }
  }
  
  #----- Output ---------------------------------------------------------------#
  resmodel <- cbind(ESTIMATE=ESTIMATE, SE=STDERR)
  rownames(resmodel) <- gsub("beta.","", rownames(resmodel))
  
  if (nRpar > 0) {
    resmodel <- cbind(resmodel, "p-val"= PVAL)
  }
  
  class(resmodel) <- c("parfm", class(resmodel))
  ### - Checks - ###############################################################
  Call <- match.call()
  if (!match("formula", names(Call), nomatch=0))
    stop("A formula argument is required")
  
  Terms <- terms(formula, data = data)
  ######################################################## - End of Checks - ###
  attributes(resmodel) <- c(attributes(resmodel), list(
    call        = Call,
    convergence = res$convergence,
    it          = it,
    extime      = extime,
    nobs        = nrow(data),
    shared      = (nrow(data) > obsdata$ncl),
    loglik      = lL,
    dist        = dist,
    cumhaz      = attributes(Mloglikelihood(p=res$par,
                                            obs=obsdata, dist=dist, 
                                            frailty=frailty,
                                            correct=correct))$cumhaz,
    cumhazT     = attributes(Mloglikelihood(p=res$par,
                                            obs=obsdata, dist=dist, 
                                            frailty=frailty,
                                            correct=correct))$cumhazT,
    di          = obsdata$di,
    dq          = obsdata$dq,
    dqi         = obsdata$dqi,
    frailty     = frailty,
    clustname   = cluster,
    stratname   = strata,
    correct     = correct,
    formula     = as.character(Call[match("formula", names(Call), nomatch=0)]),
    terms       = attr(Terms, "term.labels")
  ))
  if (frailty != "none") {
    names(attr(resmodel, "cumhaz")) <-
      names(attr(resmodel, "di")) <- unique(obsdata$cluster)
  }
  if (showtime){
    cat("\nExecution time:", extime, "second(s) \n")
  }
  
  return(resmodel)
}