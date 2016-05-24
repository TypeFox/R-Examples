##
## Functions to fit the Common Component Bivariate Model
##
## 
## cov0.pars <- list(sigma, phi0)  sigma is optional in some functions
## cov1.pars <- list(nu1, phi1)
## cov2.pars <- list(eta, nu2, phi2)
##
".dist12" <- function(geodata1, geodata2)
{
  res <- list()
  res$dist1 <- dist(geodata1$coords)
  res$dist2 <- dist(geodata2$coords)
  res$dist12 <- loccoords(geodata1$coords, geodata2$coords)
  return(res)
}

".cov012.model" <-
  function(cov0.model="matern", cov1.model="matern", cov2.model="matern",
           kappa0=0.5, kappa1=0.5, kappa2=0.5)
{
  cov0.model <- match.arg(cov0.model, choices =  .geoR.cov.models)
  cov1.model <- match.arg(cov1.model, choices =  .geoR.cov.models)
  cov2.model <- match.arg(cov2.model, choices =  .geoR.cov.models)
  if(cov0.model == "stable") cov0.model <- "powered.exponential"
  if(cov1.model == "stable") cov1.model <- "powered.exponential"
  if(cov2.model == "stable") cov2.model <- "powered.exponential"
  if(any(c(cov0.model, cov1.model, cov2.model) == "power"))
    stop("parameter estimation for power model is not implemented")
  if(cov0.model=="matern" & kappa0 == 0.5) cov0.model <- "exponential"
  if(cov1.model=="matern" & kappa1 == 0.5) cov1.model <- "exponential"
  if(cov2.model=="matern" & kappa2 == 0.5) cov2.model <- "exponential"
  return(list(cov.model=c(cov0.model=cov0.model,cov1.model=cov1.model,
              cov2.model=cov2.model),
              kappa = c(kappa0=unname(kappa0), kappa1=unname(kappa1),
                kappa2=unname(kappa2))))
}

"varcovBGCCM" <- function(dists.obj, cov0.pars, cov1.pars, cov2.pars,
                         cov0.model = "matern", cov1.model = "matern",
                         cov2.model = "matern",
                         kappa0 = 0.5, kappa1 = 0.5, kappa2 = 0.5,
                         scaled = TRUE, inv = FALSE, det = FALSE)
{
  ##  A = ( E  F )  ; D = H - G E^{-1} F
  ##      ( G  H )
  ##
  ##  A^{-1} = ( E^{-1} (I + F D^{-1} G E^{-1}  -E^{-1} F D^{-1} )
  ##           ( -D^{-1} G E^{-1}                D^{-1}          )
  ##
  ##  det(A) = det(E) * det(D)
  ##
  ##  So, for:
  ##  A = ( E  F )  ; D = H - F' E^{-1} F
  ##      ( F' H )
  ##
  ##  A^{-1} = ( E^{-1} + E^{-1} F D^{-1} F' E^{-1}   -E^{-1} F D^{-1} )
  ##           ( -D^{-1} F' E^{-1}                     D^{-1}          )
  ##
  CM <- .cov012.model(cov0.model=cov0.model,
                      cov1.model= cov1.model, cov2.model=cov2.model,
                      kappa0=kappa0, kappa1=kappa1, kappa2=kappa2)
  ##
  n1 <- nrow(dists.obj$dist12)
  n2 <- ncol(dists.obj$dist12)
  n <- n1+n2
  ##
  S <- matrix(0, nrow=n1+n2, ncol=n1+n2)
  if(scaled) fc <- 0
  else fc <- n * log(cov0.pars$sigma^2)
  ##
  if(!inv){
    S[(1:n1), (1:n1)] <-
      varcov.spatial(dists.lowertri = dists.obj$dist1,
                     cov.pars=cbind(c(1,cov1.pars$nu1^2),
                       c(cov0.pars$phi0, cov1.pars$phi1)),
                     cov.model=CM$cov.model[c("cov0.model","cov1.model")],
                     kappa=CM$kappa[c("kappa0", "kappa1")])$varcov
    S[(n1+1):n, (n1+1):n] <-
      varcov.spatial(dists.lowertri = dists.obj$dist2,
                     cov.pars=cbind(c(cov2.pars$eta^2, cov2.pars$nu2^2),
                       c(cov0.pars$phi0, cov2.pars$phi2)),
                     cov.model=CM$cov.model[c("cov0.model","cov2.model")],
                     kappa=CM$kappa[c("kappa0", "kappa2")])$varcov
    S[(1:n1), (n1+1):n] <-
      cov.spatial(dists.obj$dist12,
                  cov.pars=c(cov2.pars$eta, cov0.pars$phi0),
                  cov.model=CM$cov.model["cov0.model"],
                  kappa=CM$kappa["kappa0"])
    S[(n1+1):n, (1:n1)] <- t(S[(1:n1), (n1+1):n])
    if(det){
      Einv <- varcov.spatial(dists.lowertri = dists.obj$dist1,
                             cov.pars=cbind(c(1,cov1.pars$nu1^2),
                               c(cov0.pars$phi0, cov1.pars$phi1)),
                             cov.model=CM$cov.model[c("cov0.model",
                               "cov1.model")],
                             kappa=CM$kappa[c("kappa0", "kappa1")],
                             sqrt.inv=TRUE, det=TRUE)
      attr(S, "logdetS") <- fc + (2*Einv$log.det.to.half) +
        determinant(S[(n1+1):n, (n1+1):n] -
                crossprod(crossprod(Einv$sqrt,S[(1:n1),(n1+1):n])))$modulus
    }
    if(scaled)
      return(S)
    else
      return(S*cov0.pars$sigma^2)
  }
  else{
    Einv <- varcov.spatial(dists.lowertri = dists.obj$dist1,
                           cov.pars=cbind(c(1,cov1.pars$nu1^2),
                             c(cov0.pars$phi0, cov1.pars$phi1)),
                           cov.model=CM$cov.model[c("cov0.model","cov1.model")],
                           kappa=CM$kappa[c("kappa0", "kappa1")],
                           inv=TRUE, sqrt.inv=TRUE, det=det)
    H <- varcov.spatial(dists.lowertri = dists.obj$dist2,
                        cov.pars=cbind(c(cov2.pars$eta^2, cov2.pars$nu2^2),
                          c(cov0.pars$phi0, cov2.pars$phi2)),
                        cov.model=CM$cov.model[c("cov0.model","cov2.model")],
                        kappa=CM$kappa[c("kappa0", "kappa2")])$varcov
    F <- cov.spatial(dists.obj$dist12,
                     cov.pars=c(cov2.pars$eta, cov0.pars$phi0),
                     cov.model=CM$cov.model["cov0.model"],
                     kappa=CM$kappa["kappa0"])
    FpEinv <- crossprod(F,Einv$inverse)
    Dinv <- solve(H - FpEinv %*% F)
    S[(n1+1):n,(n1+1):n] <- Dinv
    S[(n1+1):n, (1:n1)] <- - Dinv %*% FpEinv
    S[1:n1, (n1+1):n] <- t(S[(n1+1):n, (1:n1)])
    S[1:n1,1:n1] <- Einv$inverse + crossprod(-S[(n1+1):n,(1:n1)], FpEinv)
    #if(det){
    #  detDinv <- det(Dinv)
    #  if(detDinv <= 0) attr(S, "logdetS") <- NaN
    #  else attr(S, "logdetS") <- fc + (2*Einv$log.det.to.half) - log(detDinv)
    #}
    if(det) attr(S, "logdetS") <- fc + (2*Einv$log.det.to.half) - determinant(Dinv)$modulus
    ## returning S^{-1}
    if(scaled)
      return(S)
    else
      return(S/cov0.pars$sigma^2)
  }
}

".negloglikBGCCM" <- 
  function(pars, ...)
{
  return(-loglikBGCCM(pars, ...))
}

"loglikBGCCM" <-
  function(pars, geodata1, geodata2, cov.model, kappa, envir, ...)
{
  ## pars = c(eta, nu1, nu2, phi0, phi1, phi2)
  ##
  # print(pars)
  if(any(pars[-1] < 0)) return(-.Machine$double.xmax^0.5)
#  print("-----------------------------")
  if(length(cov.model) != 1 & length(cov.model) != 3)
    stop("cov.model must have length 1 or 3")
  if(length(kappa) != 1 & length(kappa) != 3)
    stop("cov.model must have length 1 or 3")  
  if(missing(envir)) envir=sys.frame(sys.nframe())
  calcs <- mget(c("X","y","dist12"), envir=envir, ifnotfound=list(NULL))
  if(is.null(calcs$y))
    calcs$y <- c(geodata1$data, geodata2$data)
  n <- length(calcs$y)
  if(is.null(calcs$X))
    calcs$X <- rbind(kronecker(t(1:0), rep(1,nrow(geodata1$coords))),
                     kronecker(t(0:1), rep(1,nrow(geodata2$coords))))
  if(is.null(calcs$dist12))
    calcs$dist12 <- .dist12(geodata1, geodata2)
  Sinv <- varcovBGCCM(dists.obj=calcs$dist12,
                      cov0.pars=list(phi0=pars[4]),
                      cov1.pars=list(nu1=pars[2],phi1=pars[5]),
                      cov2.pars=list(eta=pars[1],nu2=pars[3],phi2=pars[6]),
                      cov0.model = cov.model[1], cov1.model = cov.model[2],
                      cov2.model = cov.model[3],
                      kappa0 = kappa[1], kappa1 = kappa[2], kappa2 = kappa[3],
                      scaled = TRUE, inv = TRUE, det=TRUE)
  # if(is.nan(attributes(Sinv)$logdetS)) return(-.Machine$double.xmax^0.5)
  SinvX <- crossprod(Sinv,calcs$X)
  QF <- sum(crossprod(calcs$y,Sinv)*calcs$y) - (crossprod(calcs$y, SinvX) %*%
     solve(crossprod(SinvX,calcs$X), crossprod(SinvX,calcs$y)))
  loglik <- drop(-0.5*(attributes(Sinv)$logdetS + n*(log(2*pi) + 1- log(n)+log(QF))))
#  if(is.nan(attributes(Sinv)$logdetS) || QF < 0 || is.nan(loglik)){  
#    print(c(pars, loglik))
#  }
  return(loglik)
}

"likfitBGCCM" <-
  function(geodata1, geodata2, ini.sigmasq, ini.phi,
           cov0.model="matern", cov1.model="matern", cov2.model="matern",
           kappa0=0.5, kappa1=0.5, kappa2=0.5,
           fc.min = c("optim", "nlminb"), ...)
{
  ##
  ## Model:
  ##    Y_1 = mu_1 + S_{01}(x) + S_1(x)
  ##        = mu_1 + \sigma_{01} R_0(x) + \sigma_1 R_1(x)
  ##    Y_2 = mu_2 + S_{02}(x) + S_2(x)
  ##        = mu_2 + \sigma_{02} R_0(x) + \sigma_2 R_1(x)
  ##
  ## Reparametrization
  ##     \sigma = \sigma_{01}  , \eta = \sigma_{02}/\sigma
  ##     \nu1 = \sigma_1/\sigma, \nu2 = \sigma_2/\sigma
  ##
  ## optimization: using "concentrated likelihood" where
  ## \mu_1, \mu_2, \sigma  are given analitically
  ## and the remaning parameters through numerical optimization 
  ## pars = c(eta, nu1, nu2, phi0, phi1, phi2)
  ##
  call.fc <- match.call()
  fc.min <- match.arg(fc.min)
  name.geodata1 <- deparse(substitute(geodata1))
  name.geodata2 <- deparse(substitute(geodata2))
  CM <- .cov012.model(cov0.model=cov0.model,
                      cov1.model= cov1.model, cov2.model=cov2.model,
                      kappa0=kappa0, kappa1=kappa1, kappa2=kappa2)
  ##
  ##
  ##
  n1 <- nrow(geodata1$coords)
  n2 <- nrow(geodata2$coords)
  ## ini.sigmasq = c(sigmasq0.1, sigmasq1, sigmasq0.2, sigmasq2)
  ## ini.phi     = c(phi0, phi1, phi2)
  ##
  ## checking initial values
  ##
  if(!missing(ini.sigmasq) && length(ini.sigmasq) != 4)
    stop("ini.sigmasq must have length 4")
  if(!missing(ini.phi) && length(ini.phi) != 3)
    stop("ini.phi must have length 3")
  ##
  ## computing default initial values
  ##
  if(missing(ini.sigmasq))
    ini.sigmasq <- c(rep(0.5*var(geodata1$data),2),
                     rep(0.5*var(geodata2$data),2))
  if(missing(ini.phi)){
    m1 <- max(dist(geodata1$coords))
    m2 <- max(dist(geodata2$coords))
    ini.phi <- c(0.1*max(c(m1,m2)), 0.15*m1, 0.15*m2)
  }
  ##
  ## reparametrising initial values for
  ##
  eta <- sqrt(ini.sigmasq[3]/ini.sigmasq[1])
  nu1 <- sqrt(ini.sigmasq[2]/ini.sigmasq[1])
  nu2 <- sqrt(ini.sigmasq[4]/ini.sigmasq[1])
  ini.pars = c(eta, nu1, nu2, ini.phi)
  ##
  ## storing useful results to prevent unnecessary (re)calculations
  ##
  likBGCCM.env <- new.env()
  assign("dist12", .dist12(geodata1, geodata2), envir=likBGCCM.env)
  assign("X", rbind(kronecker(t(1:0), rep(1,n1)),
                    kronecker(t(0:1), rep(1,n2))),
         envir=likBGCCM.env)
  assign("y", c(geodata1$data, geodata2$data), envir=likBGCCM.env)
  ##
  ## obtaining numerical estimates
  ##
  ldots <- list(...)
  if(!is.null(names(ldots))){
    names(ldots)[which(as.logical(pmatch(names(ldots), "lower", nomatch=0)))] <- "lower"
    names(ldots)[which(as.logical(pmatch(names(ldots), "method", nomatch=0)))] <- "method"
  }
  if(fc.min == "optim"){
    if(!is.null(ldots$method) && ldots$method == "L-BFGS-B" && is.null(ldots$lower))
      ldots$lower <- c(-Inf, rep(0,5))
    est <- do.call("optim", c(list(par=ini.pars, fn=.negloglikBGCCM,
                                   geodata1 = geodata1, geodata2 = geodata2,
                                   cov.model=CM$cov.model, kappa=CM$kappa,
                                   envir=likBGCCM.env), ldots))
  }
  else{
    if(is.null(ldots$lower)) ldots$lower <- rep(0,6)
    est <- do.call("nlminb", c(list(start=ini.pars, objective=.negloglikBGCCM,
                                    geodata1 = geodata1, geodata2 = geodata2,
                                    cov.model=CM$cov.model, kappa=CM$kappa,
                                    envir=likBGCCM.env), ldots))
    est$value <- est$objective
  }
  ##
  ## computing means and basic variance sigma^2
  ##
  par0 <- list(phi0=est$par[4])
  par1 <- list(nu1=est$par[2], phi1=est$par[5])
  par2 <- list(eta=est$par[1], nu2=est$par[3], phi2=est$par[6])
  calcs <- mget(c("X","y","dist12"), envir=likBGCCM.env,
                ifnotfound=list(NULL))
  Sinv <- varcovBGCCM(calcs$dist12,
                      cov0.pars = par0, cov1.pars=par1, cov2.pars=par2,
                      cov0.model=cov0.model,
                      cov1.model= cov1.model, cov2.model=cov2.model,
                      kappa0=kappa0, kappa1=kappa1, kappa2=kappa2,
                      scaled = TRUE, inv = TRUE)
  SinvX <- crossprod(Sinv,calcs$X)
  mu12 <- solve(crossprod(calcs$X,SinvX), crossprod(SinvX,calcs$y))
  yminusXB <- drop(calcs$y - calcs$X %*% mu12)
  sigma2 <- sum(crossprod(yminusXB,Sinv)*yminusXB)/length(calcs$y)
  ##
  ## preparing output
  ##
  res=list()
  res$mu <- drop(mu12)
  names(res$mu) <- c("mu1", "mu2")
  ## converting back to original parametrisation
  res$sigmasq <- c(sigmasq01 = sigma2, sigmasq1=(sigma2)*(par1$nu1^2),
                   sigmasq02 = (sigma2)*(par2$eta^2),
                   sigmasq2=(sigma2)*(par2$nu2^2))
  res$phi <- c(phi0=par0$phi0, phi1=par1$phi1, phi2=par2$phi2)
  res$cov0.pars <- c(sigma=sqrt(sigma2), par0)
  res$cov1.pars <- par1
  res$cov2.pars <- par2
  res$loglik <- (- est$value)
  res$cov.model <- CM$cov.model
  res$kappa <- CM$kappa
  res$practicalRange <- c(pr0=with(res, practicalRange(cov.model[1], phi[1], kappa[1])),
                          pr1=with(res, practicalRange(cov.model[2], phi[2], kappa[2])),
                          pr2=with(res, practicalRange(cov.model[3], phi[3], kappa[3])))
  res$n <- c(n1=n1, n2=n2)
  res$optim <- est
##  res$y <- calcs$y
  res$Sinv <- Sinv
  names(res$optim$par) <- c("eta","nu1","nu2","phi0","phi1","phi2")
  res$call <- call.fc
  attr(res, "geodata1") <- name.geodata1
  attr(res, "geodata2") <- name.geodata2
  oldClass(res) <- "BGCCM"
  return(res)
} 

"print.BGCCM" <-
  function(x, ...)
{
  cat("Mean parameters:\n")
  print(format(x$mu, ...))
  cat("Variance parameters:\n")
  print(format(x$sigmasq, ...))
  cat("Correlation parameters:\n")
  print(format(x$phi, ...))
  cat("Extra correlation parameter (fixed):\n")
  print(format(x$kappa, ...))
  cat("Practical Ranges with cor=0.05 for asymptotic range:\n")
  print(format(x$practicalRange, ...))
  cat("Reparametrised:\n")
  print(format(x$optim$par, ...))
  cat(paste("Maximised log-Likelihood:",x$loglik,"\n"))
  return(invisible())
}

".naiveLL.BGCCM" <-
  function(geodata1, geodata2, mu, cov0.pars, cov1.pars, cov2.pars,
           cov0.model="matern", cov1.model="matern", cov2.model = "matern",
           kappa0 = 0.5, kappa1 = 0.5, kappa2 = 0.5, profile=FALSE)
{
  if(!profile && missing(mu))
    stop("mu is needed if profile=FALSE")
  if(!profile && is.null(cov0.pars$sigma))
    stop("cov0.pars$sigma is needed if profile=FALSE")
  y <- c(geodata1$data, geodata2$data)
  n <- length(y)
  C <- varcovBGCCM(.dist12(geodata1,geodata2), cov0.pars=cov0.pars,
                   cov1.pars=cov1.pars, cov2.pars=cov2.pars,
                   cov0.model=cov0.model, cov1.model=cov1.model,
                   cov2.model=cov2.model, kappa0=kappa0,
                   kappa1=kappa1, kappa2=kappa2, scaled=TRUE)
  X <- cbind(rep(c(1,0), c(length(geodata1$data),length(geodata2$data))),
             rep(c(0,1), c(length(geodata1$data),length(geodata2$data))))
  if(profile){
    iC <- solve(C)
    iCX <- crossprod(iC,X)
    mu <- solve(crossprod(iCX,X),crossprod(iCX,y))
  }
  yXb <- drop(y - X %*% mu)
  Q <- sum(crossprod(yXb, solve(C))*yXb)
  if(profile)
    cov0.pars$sigma <- sqrt(Q/n)
  loglik <- -0.5 * (log(2*pi) + n*log(cov0.pars$sigma^2) +
                    determinant(C)$modulus + Q/(cov0.pars$sigma^2))
  return(loglik)
}

"predict.BGCCM" <-
  function(object, locations, borders, variable.to.predict=1, ... )
{
  call.fc <- match.call()
  ##  if(!is.vector(variable.to.predict) || !is.numeric(variable.to.predict))
  ##    stop("variable must be a scalar equals to 1 or 2,
  ## or a vector with elements 1:2")
  ##  if(any(is.na(match(variable.to.predict,1:2))))
  ##    stop("variable must be a scalar equals to 1 or 2,
  ## or a vector with elements 1:2")
  if(!is.vector(variable.to.predict) ||
     !is.numeric(variable.to.predict) || length(variable.to.predict) != 1)
    stop("variable must be a scalar equals to 1 or 2")
  ##
  locations <- .check.locations(locations)
  ##
  geodata1 <-  get(attr(object, "geodata1"))
  geodata2 <-  get(attr(object, "geodata2"))
  ##
  ## selecting locations inside the borders 
  ##
  if(missing(borders)){
    if(variable.to.predict==1) borders <- geodata1$borders
    if(variable.to.predict==2) borders <- geodata2$borders
  }
  if(!is.null(borders)){
    locations <- locations[.geoR_inout(locations, borders),]
    if(nrow(locations) == 0)
      stop("there are no prediction locations inside the borders")
#    if(messages.screen)
#      cat("results will be returned only for prediction locations inside the borders\n")
  }
  dimnames(locations) <- list(NULL, NULL)
  ##
  n0 <- nrow(locations)
  n1 <- object$n[1]
  n2 <- object$n[2]
  n <- n1+n2
  ##
  res12 <- c((geodata1$data-object$mu[1]),(geodata2$data-object$mu[2]))
  ##
  S012 <- matrix(0, nrow = length(res12), ncol=nrow(locations))
  names(object$mu) <- NULL
  if(any(variable.to.predict == 1)){
    mu0 <- rep(object$mu[1], n0)
    var0 <- sum(object$sigmasq[c("sigmasq01","sigmasq1")])
    S012[1:n1,] <-
      cov.spatial(loccoords(geodata1$coords, locations),
                  cov.pars=cbind(c(1,object$cov1.pars$nu1^2),
                    c(object$cov0.pars$phi0, object$cov1.pars$phi1)),
                  cov.model=object$cov.model[c("cov0.model", "cov1.model")],
                  kappa=object$kappa[c("kappa0", "kappa1")])
    S012[(n1+1):(n1+n2),] <-
      cov.spatial(loccoords(geodata2$coords, locations),
                  cov.pars=c(object$cov2.pars$eta,
                    object$cov0.pars$phi0),
                  cov.model=object$cov.model["cov0.model"],
                  kappa=object$kappa["kappa0"])
  }
  if(any(variable.to.predict == 2)){
    mu0 <- rep(object$mu[2], n0)
    var0 <- sum(object$sigmasq[c("sigmasq02","sigmasq2")])
    S012[1:n1,] <-
      cov.spatial(loccoords(geodata1$coords, locations),
                  cov.pars=c(object$cov2.pars$eta,
                    object$cov0.pars$phi0),
                  cov.model=object$cov.model["cov0.model"],
                  kappa=object$kappa["kappa0"])
    S012[(n1+1):(n1+n2),] <-
      cov.spatial(loccoords(geodata2$coords, locations),
                  cov.pars=cbind(c(object$cov2.pars$eta^2,
                    object$cov2.pars$nu2^2),
                    c(object$cov0.pars$phi0, object$cov2.pars$phi2)),
                  cov.model=object$cov.model[c("cov0.model", "cov2.model")])
  }
  ##
  res <- list()
  res$predict <- mu0 + drop(crossprod(S012,object$Sinv)%*%res12)
  res$krige.var <- var0 -
    .diagquadraticformXAX(S012,lowerA=object$Sinv[lower.tri(object$Sinv)],
                          diagA= diag(object$Sinv))
  res$call <- call.fc
  attr(res, "prediction.locations") <- call.fc$locations
  if(!is.null(call.fc$borders))
    attr(res, "borders") <- call.fc$borders
  oldClass(res) <- c("BGCCMpred", "kriging")
  return(res)
}

