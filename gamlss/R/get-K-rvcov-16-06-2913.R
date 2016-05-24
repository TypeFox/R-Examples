
#-------------------------------------------------------------------------------
# What I need is
# i)    y and fitted mu sigma nu and tau
# ii)   the X_i i=1,2,3,4 matrices
# iii)  the dl/dmu and dmu/deta for all parameters
# TO DO
# i)  shall we use the variance of dl/db or calculate it ourself
# ii) with what shall we devide?
# iii) do we have to multiply with weights"???
#-------------------------------------------------------------------------------
get.K <- function(object, what=c("K", "Deriv") )# a gamlss object
{
if (!is.gamlss(object)) stop("needs a gamlss object")
        what <-  match.arg(what)
         fam <-  if(is.null(object$call$family)) as.gamlss.family(NO) 
                 else as.gamlss.family(object$call$family)
       fname <- object$family[1]
       #dfun <- paste("d",fname,sep="")
       #mpdf <- eval(parse(text=dfun))
       # cdf <- eval(parse(text=pfun)) # may be we need this for censored data
       nopar <- length(object$par)
# get y and weights
          y  <- object$y
          w  <- object$weights 
          X  <- list()
          lp <- list() # predictors
        dldp <- list() # first derivarives  
        dpdh <- list() # deta/dmu
# get the parameters
  if ("mu"%in%object$par)       mu <- fitted(object, "mu") 
  if ("sigma"%in%object$par) sigma <- fitted(object, "sigma") 
  if ("nu"%in%object$par)       nu <- fitted(object, "nu") 
  if ("tau"%in%object$par)     tau <- fitted(object, "tau") 
   for (i in object$par)
   {
     if (i=="mu")    {ii <- "m"} # for the first deri
     if (i=="sigma") {ii <- "d"}
     if (i=="nu")    {ii <- "v"}  
     if (i=="tau")   {ii <- "t"}
 #dPar[[i]] <- eval(parse(text=paste(paste("object$", i, sep=""), ".fv", sep="")))
     X[[i]] <- eval(parse(text=paste(paste("object$", i, sep=""), ".x", sep="")))
    lp[[i]] <- eval(parse(text=paste(paste("object$", i, sep=""), ".lp", sep="")))
  dldp[[i]] <- eval(parse(text=paste("fam$dld", ii, sep="")), envir=environment()) # functions
  dpdh[[i]] <- eval(parse(text=paste(paste("fam$", i, sep=""), ".dr", sep="")))    # functions
# we need to get the derivative dmu/deta
# get the argument right and evaluate them     
#     eval(parse(text=paste("fam$dld", ii, sep="")), envir=environment())
     #dr <- f$dr(eta)
#     dr <- 1/dr     
#     wt <- -(d2ldp2/(dr*dr))# 
#     # wv <- (eta-os)+step*dldp/(dr*wt)
#     wv <- (eta-os)+dldp/(dr*wt)
#     eval(dummy <- Family.d(family = fname, type = type, ...))
#  eval(call("<-","invlink",dummy), envir=environment()) 
#     eval(parse(text=paste("fam$dld", ii, sep="")))
#    OffSet <- if (exitData) model.extract(model.frame(eval(parse(text=paste(paste("object$", 
# i, sep=""), ".terms", sep=""))), data=DaTa), "offset")
#               else  model.extract(model.frame(eval(parse(text=paste(paste("object$", i, 
# sep=""), ".terms", sep="")))), "offset")
#    offSet[[i]] <- if (is.null(OffSet)) rep(0, length(y)) else  OffSet 
   } 
    switch(nopar,
     { #  1 parameter
        dldm <- dldp[["mu"]](y=y, mu=mu)
       dmudh <- dpdh[["mu"]](eta=lp[["mu"]])
         wmu <- as.vector(dldm*dmudh*sqrt(w))
         Xmu <- wmu*X[["mu"]]
        AllX <- cbind(Xmu)
        #  K <- var(AllX)*length(y)
           K <- crossprod(AllX)*length(y)/(length(y)-dim(AllX)[2])
        #  K <- t(AllX)%*%AllX/(length(y)-dim(AllX)[2]) 
     },
     { # 2 parameter
        dldm <- dldp[["mu"]](y=y, mu=mu, sigma=sigma)
        dldd <- dldp[["sigma"]](y=y, mu=mu, sigma=sigma)
       dmudh <- dpdh[["mu"]](eta=lp[["mu"]])
    dsigmadh <- dpdh[["sigma"]](eta=lp[["sigma"]])
         wmu <- as.vector(dldm*dmudh*sqrt(w))
      wsigma <- as.vector(dldd*dsigmadh*sqrt(w))
         Xmu <- wmu*X[["mu"]]
      Xsigma <- wsigma*X[["sigma"]]
        AllX <- cbind(Xmu, Xsigma)
      #    K <- var(AllX) *length(y)
           K <- crossprod(AllX)*length(y)/(length(y)-dim(AllX)[2])
         # K <- t(AllX)%*%AllX/(length(y)-dim(AllX)[2]) # do we divide? by n-k
        },         
      { # 3 parameter 
        dldm <- dldp[["mu"]](y=y,    mu=mu, sigma=sigma, nu=nu)
        dldd <- dldp[["sigma"]](y=y, mu=mu, sigma=sigma, nu=nu)
        dldn <- dldp[["nu"]](y=y,    mu=mu, sigma=sigma, nu=nu)
       dmudh <- dpdh[["mu"]](eta=lp[["mu"]])
    dsigmadh <- dpdh[["sigma"]](eta=lp[["sigma"]])
       dnudh <- dpdh[["nu"]](eta=lp[["nu"]])
         wmu <- as.vector(dldm*dmudh*sqrt(w))
      wsigma <- as.vector(dldd*dsigmadh*sqrt(w))
         wnu <- as.vector(dldn*dnudh*sqrt(w))
         Xmu <- wmu*X[["mu"]]
      Xsigma <- wsigma*X[["sigma"]]
         Xnu <- wnu*X[["nu"]]
        AllX <- cbind(Xmu, Xsigma, Xnu)
        K <- crossprod(AllX)*length(y)/(length(y)-dim(AllX)[2])
        #   K <- t(AllX)%*%AllX/(length(y)-dim(AllX)[2]) # do we divide? by n-k
       },
      { # 4 parameter
        dldm <- dldp[["mu"]](y=y,    mu=mu, sigma=sigma, nu=nu, tau=tau)
        dldd <- dldp[["sigma"]](y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)
        dldn <- dldp[["nu"]](y=y,    mu=mu, sigma=sigma, nu=nu, tau=tau)
        dldt <- dldp[["tau"]](y=y,   mu=mu, sigma=sigma, nu=nu, tau=tau)
       dmudh <- dpdh[["mu"]](eta=lp[["mu"]])
    dsigmadh <- dpdh[["sigma"]](eta=lp[["sigma"]])
       dnudh <- dpdh[["nu"]](eta=lp[["nu"]])
      dtaudh <- dpdh[["tau"]](eta=lp[["tau"]])  
         wmu <- as.vector(dldm*dmudh*sqrt(w))
      wsigma <- as.vector(dldd*dsigmadh*sqrt(w))
         wnu <- as.vector(dldn*dnudh*sqrt(w))
        wtau <- as.vector(dldt*dtaudh*sqrt(w))
         Xmu <- wmu*X[["mu"]]
      Xsigma <- wsigma*X[["sigma"]]
         Xnu <- wnu*X[["nu"]]
        Xtau <- wtau*X[["tau"]]
        AllX <- cbind(Xmu, Xsigma, Xnu, Xtau)
        K <- crossprod(AllX)*length(y)/(length(y)-dim(AllX)[2])
        #   K <- t(AllX)%*%AllX/(length(y)-dim(AllX)[2]) # do we divide? by n-k
      })
if (what=="K") K else AllX
}
#-------------------------------------------------------------------------------
rvcov <- function(object, type = c("vcov", "cor", "se", "coef", "all"),  hessian.fun = c("R", "PB"))
{
    type <- match.arg(type)
    hessian.fun <- match.arg( hessian.fun)
  if (!("gamlss" %in% class(object))) 
    stop("the null model is not a gamlss model")
        V <- vcov(object, hessian.fun= hessian.fun)
        K <- get.K(object)
   varCov <- V%*%K%*%V
      se <- sqrt(diag(varCov))
    corr <- cov2cor(varCov)
  #coefBeta <- unlist(coefBeta)
  switch(type, vcov = varCov, cor = corr, se = se,# coef = coefBeta, 
         all = list( se = se, vcov = varCov, cor = corr)) # varCov
}