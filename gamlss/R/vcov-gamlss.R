#  Mikis Stasinopoulos 22-03-13
# This is the new vcov.gamlss() function 
# it was created in Jan 2013 but impemented in gamlss on the 22 March 2013  
# it is using the function gen.likelihood

#-------------------------------------------------------------------------------
vcov.gamlss <- function (object, 
                           type = c("vcov", "cor", "se", "coef", "all"),
                         robust = FALSE, 
                    hessian.fun = c("R", "PB"),   
                ...) 
{
## local function -------------------------------------------------------------
  HessianPB<-function (pars, fun, ..., .relStep = (.Machine$double.eps)^(1/3), 
                       minAbsPar = 0) 
  {
    pars <- as.numeric(pars)
    npar <- length(pars)
    incr <- ifelse(abs(pars) <= minAbsPar, minAbsPar * .relStep, 
                   abs(pars) * .relStep)
    baseInd <- diag(npar)
    frac <- c(1, incr, incr^2)
    cols <- list(0, baseInd, -baseInd)
    for (i in seq_along(pars)[-npar]) {
      cols <- c(cols, list(baseInd[, i] + baseInd[, -(1:i)]))
      frac <- c(frac, incr[i] * incr[-(1:i)])
    }
    indMat <- do.call("cbind", cols)
    shifted <- pars + incr * indMat
    indMat <- t(indMat)
    Xcols <- list(1, indMat, indMat^2)
    for (i in seq_along(pars)[-npar]) {
      Xcols <- c(Xcols, list(indMat[, i] * indMat[, -(1:i)]))
    }
    coefs <- solve(do.call("cbind", Xcols), apply(shifted, 2, 
                                                  fun, ...))/frac
    Hess <- diag(coefs[1 + npar + seq_along(pars)], ncol = npar)
    Hess[row(Hess) > col(Hess)] <- coefs[-(1:(1 + 2 * npar))]
    list(mean = coefs[1], gradient = coefs[1 + seq_along(pars)], 
         
         Hessian = (Hess + t(Hess)))
  }  
## end of local ----------------------------------------------------------------  
       type <- match.arg(type)
hessian.fun <- match.arg(hessian.fun)
  if (!is.gamlss(object)) 
     stop(paste("This is not an gamlss object", "\n", ""))
  coefBeta <- list()
  for (i in object$par) 
  {
    if (i == "mu") 
      {
      if (!is.null(unlist(attr(terms(formula(object), specials = .gamlss.sm.list), 
                               "specials")))) 
        warning("Additive terms exists in the mu formula. \n "
                ,"Standard errors for the linear terms maybe are not appropriate")
    }
    else 
    {
      if (!is.null(unlist(attr(terms(formula(object, i), 
                                     specials = .gamlss.sm.list), "specials")))) 
        warning(paste("Additive terms exists in the ", i, "formula. \n "
                ,"Standard errors for the linear terms maybe are not appropriate"))
    }
  #    parname <- paste(i, "start", sep = ".")
    nonNAcoef <- !is.na(coef(object, i))
     coefBeta <- c(coefBeta, coef(object, i)[nonNAcoef])
  }
   betaCoef <- unlist(coefBeta)      
   like.fun <- gen.likelihood(object)
## we have a problem here if the likelihood has a lot of parameter for example 
## a lot of factors with large number of levels as in BAT data where 
##  system.time(H <- HessianPB(betaCoef, like.fun))
##      user  system elapsed 
##      25.875   6.076  37.473 
##  system.time(optimHess(betaCoef, like.fun))
##      user  system elapsed 
##      194.455  48.060 242.107
## I think we should have a option to be able to use       
       hess <- if (hessian.fun=="R" ) optimHess(betaCoef, like.fun)
               else HessianPB(betaCoef, like.fun)$Hessian
     varCov <- try(solve(hess), silent = TRUE)
      if (any(class(varCov)%in%"try-error")||any(diag(varCov)<0))
        { # now if it fails try the "PB" function MS 9-2-14
        varCov <- try(solve(HessianPB(betaCoef, like.fun)$Hessian), silent = TRUE) 
        if (any(class(varCov)%in%"try-error")) # if it still fail give up
        stop("the Hessian matrix is singular probably the model is overparametrised")
        }
     rownames(varCov) <- colnames(varCov) <- rownames(hess)
         se <- sqrt(diag(varCov))
       corr <- cov2cor(varCov) # cov/(se %o% se)
   coefBeta <- unlist(coefBeta)
  #names(coefBeta) <- attr(a$se, "names")
    if (robust)
    {
      K <- get.K(object)
      varCov <- varCov%*%K%*%varCov
      se <- sqrt(diag(varCov))
      corr <- cov2cor(varCov)
    }
  switch(type, vcov = varCov, cor = corr, se = se, coef = coefBeta, 
         all = list(coef = coefBeta, se = se, vcov = varCov, 
                    cor = corr))
}