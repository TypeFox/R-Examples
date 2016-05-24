## this is a  Penalised B-spline 
## Paul Eilers and Mikis Stasinopoulos
# needs prediction
#----------------------------------------------------------------------------------------
#source('~/Documents/gamlss/projects/TESTS/PenRegQ.R')
pbq<-function(x, control=pbq.control(...), ...) 
{
# ---------------------------------------------------
            scall <- deparse(sys.call())
#----------------------------------------------------
             rexpr <- regexpr("gamlss",sys.calls())
          for (i in length(rexpr):1)
          { 
          position <- i 
           if (rexpr[i]==1) break
          }
        gamlss.env <- sys.frame(position)
#--------
## get a random name to use it in the gamlss() environment
#--------
                sl <- sample(letters, 4)
       fourLetters <- paste(paste(paste(sl[1], sl[2], sep=""), sl[3], sep=""),sl[4], sep="")
   startLambdaName <- paste("start.Lambda",fourLetters, sep=".")
## put the starting values in the gamlss()environment
lambda <- 1
assign(startLambdaName, lambda, envir = gamlss.env)
#--------   
                             xvar <- x #  set x to zero
      attr(xvar, "control")       <- control
      attr(xvar, "call")          <- substitute(gamlss.pbq(data[[scall]], z, w)) 
      attr(xvar, "gamlss.env")    <- gamlss.env
      attr(xvar, "NameForLambda") <- startLambdaName
      attr(xvar, "class")         <- "smooth"
      xvar
}
#----------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
pbq.control <- function(order = 2,  plot = FALSE, ...) # 
{                         
        if(order < 0) {
warning("the value of order supplied is zero or negative the default value of 1 was used instead")
                order <- 1}                      
        list(order = order,  plot=FALSE)
}
#----------------------------------------------------------------------------------------
gamlss.pbq <- function(x, y, w, xeval = NULL, ...)
{
# -------------------------------------------------- 
# the main function starts here
# get the attributes
        control <- as.list(attr(x, "control")) 
     gamlss.env <- as.environment(attr(x, "gamlss.env"))
startLambdaName <- as.character(attr(x, "NameForLambda")) 
          order <- control$order # the order of the penalty matrix
              N <- length(y) # the no of observations
          tau2  <- sig2 <- NULL
# now the action depends on the values of lambda and df
#--------------------------------------------------------------------  
  
     lambda <-  get(startLambdaName, envir=gamlss.env) ## geting the starting value
 # cat("lambda", lambda, "\n")
        if (lambda>=1e+07) lambda <- 1e+07 # MS 19-4-12
        if (lambda<=1e-07) lambda <- 1e-07 # MS 19-4-12  
        fit <- if ((lambda==1e-07)||(lambda==1e+07)) penRegQ(y=y, x=x, weights=w,  lambda=lambda)
              else   penRegQ(y=y, x=x, weights=w,  start=lambda)
        
        lambda <- c(fit$lambda)
         assign(startLambdaName, lambda, envir=gamlss.env)
#  }
        if (is.null(xeval)) # if no prediction 
    {
     list(fitted.values=fitted(fit), residuals=y-fitted(fit), var=fitted(fit), nl.df =fit$edf-2,
          lambda=fit$sig2/fit$tau2, coefSmo=fit )
    }                            
    
}
#-------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
# this function from nlme of Pinheiro and Bates 
# HessianPB<-function (pars, fun, ..., .relStep = (.Machine$double.eps)^(1/3), 
#    minAbsPar = 0) 
#{
#    pars <- as.numeric(pars)
#    npar <- length(pars)
#    incr <- ifelse(abs(pars) <= minAbsPar, minAbsPar * .relStep, 
#        abs(pars) * .relStep)
#    baseInd <- diag(npar)
#    frac <- c(1, incr, incr^2)
#    cols <- list(0, baseInd, -baseInd)
#    for (i in seq_along(pars)[-npar]) {
#        cols <- c(cols, list(baseInd[, i] + baseInd[, -(1:i)]))
#        frac <- c(frac, incr[i] * incr[-(1:i)])
#    }
#    indMat <- do.call("cbind", cols)
#    shifted <- pars + incr * indMat
#    indMat <- t(indMat)
#    Xcols <- list(1, indMat, indMat^2)
#    for (i in seq_along(pars)[-npar]) {
#        Xcols <- c(Xcols, list(indMat[, i] * indMat[, -(1:i)]))
#    }
#    coefs <- solve(do.call("cbind", Xcols), apply(shifted, 2, 
#        fun, ...))/frac
#    Hess <- diag(coefs[1 + npar + seq_along(pars)], ncol = npar)
#    Hess[row(Hess) > col(Hess)] <- coefs[-(1:(1 + 2 * npar))]
#    list(mean = coefs[1], gradient = coefs[1 + seq_along(pars)], 
#    
#        Hessian = (Hess + t(Hess)))
#}
#------------------------------------------------------------------------

