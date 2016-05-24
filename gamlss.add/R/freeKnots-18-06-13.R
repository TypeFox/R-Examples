#----------------------------------------------------------------------------------------
# for fixed knots it fits a piecewise beta spline
# FUNCTIONS
#  i) fitFixedKnots()
# ii) Methods for FixBreakPointsReg objects
#      1) print()
#      2) fitted()
#      3) redid()
#      4) coef()
#      5) knots()
#      6) predict()
#  iii) fitFreeKnots()
#   iv) Methods for FreeBreakPointsReg
#       1) print()
#       2) fitted()
#       3) resid()
#       4) coef()
#       5) summary()
#       6) vcov()
# created by MS Tuesday, July 7, 2009 
# revised by MS Thusday, Aug 18, 2011
# the truncted basis is introdused  plus the Hessian so standard errors can be displaied for 
# "FreeBreakPointsReg" objects
# TO DO:
# i) no summary function for c("FixBreakPointsReg")
# ii) maybe need genoud for more complicated problems
# iii) use optimHess()
# iv) nlminb() instead of optim?
# v) change names to be consistent
# functions
# i)   fitFixedKnots
# ii)
#--------------------------------------------------------------------------------
fitFixedKnots <- function(y, x,  
                   weights = NULL, 
                     knots = NULL,  
                      data = NULL, 
                    degree = 3, 
                     fixed = NULL, 
                      base = c("trun","Bbase"), ...)
{
# fixed argument is not useful here (since all BP are fixed) unless is used in fitFreeKnots()
#-------------------------------------------------------------------------------
# local functions ---
# A least square fit---
# not used at the moment since lm.wfit is more stable
#ls.wr <- function(y, X, w)
#  {
#            XW <- w * X
#           XWX <- t(XW) %*% X
#          beta <- as.vector(solve(XWX , crossprod(XW, y)))
#            fv <- as.vector(X %*%beta)
#           rss <- sum((w*(y-fv))^2)
#           fit <- list(coefficients = beta, fv=fv, rss=rss, df=length(beta) ) #edf = sum(diag(H)))
#  return(fit)
#  }
#-------------------------------------------------------------------------------
# Bbase base
 get.X.beta <- function(x, knots)  
   {
            # lx <- length(x)
             xl <- min(x)
             xr <- max(x)
           xmax <- xr + 0.01 * (xr - xl)
           xmin <- xl - 0.01 * (xr - xl)
         if (is.null(knots)) warning("No knots (break points) are declared")   
             kn <-  sort(c(xmin, xmax, fixed, knots)) 
       # splineDesign(knots, x, ord = 4, derivs, outer.ok = FALSE)
             X2 <- splineDesign(knots=kn, x, ord= degree + 1, derivs=0 * x, outer.ok=TRUE)
   colnames(X2) <- paste("XatBP",1:dim(X2)[2], sep="")
                   #bbase(x=x, xl=xmin, xr=xmax, knots=knots, deg=degree) # the problem here is the difference between break points
           form <- switch(as.character(degree), "0"=~1, "1"=~x, "2"=~x+I(x^2), "3"=~x+I(x^2)+I(x^3), 
                           "4"=~x+I(x^2)+I(x^3)+I(x^4), "5"=~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)) 
             X1 <- model.matrix(form) # model.matrix(~poly(x,degree)) 
              X <-cbind(X1,X2)
             # Xbe <<- X
    X
    } 
#-------------------------------------------------------------------------------
#-trunc base
 get.X.trun <- function(x, knots)  
   {
         tpower <- function(x, t, p)(x - t) ^ p * (x > t)      
          #  lx <- length(x)
          #    xl <- min(x)
          #   xr <- max(x)
          # xmax <- xr + 0.01 * (xr - xl)
          # xmin <- xl - 0.01 * (xr - xl)
         if (is.null(knots)) warning("No knots (break points) are declared")   
          kn <-  sort(c(fixed, knots))
       # splineDesign(knots, x, ord = 4, derivs, outer.ok = FALSE)
             X2 <- outer(x, kn, tpower, degree)# calculate the power in the knots
   colnames(X2) <- paste("XatBP",1:dim(X2)[2], sep="")
           form <- switch(as.character(degree), "0"=~1, "1"=~x, "2"=~x+I(x^2), "3"=~x+I(x^2)+I(x^3), 
                           "4"=~x+I(x^2)+I(x^3)+I(x^4), "5"=~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)) 
             X1 <- model.matrix(form) # model.matrix(~poly(x,degree)) 
              X <-cbind(X1,X2)
            # Xtr <<- X
    X
    } 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# main function starts here
# if data exit attach them
#       if (is.data.frame(data)) { attach(data); on.exit(detach(data))}
         ylab <- deparse(substitute(y))
         xlab <- deparse(substitute(x))
            y <- if (!is.null(data)) get(deparse(substitute(y)), envir=as.environment(data)) else y
            x <- if (!is.null(data)) get(deparse(substitute(x)), envir=as.environment(data)) else x
    w <- if(is.null(weights))  rep(1,length(y)) 
         else {
                 if(!is.null(data))  get(deparse(substitute(weights)), envir=as.environment(data))
                  else weights
               }
              lx <- length(x)
        whatbase <- match.arg(base)
               X <- if (whatbase=="Bbase") get.X.beta(x, knots) else get.X.trun(x, knots)
             fit <- lm.wfit(X,y,w) 
           sigma <- sum(resid(fit)^2)/length(y) 
    names(sigma) <-"sigma"
            out  <- list(call = sys.call(), 
                fitted.values = fitted(fit), 
                    residuals = resid(fit), 
                           df = fit$rank, 
                          rss = sum(resid(fit)^2),  
                         coef = fit$coefficients, sigma = sigma, 
                        fixed = fixed, 
                  breakPoints = knots,  
                       degree = degree, 
                         base = whatbase, 
                            y = y, 
                            x = x, 
                            w = w, 
                            X = X,  
                           qr = fit$qr,  
                     deviance = sum(-2*dNO(y, mu=fitted(fit), 
                        sigma = sigma)))
      class(out) <- c("FixBreakPointsReg")
         out
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# methods
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# print 
print.FixBreakPointsReg<-function(x, digits=max(3, getOption("digits") - 3), ...) 
{
 cat("\nCall: ", deparse(x$call),  "\n", fill=TRUE)
 #cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("Fixed Knots: \n")
      print.default(format(knots(x), digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
    invisible(x)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# fitted 
fitted.FixBreakPointsReg<-function(object,...) 
{
object$fitted.values
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
residuals.FixBreakPointsReg<-function(object,...) 
{
object$residuals
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
coef.FixBreakPointsReg<-function(object,...) 
{
object$coef
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
knots.FixBreakPointsReg<-function(Fn,...)
{
BP<-Fn$breakPoints
names(BP)<-paste("BP",1:length(Fn$breakPoints), sep="")
BP
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
predict.FixBreakPointsReg<-function(object, newdata=NULL, old.x.range=TRUE,...)
{
#-------------------------------------------------------------------------
 # Bbase
 get.X.beta <- function(x, knots, xmin=NULL, xmax=NULL)  
   {
          #   lx <- length(x)
             xl <- min(x)
             xr <- max(x)
          if (is.null(xmax)) xmax <- xr + 0.01 * (xr - xl) else xmax
          if (is.null(xmin))  xmin <- xl - 0.01 * (xr - xl) else xmin
         if (is.null(knots)) warning("No knots (break points) are declared")   
             kn <-  sort(c(xmin, xmax, fixed, knots)) 
       # splineDesign(knots, x, ord = 4, derivs, outer.ok = FALSE)
             X2 <- splineDesign(knots=kn, x, ord= degree + 1, derivs=0 * x, outer.ok=TRUE)
   colnames(X2) <- paste("XatBP",1:dim(X2)[2], sep="")
                   #bbase(x=x, xl=xmin, xr=xmax, knots=knots, deg=degree) # the problem here is the difference between break points
           form <- switch(as.character(degree), "0"=~1, "1"=~x, "2"=~x+I(x^2), "3"=~x+I(x^2)+I(x^3), 
                           "4"=~x+I(x^2)+I(x^3)+I(x^4), "5"=~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)) 
             X1 <- model.matrix(form) # model.matrix(~poly(x,degree)) 
              X <-cbind(X1,X2)
    X
    } 
#--------------------------
#-trunc base
 get.X.trun <- function(x, knots, xmin=NULL, xmax=NULL)  
   {
         tpower <- function(x, t, p)(x - t) ^ p * (x > t)      
            # lx <- length(x)
             xl <- min(x)
             xr <- max(x)
          if (is.null(xmax)) xmax <- xr + 0.01 * (xr - xl) else xmax
          if (is.null(xmin))  xmin <- xl - 0.01 * (xr - xl) else xmin
          if (is.null(knots)) warning("No knots (break points) are declared")   
             kn <-  sort(c(fixed, knots)) 
       # splineDesign(knots, x, ord = 4, derivs, outer.ok = FALSE)
             X2 <- outer(x, kn, tpower, degree)# calculate the power in the knots
   colnames(X2) <- paste("XatBP",1:dim(X2)[2], sep="")
           form <- switch(as.character(degree), "0"=~1, "1"=~x, "2"=~x+I(x^2), "3"=~x+I(x^2)+I(x^3), 
                           "4"=~x+I(x^2)+I(x^3)+I(x^4), "5"=~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)) 
             X1 <- model.matrix(form) # model.matrix(~poly(x,degree)) 
              X <-cbind(X1,X2)
    X
    } 
#----------------------------------------------------------------------- -------
# function start here 
if (is.null(newdata))  #
    {
    predictor<- fitted(object)
    return(predictor)
    }
    lox <- length(object$x)
      x <- c(object$x, newdata)    
     lx <- length(x)
  knots <- knots(object)
  fixed <- object$fixed
 degree <- object$degree
# here we have two choices
# the first is to create new end points including the newdata
# the second is to only include the old data end points
# If we take the fist choice old.x.range=TRUE
# the prediction could be possible better outside the x range
# but would not coincide with the original predictions i.e. fitted(model)
# here the default includes the old data point range with the consequence that predictions 
# outside the original x-range will linear, quadratic, cubic  etc, depending on the 
# value of the degree of the polynomial fitted
if (old.x.range)
 {# TRUE 
     xl <- min(object$x)
     xr <- max(object$x)
   xmax <- xr + 0.01 * (xr - xl)
   xmin <- xl - 0.01 * (xr - xl)
      X <- if (object$base=="Bbase") get.X.beta(x, knots, xmin-xmin, xmax=xmax) 
           else get.X.trun(x, knots, xmin-xmin, xmax=xmax)
   #X2 <- splineDesign(knots=object$knots, x, ord= object$degree + 1, derivs=0 * x, outer=T)#
 }
 else
 {
   X <-  if (object$base=="Bbase") get.X.beta(x, knots) else get.X.trun(x, knots)
     #splineDesign(knots=kn, x, ord= object$degree + 1, derivs=0 * x, outer=T)#
 }   
# browser()
  #form <- switch(as.character(object$degree), "0"=~1, "1"=~x, "2"=~x+I(x^2), "3"=~x+I(x^2)+I(x^3), 
  #                         "4"=~x+I(x^2)+I(x^3)+I(x^4), "5"=~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)) 
   #X1 <- model.matrix(form) # model.matrix(~poly(x,degree))
    #X <- cbind(X1,X2) 
   fv <- X%*%coef(object) 
   if (any(abs(fv[1:lox]-fitted(object))>1e-005)) 
warning(paste("There is a discrepancy  between the original prediction and the re-fit"," \n used to achieve 'safe' predictions \n ", sep = "" ))  
 pred <-  fv[(lox+1):lx]
 pred
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#****************************************************************************************
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
fitFreeKnots <- function(y,x, 
                         weights = NULL, 
                           knots = NULL, 
                          degree = 3, 
                           fixed = NULL, 
                           trace = 0, 
                            data = NULL, 
                            base = c("trun","Bbase"), ...)
{
#-------------------------------------------------------------------------------
    penalty.opt <- function(kn, x, y, w, k, fixed = NULL, degree, ...) 
        {
       
       # kn <- sort(c(kn, fixed))
        if (length(u <- unique(knots)) < length(knots)) 
            stop(sprintf("%d coincident knot(s) detected", length(kn) - length(u)))
        #       fitFixedKnots(x, y, knots = knots, degree=degree, fixed = fixed)
        #cat("knots", kn, "\n")
        sp <-  fitFixedKnots(x=x, y=y, weights=w, knots = kn, degree=degree, base=base, ...)
        sp$rss
       }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
          ylab <- deparse(substitute(y))
          xlab <- deparse(substitute(x))
             y <- if (!is.null(data)) get(deparse(substitute(y)), envir=as.environment(data)) else y
             x <- if (!is.null(data)) get(deparse(substitute(x)), envir=as.environment(data)) else x
            w <- if(is.null(weights))  rep(1,length(y)) 
                 else 
                   {
                   if(!is.null(data))  get(deparse(substitute(weights)), envir=as.environment(data))
                   else weights
                   }
#if (is.data.frame(data)) { attach(data); on.exit(detach(data))}
      #  lx <- length(x)
       xmin <- min(x)
       xmax <- max(x)
    #     w <- if (is.null(w))  rep(1, lx)
    #        else rep(w, length = lx)
       eps <- 5e-04
     shift <- .Machine$double.eps^0.25
        # g <- length(knots)
    # initial fit
    #   sp0 <- fitFixedKnots(x=x, y=y, w=w, knots = knots, degree=degree, fixed = fixed, base=base)
    #sigma0 <-sp0$rss
    lambda <- if (length(knots) > 1) 
        {
            optim(knots, penalty.opt, #if (is.null(fixed)) penalty.gr, 
                  x = x, y = y, w = w, degree = degree, 
                lower = xmin + shift, upper = xmax - shift, method = "L-BFGS-B", 
                eps = eps, fixed = fixed, )$par #, control = list(trace = trace - 1))$par
        }
        else
        {
            optimize(penalty.opt, c(xmin, xmax), x = x, y = y, w=w,  degree=degree, fixed = fixed)$minimum
        }
          fit <- fitFixedKnots(x = x, y = y, weights = w, knots = lambda, degree = degree, fixed = fixed, base=base)
         out  <- list(call = sys.call(), fitted.values = fitted(fit), residuals = resid(fit), df = fit$df+length(lambda), 
                        rss = fit$rss, knots = fit$knots, coef = coef(fit), sigma = fit$sigma, 
                        fixed = fit$fixed, breakPoints = lambda,  degree = degree, base=fit$base,  y = y, x = x, X=fit$X, w = w, qr =fit$qr,  
                        deviance = sum(-2*dNO(y, mu=fitted(fit), sigma=fit$sigma)))
 class(out) <- c("FreeBreakPointsReg", "FixBreakPointsReg")
         out
 }
#-------------------------------------------------------------------------------
# print 
print.FreeBreakPointsReg<-function(x, digits=max(3, getOption("digits") - 3), ...) 
{
	cat("\nCall: ", deparse(x$call),  "\n", fill=TRUE)
 #cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n")
  cat("Estimated Knots: \n")
      print.default(format(knots(x), digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
    invisible(x)
}
#----------------------------------------------------------------------------------------
fitted.FreeBreakPointsReg<-function(object,...) 
{
object$fitted.values
}
#----------------------------------------------------------------------------------------
residuals.FreeBreakPointsReg<-function(object,...) 
{
object$residuals
}
#----------------------------------------------------------------------------------------
coef.FreeBreakPointsReg<-function(object,...) 
{
object$coef
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# SUMMARY
summary.FreeBreakPointsReg <- function(object, ...)
{
#--------------------------
#------------------
  #       x <- object$x
  #       y <- object$y 
  #       X <- object$X
  #  degree <- object$degree   
  #   fixed <- object$fixed
  #whatbase <- object$base   
  #  Nknots <- length(knots(object))
  #   Nbeta <- length(coef(object))
       par <- c(coef(object), knots(object), object$sigma)
  # Hessian <- hessian(par) 
        se <- vcov(object, type="se")
      tval <- par/se
   matcoef <- cbind(par, se, tval, 2*(1-pnorm(abs(tval))))
dimnames(matcoef) <- list(names(tval), c(" Estimate"," Std. Error", " t value", "Pr(>|t|)"))
cat("\nCall:\n", deparse(object$call), "\n\n", sep = "")
matcoef           
}
#----------------------------------------------------------------------------------------

vcov.FreeBreakPointsReg <- function (object, type=c("vcov", "cor", "se"), ...)
{
# Bbase base
 get.X.beta <- function(x, knots)  
   {
            # lx <- length(x)
             xl <- min(x)
             xr <- max(x)
           xmax <- xr + 0.01 * (xr - xl)
           xmin <- xl - 0.01 * (xr - xl)
         if (is.null(knots)) warning("No knots (break points) are declared")   
             kn <-  sort(c(xmin, xmax, fixed, knots)) 
       # splineDesign(knots, x, ord = 4, derivs, outer.ok = FALSE)
             X2 <- splineDesign(knots=kn, x, ord= degree + 1, derivs=0 * x, outer.ok=TRUE)
   colnames(X2) <- paste("XatBP",1:dim(X2)[2], sep="")
                   #bbase(x=x, xl=xmin, xr=xmax, knots=knots, deg=degree) # the problem here is the difference between break points
           form <- switch(as.character(degree), "0"=~1, "1"=~x, "2"=~x+I(x^2), "3"=~x+I(x^2)+I(x^3), 
                           "4"=~x+I(x^2)+I(x^3)+I(x^4), "5"=~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)) 
             X1 <- model.matrix(form) # model.matrix(~poly(x,degree)) 
              X <-cbind(X1,X2)
    X
    } 
#--------------------------
#-trunc base
 get.X.trun <- function(x, knots)  
   {
         tpower <- function(x, t, p)(x - t) ^ p * (x > t)      
           #  lx <- length(x)
           #  xl <- min(x)
           #  xr <- max(x)
           #xmax <- xr + 0.01 * (xr - xl)
           #xmin <- xl - 0.01 * (xr - xl)
         if (is.null(knots)) warning("No knots (break points) are declared")   
             kn <-  sort(c(fixed, knots)) 
       # splineDesign(knots, x, ord = 4, derivs, outer.ok = FALSE)
             X2 <- outer(x, kn, tpower, degree)# calculate the power in the knots
   colnames(X2) <- paste("XatBP",1:dim(X2)[2], sep="")
           form <- switch(as.character(degree), "0"=~1, "1"=~x, "2"=~x+I(x^2), "3"=~x+I(x^2)+I(x^3), 
                           "4"=~x+I(x^2)+I(x^3)+I(x^4), "5"=~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)) 
             X1 <- model.matrix(form) # model.matrix(~poly(x,degree)) 
              X <-cbind(X1,X2)
    X
    } 
#--------------------------
#------------------
# the log-likelihood function
LL <- function(par) 
   {   
       # this part depends on the parametrization
       # and also on the number of break points
        knots <- par[(Nbeta+1):(Nbeta+Nknots)]
            X <- if (whatbase=="Bbase") get.X.beta(x, knots) else get.X.trun(x, knots) #get.X(x, knots) # as a function of Break points 
         beta <- par[1:Nbeta] 
           mu <- X%*%beta
        sigma <- par[length(par)]
           mu <- X%*%beta
          llh <- -sum((dNO(x=y, mu=mu, sigma=sigma, log=TRUE))) # -logL
        llh 
    }
#------------------
#------------------
# the Hessian
hessian<-function(par)
  {
      npar <- length(par)
   epsilon <- 0.0001 * par
Hessian = matrix(0, ncol = npar, nrow = npar)
for (i in 1:npar) 
     {
     for (j in 1:npar) 
       {
                x1 <- x2 <- x3 <- x4 <- par
         x1[i] <- x1[i] + epsilon[i]; x1[j] <- x1[j] + epsilon[j] 
         x2[i] <- x2[i] + epsilon[i]; x2[j] <- x2[j] - epsilon[j]
         x3[i] <- x3[i] - epsilon[i]; x3[j] <- x3[j] + epsilon[j]
         x4[i] <- x4[i] - epsilon[i]; x4[j] <- x4[j] - epsilon[j]
 Hessian[i, j] <- (LL(x1)-LL(x2)-LL(x3)+LL(x4))/(4*epsilon[i]*epsilon[j])
       }
     }
  Hessian
   }
#------------------
      type <- match.arg(type)
         x <- object$x
         y <- object$y 
        # X <- object$X
    degree <- object$degree   
     fixed <- object$fixed
  whatbase <- object$base   
    Nknots <- length(knots(object))
     Nbeta <- length(coef(object))
       par <- c(coef(object), knots(object), object$sigma)
   Hessian <- hessian(par) 
       cov <- solve(Hessian)
  switch(type,"vcov" = cov,
               "cor" = cov2cor(cov), 
                "se" = sqrt(diag(cov)))
}
