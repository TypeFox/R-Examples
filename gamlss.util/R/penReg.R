#---------------------------------------------------------------------------------------
# Mikis Stasinopoulos 2-11-09
# revised 23-04-12
# TO DO : need predict()
#----------------------------------------------------------------------------------------
# the penalised regression function is ckecked against gamlss+pb()
# the results are identical or very similar 
# the only think I am concern is EM (ML-1)
#----------------------------------------------------------------------------------------
# this is a simple smoother using P-splines
# Paul Eilers and Mikis Stasinopoulos 
penReg <- function(y, x, 
            weights = rep(1,length(y)), 
                 df = NULL, 
             lambda = NULL,  
              start = 10,
              inter = 20, 
              order = 2,
             degree = 3,
               plot = FALSE,
             method = c("ML","ML-1","GAIC", "GCV", "EM"), 
                  k = 2,
                  ...)
{
#----------------------------------------------------
# inter : is the number of equal space intervals in x
# degree: is the degree of the polynomial
# order refers to differences in the penalty matrix
# order = 0 : white noise
# order = 1 : random effect
# order = 2 : random walk
# order = 3 : random walk of order 2
# lambda : the smoothing parameter
# df : the effective df's  
# if both lambda=NULL  and df=NULL then lambda is estimated using the Schall method
# if df is not NULL but lambda is NULL then df are used for smoothing
# if lambda is not NULL (whether df=NULL  or not) lambda is used for smoothing
# ---------------------------------------------------
# ---------------------------------------------------  
# creates the basis for p-splines
 bbase <- function(x, xl, xr, ndx, deg, quantiles=FALSE)
  {
 tpower <- function(x, t, p)
# Truncated p-th power function
    (x - t) ^ p * (x > t)
# DS xl= min, xr=max, ndx= number of points within 
# Construct B-spline basis
     dx <- (xr - xl) / ndx # DS increment 
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx) 
      P <- outer(x, knots, tpower, deg)# calculate the power in the knots
      n <- dim(P)[2]
      D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg) # 
      B <- (-1) ^ (deg + 1) * P %*% t(D) 
      B 
  }
# ---------------------------------------------------
# ---------------------------------------------------
# a siple penalized regression
regpen <- function(y, X, w, lambda, D)
  {
             G <- lambda * t(D) %*% D
            XW <- w * X
           XWX <- t(XW) %*% X
          beta <- solve(XWX + G, t(XW) %*% y)
            fv <- X %*%beta
             H <- solve(XWX + G, XWX)
         #  edf <- sum(diag(H))
           fit <- list(beta = beta, edf = sum(diag(H)))
  return(fit)  
  }
 # ---------------------------------------------------
 # ---------------------------------------------------
# a similar as obove but extra saving
regpenEM <- function(y, X, w, lambda, order, D)
  {
             G <- lambda * t(D) %*% D
            XW <- w * X
           XWX <- t(XW) %*% X
          beta <- solve(XWX + G, t(XW) %*% y)
            fv <- X %*%beta
             H <- solve(XWX + G, XWX)
             V <- solve(XWX + G)
           fit <- list(beta = beta, edf = sum(diag(H)), V=V)
  return(fit)  
  }
 # ---------------------------------------------------
 # ---------------------------------------------------
## function to find lambdas miimizing the local GAIC        
     fnGAIC <- function(lambda, k)
    {
    #cat("lambda, k", lambda, k, "\n")
       fit <- regpen(y=y, X=X, w=weights, lambda=lambda, D)
        fv <- X %*% fit$beta
       sig2 <- sum(weights*(y-fv)^2)/(length(y))
       NOd <- -2*dNO(y, mu=fv, sigma=sqrt(sig2), log=TRUE)    
      GAIC <- sum(weights*NOd)+k*fit$edf 
    #cat("GAIC", GAIC, "\n")
      GAIC   
    }
 # ---------------------------------------------------
 # ---------------------------------------------------
## function to find the lambdas wich minimise the local GCV 
      fnGCV <- function(lambda, k)
           {
    I.lambda.D <- (1+lambda*UDU$values)
           edf <- sum(1/I.lambda.D)
         y_Hy2 <- y.y-2*sum((yy^2)/I.lambda.D)+sum((yy^2)/((I.lambda.D)^2))
           GCV <- (n*y_Hy2)/(n-k*edf)
           GCV
           }  
 # ---------------------------------------------------
 # ---------------------------------------------------
## local function to get edf from lambda 
#   edf_df <- function(lambda)
#         {
#             G <- lambda * t(D) %*% D
#             H <- solve(XWX + G, XWX)
#           edf <- sum(diag(H))
#          # cat("edf", edf, "\n")
#           (edf-df)
#          }
## local function to get df using eigen values
    edf1_df <- function(lambda)
           {
           edf <-  sum(1/(1+lambda*UDU$values))
           (edf-df)
           }  
#------------------------------------------------------------------
# the main function starts here
#------------------------------------------------------------------
         scall <- deparse(sys.call())
        method <- match.arg(method)
            lx <- n <- length(x)
         inter <- if (lx<100) 10 else inter # this is to prevent singularities when length(x) is small
            xl <- min(x)
            xr <- max(x)
          xmax <- xr + 0.01 * (xr - xl)
          xmin <- xl - 0.01 * (xr - xl)   
             X <- bbase(x, xmin, xmax, inter, degree) # create the basis
             r <- ncol(X)
             D <- if(order==0) diag(r) else diff(diag(r), diff=order) # the penalty matrix
             if(!is.null(df)) # degrees of freedom
             {
             if (df>(dim(X)[2]-2)) 
              {df <- 3;  
              warning("The df's exceed
                      the number of columns of the design matrix", "\n",  "   they are set to 3") }
              df <- if (df < 1)  1  else  df+2
              if (df < 1)  warning("the df are set to 1")    
             }
#---------------------------------------------------
# case 1: if lambda is known just fit
 if (is.null(df)&&!is.null(lambda)||!is.null(df)&&!is.null(lambda))
 {
          fit <- regpen(y, X, weights, lambda,  D)
           fv <- X %*% fit$beta
         sig2 <- sum(weights*(y-fv)^2)/(length(y))     
 } # case 2: if lambda is estimated ------------------------------------------- 
 else if (is.null(df)&&is.null(lambda)) 
 { #   
  # cat("----------------------------","\n")
        lambda <- start# get(startLambdaName, envir=gamlss.env) ## geting the starting value
  # if ML ----------------------
  switch(method,
  "ML"={ 
       for (it in 1:50) 
         {
           fit  <- regpen(y, X, weights, lambda, D) # fit model
         gamma. <- D %*% as.vector(fit$beta)  # get the gamma differences
             fv <- X %*% fit$beta             # fitted values
           sig2 <- sum(weights * (y - fv) ^ 2) / (n - fit$edf)
           tau2 <- sum(gamma. ^ 2) / (fit$edf-order)# Monday, March 16, 2009 at 20:00 see LNP page 279
           if(tau2<1e-7) tau2 <- 1.0e-7 # MS 19-4-12
     lambda.old <- lambda
         lambda <- sig2 / tau2 # maybe only 1/tau2 will do since it gives exactly the EM results see LM-1
     if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
     if (lambda>1.0e+7) lambda<-1.0e+7 # DS Wednesday, July 29, 2009 at 19:00
      #    cat("iter tau2 sig2",it,tau2, sig2, '\n')
     if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e20) break
   #   assign(startLambdaName, lambda, envir=gamlss.env)
     #cat("lambda",lambda, '\n')
         }
        sig2 <- sum(weights*(y-fv)^2)/(length(y))  
       },
  "ML-1"={
       for (it in 1:50) 
         {
           fit  <- regpen(y, X, weights, lambda, D) # fit model
         gamma. <- D %*% as.vector(fit$beta)  # get the gamma differences
             fv <- X %*% fit$beta             # fitted values
         #  sig2 <- 1 # sum(w * (y - fv) ^ 2) / (n - fit$edf)
           tau2 <- sum(gamma. ^ 2) / (fit$edf-order)# Monday, March 16, 2009 at 20:00 see LNP page 279
           if(tau2<1e-7) tau2 <- 1.0e-7          
      lambda.old <- lambda
         lambda <- 1 / tau2 # 1/tau2 
     if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009
     if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e20) break
    #  assign(startLambdaName, lambda, envir=gamlss.env)
         }
       },
  "EM"={
      for (it in 1:500) 
         {
             fit  <- regpenEM(y, X, weights, lambda, order, D)
           gamma. <- D %*% as.vector(fit$beta)
           vgamma <- sum(diag(D%*%fit$V%*%t(D))) # this is crucial for estimating the variance of gamma Monday, March 23, 2009
               fv <- X %*% fit$beta
             tau2 <- ((sum(gamma.^ 2))+vgamma)/length(gamma.) 
             if(tau2<1e-7) tau2 <- 1.0e-7
       lambda.old <- lambda
           lambda <- 1 / tau2
         if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009    
       #    cat("iter sigma_t^2",it, tau2, "lambda",lambda, '\n')
       if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e20) break
         }
          sig2 <- sum(weights*(y-fv)^2)/(length(y))
    #cat("lambda",lambda, '\n')
     # assign(startLambdaName, lambda, envir=gamlss.env)
       },
  "GAIC"=
       {
        lambda <- nlminb(lambda, fnGAIC,  lower = 1.0e-7, upper = 1.0e20, k=k)$par 
           fit <- regpen(y=y, X=X, w=weights, lambda=lambda, D)
            fv <- X %*% fit$beta 
          sig2 <- sum(weights*(y-fv)^2)/(length(y))    
       },
  "GCV"={
  # 
           QR <-qr(sqrt(weights)*X)
           wy <- sqrt(weights)*y
          y.y <- sum(wy^2)
         Rinv <- solve(qr.R(QR))
            S <- t(D)%*%D
          UDU <- eigen(t(Rinv)%*%S%*%Rinv)
           yy <- t(UDU$vectors)%*%t(qr.Q(QR))%*%wy
       lambda <- nlminb(lambda, fnGCV,  lower = 1.0e-7, upper = 1.0e20, k=k)$par
          fit <- regpen(y=y, X=X, w=weights, lambda=lambda, D)
           fv <- X %*% fit$beta
         sig2 <- sum(weights*(y-fv)^2)/(length(y))     
       })
  }
  else # case 3 : if df are required---------------------------------
  { 
      #method 1
      #      XW <- w * X
      #     XWX <- t(XW) %*% X
      #  lambda <- if (sign(edf_df(0))==sign(edf_df(100000))) 100000  # in case they have the some sign
      #            else  uniroot(edf_df, c(0,100000))$root
      #method 2 from Simon Wood (2006) pages 210-211, and 360 
           QR <- qr(sqrt(weights)*X)
         Rinv <- solve(qr.R(QR))
          S   <- t(D)%*%D
          UDU <- eigen(t(Rinv)%*%S%*%Rinv)           
       lambda <- if (sign(edf1_df(0))==sign(edf1_df(100000))) 100000  # in case they have the some sign
                 else  uniroot(edf1_df, c(0,100000))$root
      # if (any(class(lambda)%in%"try-error")) {lambda<-100000}   
           fit <- regpen(y, X, weights, lambda, D)
            fv <- X %*% fit$beta
          sig2 <- sum(weights*(y-fv)^2)/(length(y))
  }#--------------------------------------------------------------------------end of case 3
#---------
         pfit <- list(coefficients = as.vector(fit$beta),
                            fitted = as.vector(fv), 
                                 y = y, 
                            ylabel =  substitute(y),
                                 x = x,
                           weights = weights,
                            xlabel =   substitute(x),
                            lambda = lambda, 
                            method = method,
                              call =  scall,
                                edf = fit$edf,
                             sigma = sqrt(sig2),
                              sig2 = sig2,
                                 N = n,
                               rss = sum(weights*(y-fv)^2),
                               aic = sum(weights*(-2*dNO(y, mu=fv, sigma=sqrt(sig2), log=TRUE)))+2*(fit$edf+1) , 
                               sbc = sum(weights*(-2*dNO(y, mu=fv, sigma=sqrt(sig2), log=TRUE)))+log(n)*(fit$edf+1),
                          deviance = sum(weights*(-2*dNO(y, mu=fv, sigma=sqrt(sig2), log=TRUE))))
  class(pfit) <- "penReg"
  if (plot==TRUE)
  {
  plot(y~x, col = gray(0.7))
   lines(fv[order(x)]~x[order(x)], col="red")
  }
  return(pfit)
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# methods
#----------------------------------------------------------------------------------------
fitted.penReg<-function(object,...) 
{
object$fitted
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
coef.penReg<-function(object,...) 
{
object$coefficients
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
residuals.penReg<-function(object,type = c("z-scores", "simple"), ...) 
{
  type <- match.arg(type)
  res <- if (type == "z-scores") 
  {qNO(object$weights*pNO(object$y, mu=object$fitted, sigma=object$sigma))}
  else
  {object$y-fitted(object)}
  res <- ifelse(res==-Inf, NA, res)
  res <- res[!is.na(res)]
  attr(res, "class") <- attr(object$y, "class")
  res
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#AIC.penReg <- function(object, ...,k=2)
#{
# val <- if (is(object, "penReg")) 
#            object$deviance + (object$edf+1) * k
#        else stop(paste("this is not a penReg object"))
#val
#}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
AIC.penReg <- function(object, ...,k=2)
{
  if (length(list(...))) {
    object <- list(object, ...)
    df <- 1+as.numeric(lapply(object, function(x) x$edf))
    AIC <- as.numeric(lapply(object, function(x) deviance(x) + 
      (x$edf+1) * k))
    val <- as.data.frame(cbind(df, AIC))
    Call <- match.call()
    Call$k <- NULL
    row.names(val) <- as.character(Call[-1])
    val <- val[order(AIC), ]
    val
  }
  else 
  {
    val <- if (is(object, "penReg")) 
      object$deviance + (object$edf+1) * k
    else stop(paste("this is not a penReg object"))
  }
  val
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
deviance.penReg<-function(object,...) 
{
object$deviance
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
plot.penReg <- function (x, xvar=NULL, parameters=NULL, ts=FALSE, summaries=TRUE, ...) 
{
  residx <- resid(x) # get the residuals 
  w <- x$weights
  x$x <- if (is(x$y , "ts"))  as.numeric(time(x$y))
  else                 1:length(x$y)
  xlabel <- if(!missing(xvar)) deparse(substitute(xvar)) else deparse(substitute(index))
  ## plotting parameters
  if(is.null(parameters))
    op <- par(mfrow=c(2,2), mar=par("mar")+c(0,1,0,0), col.axis="blue4", col.main="blue4", col.lab="blue4",  col="darkgreen", bg="beige" )
  else  op <- parameters
  ## now the two top  figures 
  ## if time series plot acf and pacf  
  if(identical(ts, TRUE))
  {  # get the acf and pacf
    acf.new<-acf(residx,plot=FALSE)
    plot(acf.new,xlim=c(2,length(acf.new$acf)),ylim=range(acf.new$acf[-1]))   # ms Tuesday, August 19, 2003 at 11:04
    pacf(residx)
  }
  else 
  {# otherwise 
    ## I am assuming that is x$noObs!=x$N then we have weights (with frequencies)
    if (length(residx)==x$N)
    {
      fittedvalues <- if(is.null(fitted(x))) fitted(x,"sigma") else fitted(x) # MS Wednesday, September 10, 2003 at 21:20
      ## whether index or x-variable
      if(is.null(xvar))     xvar <- seq(1,length(residx),1) # MS
    }
    else
    { # if weights
      fittedvalues <- rep( if(is.null(fitted(x))) fitted(x,"sigma") else fitted(x), w)
      xvar <- if(is.null(xvar))  seq(1,length(residx),1) else rep(xvar,w)
    } 
    # top left
    plot(fittedvalues , residx,
         xlab = "Fitted Values",  
         ylab = "Quantile Residuals", 
         main = "Against Fitted Values",
         frame.plot = TRUE) 
    # top right  
    plot(xvar, residx, 
         ylab = "Quantile Residuals",
         xlab = xlabel, 
         main = paste("Against ", xlabel), 
         frame.plot = TRUE) #  points(par(col="blue4"))
  }    
  plot(density(residx), 
       xlab = "Quantile. Residuals", 
       ylab = "Density", 
       main = "Density Estimate",
       frame.plot = TRUE, 
       col="black", 
       lwd=0.4 ) #col="deepskyblue4", col="darkgreen", 
  rug(residx, col="red")
  
  qqnorm(residx, main = "Normal Q-Q Plot",
         xlab = "Theoretical Quantiles",
         ylab = "Sample Quantiles", 
         plot.it = TRUE, 
         frame.plot = TRUE, 
         col="darkgreen")
  lines(as.numeric(residx), as.numeric(residx), col="red" , lwd=.4, cex=.4 )
  
  if ( identical(summaries, TRUE))
  { 
    qq <- as.data.frame(qqnorm(residx, plot = FALSE))
    Filliben <- cor(qq$y,qq$x)
    # mr <- as.matrix(residx)
    m.1 <- mean(residx)
    m.2 <- var(residx) # cov.wt(mr,w)$cov
    n.obs <- sum(w) 
    m.3 <- sum((residx-m.1)**3)/n.obs 
    m.4 <- sum((residx-m.1)**4)/n.obs 
    b.1 <- m.3^2/m.2^3
    sqrtb.1 <- sign(m.3)*sqrt(abs(b.1))
    b.2 <- m.4/m.2^2 
    cat("*******************************************************************")
    cat("\n")
    if (identical(x$type,"Continuous")) 
    {cat("\t","     Summary of the Quantile Residuals")}
    else{cat("\t","Summary of the Randomised Quantile Residuals")}    
    cat("\n")
    cat("                           mean   = ", m.1, "\n")
    cat("                       variance   = ", m.2, "\n")
    cat("               coef. of skewness  = ", sqrtb.1, "\n")
    cat("               coef. of kurtosis  = ", b.2, "\n")
    cat("Filliben correlation coefficient  = ", Filliben, "\n")
    cat("*******************************************************************")
    cat("\n")
    
  }
  par(op)
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
print.penReg  <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{   
  cat("\nPenalised Regression fit")
  cat("Fitting method:", deparse(x$method), "\n")
  cat("\nCall: ", deparse(x$call), "\n", fill = TRUE)
 # cat("Coefficients")
  #co <- coef(x)
  #cat("\n  ", names(co), "\n")
  #cc <-simplify2array(co, higher=TRUE)
  #cat(cc, " \n")
  cat("\n Degrees of Freedom for the fit:", x$df, "Residual Deg. of Freedom  ", 
      x$N-x$df, "\n")
  cat("Global Deviance:    ", format(signif(x$deviance)), 
      "\n            AIC:    ", format(signif(x$aic)), "\n            SBC:    ", 
      format(signif(x$sbc)), "\n")
  invisible(x)
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
summary.penReg  <- function (object, digits = max(3, getOption("digits") - 3), ...) 
{   
  cat("\nDiscrete Smoothing fit")
  cat("Fitting method:", deparse(object$method), "\n")
  cat("\nCall: ", deparse(object$call), "\n", fill = TRUE)
  
  
  cat("\n Degrees of Freedom for the fit:", object$df, "Residual Deg. of Freedom  ", 
      object$N-object$df, "\n")
  cat("Global Deviance:    ", format(signif(object$deviance)), 
      "\n            AIC:    ", format(signif(object$aic)), "\n            SBC:    ", 
      format(signif(object$sbc)), "\n")
  invisible(object)
}

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------