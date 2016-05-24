#---------------------------------------------------------------------------------------
# Mikis Stasinopoulos 26-11-11
# functions 
#  i) penLags()
# ii) checkStationarity()
# iii) METHODS FOR penlags
#     1) fitted()
#     2) coef()
#     3) resid()
#     4) AIC()
#     5) deviance()
#     6) print()
#     7) summary()
#     8) plot()
#     9) predict()
#
# in this version the constant is inclueded in the model
# a small constant was added for the constant after suggestion from Paul
# there is no GCV anymore
# and also if you fix the df the resulting model has edf similat but not identical
# this is due t the way the df's are calculated using S. Wood method which do not all
# We need to add a delay option in case that there is delay in reaction
# TO DO
#   i)  standard  errors for the beta's' at the moment it produces rubbish
#  ii) se for the fitted values 
# iii)  predict   OK
#  iv)  forecast  
#   v) the Q function implementation 
#  we have a problem that we can not pass no.x as argument in case of x not identical to y
# the above is fixed
#-------------------------
# TO DO 
# 1) blag should be able to start from a different point rathen that 0 or 1 This involves of changing blag with a new argument say from.lag: OK done
# the above makes the no.x redundand
# 2) check for AR modeling in la
#-------------------------------------------------------------------------------
# the penalised lags regression
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
penLags <- function(y, x, 
               lags = 10,
           from.lag = 0,
            weights = NULL, 
               data = NULL,    
                 df = NULL, 
             lambda = NULL,  
       start.lambda = 10,
              order = 1,
               plot = FALSE,
             method = c("ML", "GAIC"), 
                  k = 2,
                  ...)
{
#----------------------------------------------------
# local functions
# order refers to differences in the penalty matrix
# order = 0 : white noise
# order = 1 : random effect
# order = 2 : random walk
# order = 3 : random walk of order 2
# lambda : the smoothing parameter
# df : the effective df's  
# if both lambda=NULL  and df=NULL then lambda is estimated using local ML (PQL)
# if df is not NULL but lambda is NULL then df are used for smoothing
# if lambda is not NULL (whether df=NULL  or not) lambda is used for smoothing
# ------------------------------------------------------------------------------
# a simple penalized regression where the first column is not penalised  
regpen <- function(y, X, w, lambda, D)
  {
  	         DD <- t(D) %*% D # MS 26-11-11
            # G <- lambda * t(D) %*% D
             rn <- dim(X)[2] # this to add zeros for the constant
         topRow <- rep(0,(rn-1)) # first row
      #      rG <- rbind(topRow, G)
            rDD <- rbind(topRow, DD)    
       firstCol <- rep(0,rn)    # first col
    firstCol[1] <- 0.1e-10 #  this is sugessted by Paul rather than zero
            #NG <- cbind(firstCol,rG)      
             NG <- cbind(firstCol,rDD)
             NG <- lambda * NG      
             XW <- w * X
            XWX <- t(XW) %*% X
             P0 <- 1e-3 * diag( dim(NG)[1]) # new MIkis
           beta <- solve(XWX + NG+P0, t(XW) %*% y)
            #fv <- X %*%beta
              H <- solve(XWX + NG+P0, XWX)
         #  edf <- sum(diag(H))
           fit <- list(beta = beta, edf = sum(diag(H)), trace=diag(H))
  return(fit)  
  }
#--------------------------------------------------#---------------------------
## function to find lambdas miimizing the local GAIC        
     fnGAIC <- function(lambda, k)
    {
    #cat("lambda, k", lambda, k, "\n")
         fit <- regpen(y=y, X=X, w=weights*pw, lambda=lambda, D)
          fv <- X %*% fit$beta
        sig2 <- sum(weights*(y-fv)^2)/(n - fit$edf)
         NOd <- -2*dNO(y, mu=fv, sigma=sqrt(sig2), log=TRUE)    
        GAIC <- sum(weights*NOd)+k*fit$edf 
    #cat("GAIC", GAIC, "\n")
      GAIC   
    }
#-------------------------------------------------------------------------------
## function to find the lambdas wich minimise the local GCV 
#      fnGCV <- function(lambda, k)
#           {
#    I.lambda.D <- (1+lambda*UDU$values)
#           edf <- sum(1/I.lambda.D)
#         y_Hy2 <- y.y-2*sum((yy^2)/I.lambda.D)+sum((yy^2)/((I.lambda.D)^2))
#           GCV <- (n*y_Hy2)/(n-k*edf)
#           GCV
#           }  
#-------------------------------------------------------------------------------
    edf1_df <- function(lambda)
           {
           edf <-  sum(1/(1+lambda*UDU$values))
           (edf-df)
           }  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# the main function starts here
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# we only need y  x and w
if (!is.null(data)) 
 {
  y <- with(data, y)
  x <- with(data, x)
 }
weights <- if(is.null(weights))  rep(1,length(y)) 
          else {
                if(!is.null(data))  with(data, weights)
                else weights
                }
 #          if (is.data.frame(data)) { attach(data); on.exit(detach(data))} 
           if (lags <= order+1) stop("Increase the number of lags") 
          scall <- deparse(sys.call())
         method <- match.arg(method)
      Extra.Arg <- list(...)
             lx <- n <-  length(x)  
           lags <- if (lx<99) 10 else lags # this is to prevent singularities when length(x) is small
            # I need better statement here 
       from.lag <- if (identical(y,x)&&from.lag==0)  1 else from.lag
#browser()
            XX  <-  blag(x=x, lags=lags, value=x[1], from.lag=from.lag) # note that we condition onthe  first ellement of x[1], This do not effect the fit but effects the 
#                  else  # fitted values 
#                        { 
#                        if (is.null(Extra.Arg$x.no))
#                          blag(x=x, lags=lags, value=x[1], no.x=FALSE)
#                        else blag(x=x, lags=lags, value=x[1], no.x=Extra.Arg$x.no)
#                        }
            pw <- wlag(XX) # get the lag weights
             r <- ncol(XX)
             X <- cbind(rep(1, lx), XX)  #  put one in the fist colunm  
   colnames(X) <- unlist(c("1", attr(XX,"dimnames")[2]))
    attr(X,"no.lags") <- attr(XX,"no.lags")
       D <- if(order==0) diag(r) else diff(diag(r), diff=order) # the penalty matrix
             if(!is.null(df)) # degrees of freedom
             {
             if (df>(dim(X)[2])) 
              {df <- 3;  
              warning("The df's exceed the number of columns of the design matrix", "\n",  "   they are set to 3") }
              df <- if (df < 1)  1  else  df
              if (df < 1)  warning("the df are set to 1")    
             }
          par <- rep(0,2)
   allweights <- weights*pw  
#---------------------------------------------------
# case 1: if lambda is known just fit
 if (is.null(df)&&!is.null(lambda)||!is.null(df)&&!is.null(lambda))
 {
          fit <- regpen(y, X, allweights, lambda,  D)
           fv <- X %*% fit$beta
         sig2 <- sum(allweights*(y-fv)^2)/n     
 } # case 2: if lambda is estimated ------------------------------------------- 
 else if (is.null(df)&&is.null(lambda)) 
 { #   
  # cat("----------------------------","\n")
        lambda <- start.lambda# get(startLambdaName, envir=gamlss.env) ## geting the starting value
  # if ML ----------------------
  switch(method,
  "ML"={ 
       for (it in 1:50) 
         {
           fit  <- regpen(y, X, allweights, lambda, D) # fit model
         gamma. <- D %*% as.vector(fit$beta[-1])  # get the gamma differences
             fv <- X %*% fit$beta             # fitted values
           sig2 <- sum(allweights*(y - fv) ^ 2) / (n - fit$edf) # or (n - fit$edf)
           tau2 <- sum(gamma. ^ 2) / (fit$edf-order)# Monday, March 16, 2009 at 20:00 see LNP page 279
     lambda.old <- lambda
         lambda <- sig2 / tau2 # maybe only 1/tau2 will do since it gives exactly the EM results see LM-1
     if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
     if (lambda>1.0e+7) lambda<-1.0e+7 # DS Wednesday, July 29, 2009 at 19:00
   #      cat("iter tau2 sig2",it,tau2, sig2, '\n')
     if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e20) break
   #   assign(startLambdaName, lambda, envir=gamlss.env)
   #  cat("lambda",lambda, '\n')
          }
        sig2 <- sum(weights*(y-fv)^2)/(n - fit$edf) # see above whether n-fir$edf or even n 
        par <- c(sig2, tau2)
   names(par) <- c("sig^2", "tau^2")
         },
#  "MLQ"={ 
#           fit  <- regpen(y, X, weights*pw, lambda, D)
#           fv <- X%*%fit$beta
#       params <- c(sige = (sum(weights*(y-fv)^2)/(sum(weights)-fit$edf)), sigb=sum(fit$beta^2)/(fit$edf-order))        
#            G <- t(D) %*% D
#            W <- diag(weights*pw)
#        Qf<- function(par) 
#     { # maybe we reparametrice to work on the log scale of the parameters 
#      # cat("par", par, "\n")
#       
#       lambda <- par[1]/par[2]
#           fv <- solve(W+lambda*G,  weights*pw*y)   
#           D3 <- determinant((1/par[1])*W+(1/par[2])*G)$modulus
#            f <- -(N/2)*log(2*pi*par[1]) -sum(weights*(y-fv)^2)/(2*par[1])-((N-order)/2)*log(2*pi*par[2])-sum(b^2)/(2*par[2]) -.5*D3+(N/2)*log(2*pi)
##            b <- diff(fv, differences=order)
#           -f
#     }
 #switch(optim.p, "nlminb"={# this do not allow se's'
#            out <- nlminb(start = params, objective = Qf, lower = c(1e-10, 1e-10),  upper = c(Inf, Inf))
   # out$hessian <- HessianPB(pars=out$par, fun=Qf)$Hessian
#     value.of.Q <- out$objective
  #                       },
  #               "optim"={
  #          out<- optim(par= params, fn = Qf, lower = c(1e-10, 1e-10), upper = c(Inf, Inf),method= "L-BFGS-B",
  #                 hessian=TRUE)
                  # geting the hessian using nlme    fdHess(out$par, Qf)
  #          value.of.Q <- out$value
  #                         },
   #               "genoud"={
   #                 require(rgenoud)
   #                 B <- matrix(c(1e-10,1e-10, 1e10, 1e10), ncol=2)
    #                # spead is INCREASED  by reducing  pop.size and increasing BFGSburnin
    #          out <- genoud( Qf, nvar=2, max=FALSE, pop.size=20,  Domains=B, boundary.enforcement=2, gradient.check=TRUE, hessian=TRUE,BFGSburnin=20, max.generations=50,print.level=0)
    #        # out <- genoud( Qf, nvar=2, max=FALSE, pop.size=10,  Domains=B, boundary.enforcement=2, hessian=TRUE, BFGSburnin=50, P1=0, P2=50, P3=50, P4=0, P5=0, P6=0, P7=0, P8=0, P9=0)
    #        #out <- genoud( Qf, nvar=2, max=FALSE, pop.size=5000, Domains=B, boundary.enforcement=2, hessian=TRUE)
     #               value.of.Q <- out$value   
     #               })
#        par <- c(sig2, tau2)
#   names(par) <- c("sig^2", "tau^2")
#         },
  "GAIC"={
        lambda <- nlminb(lambda, fnGAIC,  lower = 1.0e-7, upper = 1.0e20, k=k)$par 
           fit <- regpen(y=y, X=X, w=weights*pw, lambda=lambda, D)
            fv <- X %*% fit$beta 
         sig2 <- sum(weights*(y-fv)^2)/(n - fit$edf)   
        #assign(startLambdaName, lambda, envir=gamlss.env)
        } 
#        ,
#  "GCV"={
#           QR <-qr(sqrt(weights)*X)
#           wy <- sqrt(weights)*y
#          y.y <- sum(wy^2)
#         Rinv <- solve(qr.R(QR))
#            S <- t(D)%*%D
#          UDU <- eigen(t(Rinv)%*%S%*%Rinv)
#        browser()
#           yy <- t(UDU$vectors)%*%t(qr.Q(QR))%*%wy
#       lambda <- nlminb(lambda, fnGCV,  lower = 1.0e-7, upper = 1.0e20, k=k)$par
#          fit <- regpen(y=y, X=X, w=weights*pw, lambda=lambda, D)
#           fv <- X %*% fit$beta
#        sig2 <- sum(weights*(y-fv)^2)/(sum(weights) - fit$edf)    
      #  assign(startLambdaName, lambda, envir=gamlss.env) 
#       }
  )
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
           G  <- t(D)%*%D
           rn <- dim(X)[2] # this to add zeros for the constant
       topRow <- rep(0,(rn-1)) # first row
           rG <- rbind(topRow, G) 
     firstCol <- rep(0,rn)    # first col
            S <- cbind(firstCol,rG)
          UDU <- eigen(t(Rinv)%*%S%*%Rinv)           
       lambda <- if (sign(edf1_df(0))==sign(edf1_df(100000))) 100000  # in case they have the same sign
                 else  uniroot(edf1_df, c(0,100000))$root
      # if (any(class(lambda)%in%"try-error")) {lambda<-100000}   
           fit <- regpen(y, X, weights*pw, lambda, D)
            fv <- X %*% fit$beta
         sig2 <- sum(weights*(y-fv)^2)/(n - fit$edf)
  }#--------------------------------------------------------------------------end of case 3
 # the variance of the fitted values?
 #  browser()
   #   waug <- as.vector(c(weights, rep(1,nrow(D))))
   #   xaug <- as.matrix(rbind(X,sqrt(lambda)*D))
   #    lev <- hat(sqrt(waug)*xaug,intercept=FALSE)[1:n] # get the hat matrix
     # DO I NEED THIS? lev <- (lev-.hat.WX(weights,x)) # subtract  the linear since is already fitted 
   #    var <- lev/weights              # the variance of the smoother
#---------
#cat("sig2", sig2, "\n")
                               fv <- as.vector(fv)
     if (is.ts(y)) attributes(fv) <- attributes(y)  
         pfit <- list(coefficients = as.vector(fit$beta),
                          var.coef = as.vector(fit$trace),
                     fitted.values = fv,
                  #             var = var,
                              call = scall, 
                            method = method,
                                 y = y, 
                           weights = allweights,
                            ylabel =  substitute(y),
                                 x = x,
                                 X = X,  
                                 N = length(y)-lags,
                            xlabel =   substitute(x),
                            lambda = lambda, 
                             mu.df = fit$edf,
                            df.fit = fit$edf+1,
                             sigma = sqrt(sig2),
                               par = par,
                               rss = sum(allweights*(y-fv)^2),
                               aic = sum(allweights*(-2*dNO(y, mu=fv, sigma=sqrt(sig2), log=TRUE)))+2*(fit$edf+1) , 
                               sbc = sum(allweights*(-2*dNO(y, mu=fv, sigma=sqrt(sig2), log=TRUE)))+log(n)*(fit$edf+1),
                          deviance = sum(allweights*(-2*dNO(y, mu=fv, sigma=sqrt(sig2), log=TRUE))))
  class(pfit) <- "penlags"
  if (plot==TRUE)
  {
  plot(y, col = gray(0.7), type = 'l')
   lines(fv~as.numeric(time(y)), col = 'red', lwd = 2)
  }
  return(pfit)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
fitted.penlags<-function(object,...) 
{
object$fitted.values
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
coef.penlags<-function(object, what=c("All", "AR", "varComp"), ...) 
{
 what <- match.arg(what) 
  coef <- switch(what,"All"= object$coefficients,
                       "AR"= object$coefficients[-1],
                  "varComp"= object$par)
  coef
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
residuals.penlags<-function(object,...) 
{
#object$y-object$fitted
res<-qNO(object$weights*pNO(object$y, mu=object$fitted, sigma=object$sigma))
res <- ifelse(res==-Inf, NA, res)
res <- ifelse(res==Inf, NA, res)
res[!is.na(res)]  
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
AIC.penlags <- function(object, ...,k=2)
{
    if (length(list(...))) {
        object <- list(object, ...)
        ispenlags <- unlist(lapply(object, is, "penlags"))
        if (!any(ispenlags)) 
            stop("some of the objects are not penlags")
        df <- as.numeric(lapply(object, function(x) x$df))
        AIC <- as.numeric(lapply(object, function(x) x$deviance + 
            x$df * k))
        val <- as.data.frame(cbind(df, AIC))
        Call <- match.call()
        Call$k <- NULL
        row.names(val) <- as.character(Call[-1])
        val <- val[order(AIC), ]
        val
    }
    else {
        val <- if (is(object, "penlags")) 
            object$deviance + object$df * k
        else stop(paste("this is not a gamlss object"))
        val
    }
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
deviance.penlags<-function(object,...) 
{
object$deviance
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

print.penlags  <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{   
    cat("\nPenalised Lags fit")
    cat("Fitting method:", deparse(x$method), "\n")
    cat("\nCall: ", deparse(x$call), "\n", fill = TRUE)
    cat("Coefficients")
    co <- coef(x)[1:5]
    cat("\n  ", names(co),"\n")
    cc <-simplify2array(co, higher=TRUE)
    cat(cc, "... \n")
    cat("\n Degrees of Freedom for the fit:", x$df, "Residual Deg. of Freedom  ", 
        x$N-x$df, "\n")
    cat("Global Deviance:    ", format(signif(x$deviance)), 
        "\n            AIC:    ", format(signif(x$aic)), "\n            SBC:    ", 
        format(signif(x$sbc)), "\n")
    invisible(x)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
summary.penlags  <- function (object, digits = max(3, getOption("digits") - 3), ...) 
{   
    cat("\nDiscrete Smoothing fit")
    cat("Fitting method:", deparse(object$method), "\n")
    cat("\nCall: ", deparse(object$call), "\n", fill = TRUE)
            coef  <- coef(object)
          se.coef <- sqrt(object$var.coef)
             tval <- coef/se.coef
          matcoef <- cbind(coef, se.coef, tval, 2*(1-pnorm(abs(tval))))
dimnames(matcoef) <- list(names(tval), c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
cat("\nCoefficient(s):\n")
printCoefmat(matcoef, digits = 6, signif.stars = TRUE)

   
    cat("\n Degrees of Freedom for the fit:", object$df, "Residual Deg. of Freedom  ", 
        object$N-object$df, "\n")
    cat("Global Deviance:    ", format(signif(object$deviance)), 
        "\n            AIC:    ", format(signif(object$aic)), "\n            SBC:    ", 
        format(signif(object$sbc)), "\n")
    invisible(object)
}
# methods are finish here    
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#plot.penlags <- function(x, ...)
#{
#  residx <- (x$y-fitted(x))/x$sigma
# op <- par(mfrow=c(2,2), mar=par("mar")+c(0,1,0,0), col.axis="blue4", col.main="blue4", col.lab="blue4",  col="darkgreen", bg="beige" ) 
# plot(fitted(x) , residx,
#         xlab = "Fitted Values",  
#         ylab = "Residuals", 
#         main = "Against Fitted Values",
#         frame.plot = TRUE) 
#    # top right  
#    plot(x$x, residx, 
#         ylab = "Residuals",
#         xlab = x$xlabel, 
#         main = paste("Against ", x$xlabel), 
#         frame.plot = TRUE) #  points(par(col="blue4"))  
#  plot(density(resid(x)), 
#         xlab = "Residuals", 
#         ylab = "Density", 
#         main = "Density Estimate",
#         frame.plot = TRUE, 
#         col="black", 
#         lwd=0.4 ) #col="deepskyblue4", col="darkgreen", 
#   qqnorm(residx, main = "Normal Q-Q Plot",
#          rug(residx, col="red")
#            xlab = "Theoretical Quantiles",
#            ylab = "Sample Quantiles", 
#            plot.it = TRUE, 
#            frame.plot = TRUE, 
#            col="darkgreen")
#     lines(as.numeric(residx), as.numeric(residx), col="red" , lwd=.4, cex=.4 )
#par(op)
#}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
plot.penlags <- function (x, xvar=NULL, parameters=NULL, ts=FALSE, summaries=TRUE, ...) 
{
    residx <- resid(x) # get the residuals 
         w <- x$weights
       x$x <- if (is(x$y , "ts"))  as.numeric(time(x$y))
              else                 1:length(x$y)
       x$x <- x$x[w>0]       
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
     #browser()
     if (length(residx)==x$N)
        {
         fittedvalues <- fitted(x)[w>0] # MS Wednesday, September 10, 2003 at 21:20
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
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
predict.penlags <- function(object, how.many=NULL,...)
{
   if (is.null(how.many))
   {
   return(fitted(object))
   }
   else
   {

  no.lags <-  attr(object$X, "no.lags")
     pred <- rep(0, how.many)
   #   lenX <- dim(object$X)[1]
   last.X <- tail(object$X, 1)
   for (i in 1:how.many)
   {                     
      if (i <= no.lags)
       {
           nval <- sum(last.X*coef(object))
      last.X[i+1] <- nval
        pred[i] <- nval
       }
    else
       {
       pred[i] <- nval 
       }
   }
  return(pred) 
   } 
  }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
checkStationarity <- function(par)
{
minroots <- min(Mod(polyroot(c(1, - par))))
        if (minroots <= 1) 
            cat("'ar' part of model is not stationary")
        else cat("OK")
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------     
#polyroot(c(1, 2, 1, 0, 0))