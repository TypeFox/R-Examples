# Fernanda De Bastiani    
# there are two main functions here 
#  i) MRFA
# ii) MRF   
# those functions allow to fit a simple normal mrf model for a response variable where the only explanatory factor is the area/district
# this could be then used withing GAMLSS to fitted mrf models
# ------------------------------------------------------------------------------
# Function MRFA
# To do  
# i)  for the alternate method MRFA() we may include fixing the df's that will be useful for checking 
# ii) There is a problem with fixing lambda: for example if we fix lambda to be 10 we should have sig2e/sig2b=10 but this is not the  case.  It looks that this is because the heequation lambda=sig2e/sig2b is not true  when lambda is not on its maximum ()
#  iii) check levels(X) and levels(k) that is whether there are more levels  AREA than the explanatory factor        OK fixed
#      new more examples
# # vi) check the weights OK
#-------------------------------------------------------------------------------
# Function MRF()
#  i) check levels(X) and levels(k) that is whether there are more levels in AREA than the explanatory factor OK
#     new more examples
##------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#library(spam) # this should be taken out after testing MS
#library(mgcv)
#library(gamlss)


#-------------------------------------------------------------------------------
# the Q-function methods
#-------------------------------------------------------------------------------
MRF <-     function(y, x,  
                    precision = NULL,
                    neighbour = NULL,
                    polys = NULL,
                    area = NULL, 
                    weights = rep(1,length(y)), # for weighted observations 
                    sig2e = 1,
                    sig2b = 1,
                    #plot = FALSE,
                    sig2e.fix = FALSE,
                    sig2b.fix = FALSE,
                    penalty = FALSE,
                    delta = c(0.01, 0.01),
                    shift = c(0,0))                     
 {
  
  ###########################################################
  ## local 2 : creating the precision matrix from the neighbour information
  ##          same as the external marix with the same name but here takes 
  ##          local information
  nb2prec <- function(neighbour)
  {
    a.name <- names(neighbour$nb)
    if (all.equal(sort(a.name), sort(levels(k))) != TRUE) 
      stop("mismatch between neighbour/polys supplied area names and data area names")
    np <- nfv
    G <- matrix(0, np, np)
    rownames(G) <- colnames(G) <- levels(k)
    for (i in 1:np) 
    {
      ind <- neighbour$nb[[i]]
      lind <- length(ind)
      G[a.name[i], a.name[i]] <- lind
      if (lind > 0) 
        for (j in 1:lind) G[a.name[i], a.name[ind[j]]] <- -1
    }
    if (sum(G != t(G)) > 0) 
      stop("Something wrong with auto- penalty construction")
    if (all.equal(sort(a.name),a.name)!=TRUE) 
    { ## re-order penalty to match X
      G <- G[levels(x),] 
      G <- G[,levels(x)]
    }
    G
  }
  
#------------------------------------------------------------------------------
###########################################################
#----------------------------------------------
# The Q function = Marginal likelihood
    Qf <- function(par)
      {  
        #    cat("par", exp(par), "\n")
          beta <- solve(exp(-par[1])*XWX + exp(-par[2])*G, (exp(-par[1])*XWy)) 
            fv <- beta[x]
            DT <- determinant((exp(-par[1])*XWX) + (exp(-par[2])*G))$modulus
             f <- -(nobs/2)*log(2*pi*exp(par[1])) -           # 1 
                   .5*sum(log(ifelse(weights==0,0.0001,weights)))-                  # 2 this do not depend on parameters
                   sum(weights*(y-fv)^2)/(2*exp(par[1])) + # 3 
                   (nfv/2)*log(2*pi) -                    # 4 ?
                   ((nfv)/2)*log(2*pi*exp(par[2]))-       # 5 
                   sum(t(beta)%*%G%*%beta)/(2*exp(par[2])) -           # 6 
                   .5*DT								   # 7	
attributes(f) <- NULL
#cat("f", f, "\n")
      Pen <- if (penalty==TRUE) delta[1]*(par[1]-shift[1])^2 + delta[2]*(par[2]-shift[2])^2 else 0       
      -f+Pen        # no penalty  yet    
      }
#-------------------------------------------------
########################################################### 
#  require(spam); require(gamlss)
  scall <- deparse(sys.call(), width.cutoff = 500L)
  if (!is(x, "factor")) stop("x must be a factor")
  #x <- as.factor(data[[object$term]])
         N <- length(y)
  if (any(is.na(weights))) weights[is.na(weights)] <- 0
      nobs <- sum(!weights==0)   
         k <- area
  if (is.null(k))
         k <- factor(levels(x),levels=levels(x)) # default knots = all regions in the data
  else   k <- as.factor(k)
  if (length(levels(x))==length(levels(k))) 
  if (length(levels(x))>length(levels(k))) stop("MRF basis dimension set too high")
  if (sum(!levels(x)%in%levels(k)))
    stop("data contain regions that are not contained in the area specification") 
          x <- factor(x,levels=levels(k))
          X <- model.matrix(~x-1,)        
      nfv <- nlevels(x)
  if (is.null(precision)&&is.null(neighbour)&&is.null(polys))
    stop("precision matrix, boundary polygons and/or neighbours list must be supplied")
  if (!is.null(precision))
  { 
    if (!is.matrix(precision)||dim(precision)[1]!=nfv||dim(precision)[2]!=nfv) 
      stop("the precision matrix is not suitable")
        G <- as.spam(precision) 
  } 
  # check the precision matrix
  if (!is.null(neighbour)&&is.null(precision))
  { # if neighbour exits then calculate the precision
     G   <- as.spam(nb2prec(neighbour))
  }
  if (!is.null(polys)&&is.null(neighbour)&&is.null(precision))  
  { # if polys exits then calculate the precision
    a.name <- names(polys)
    d.name <- unique(a.name[duplicated(a.name)])
    if (length(d.name)) 
    {
      for (i in 1:length(d.name)) 
      {
        ind <- (1:length(a.name))[a.name == d.name[i]]
        for (j in 2:length(ind)) polys[[ind[1]]] <- rbind(polys[[ind[1]]], 
                                                          c(NA, NA), polys[[ind[j]]])
      }
      #now delete the un-wanted duplicates
      ind <- (1:length(a.name))[duplicated(a.name)]
      if (length(ind) > 0) 
        for (i in length(ind):1) polys[[ind[i]]] <- NULL
    }#polygon list in correct format
    neighbour <- polys2nb(polys)
    G   <- as.spam(nb2prec(neighbour))
  }
         XWX <- diag.spam(diag(t(X)%*%diag(weights)%*%X))
    #    XWX <- diag.spam(tapply(weights, x, sum)) 
         XWy <- as.vector(t(X)%*%diag(weights)%*%matrix(y, nrow=length(y), ncol=1))
    #   XWy <- tapply(weights * y, x, sum)
  # we need to check G and transfere it to a sparse matrix 
         if (is(G,"gra")) 
         {
         	if (nfv!=dim(G)[1]) stop("the dimesions of G should be equal with the levels of x")
         	G <- spam(as.numeric(G), nrow=nfv)
         }
      logsig2e <- log(sig2e) # getting logs
      logsig2b <- log(sig2b)
#------------- getting starting values -----------    
            beta <- solve(exp(-logsig2e)*XWX + exp(-logsig2b)*G, (exp(-logsig2e)*XWy)) 
          lambda <- sig2e/sig2b
              fv <- beta[x]
#-------------------------------------------------
#----------------- fitting -----------------------
# ------- if both fix just get Q -----------------
            if (sig2e.fix==TRUE &&sig2b.fix==TRUE)
                {
                out <- list()	
                par <- c(logsig2e, logsig2b)
            out$par <- par
         value.of.Q <- -Qf(par)
                 se <- NULL  	
                }
#------------ only sig^2_e is fixed ---------------
              if (sig2e.fix==TRUE &&sig2b.fix==FALSE) 
                {
                par <- log(sum((beta)^2)/N) # dO we divide by N?? 	
                Qf1 <- function(par)  Qf(c(logsig2e, par))  	
                out <- nlminb(start = par, objective = Qf1, lower = c(-20),  upper = c(20))
        out$hessian <- optimHess(out$par, Qf1)
         value.of.Q <- -out$objective
               shes <- 1/out$hessian
                se1 <- ifelse(shes>0, sqrt(shes), NA) 
                par <- c(logsig2e, out$par) 
         names(par) <- c("log(sige^2)", "log(sigb^2)")  
            out$par <- par    
                 se <- c(NA, se1)
                }
#------------ only sig^2_b is fixed ---------------
                 if (sig2e.fix==FALSE && sig2b.fix==TRUE)# sig^2_b fixed 
                {
               par <- log(sum(weights*(y-fv)^2)/N)   
               Qf2 <- function(par)  Qf(c(par, logsig2b))  	
               out <- nlminb(start = par, objective = Qf2, lower = c(-20),  upper = c(20))
       out$hessian <- optimHess(out$par, Qf2)
        value.of.Q <- -out$objective
              shes <- 1/out$hessian
               se1 <- ifelse(shes>0, sqrt(shes), NA) 
               par <- c( out$par, logsig2b) 
        names(par) <- c("log(sig2e)", "log(sig2b)")  
           out$par <- par    
                se <- c(se1, NA)               	
                }
#------------ both sig^2_e and sig^2_b are estimated ---------------
            if (sig2e.fix==FALSE && sig2b.fix==FALSE) # both estimated 
                {
               par <- c(logsig2e <- log(sum(weights*(y-fv)^2)/N), logsig2b <-log(sum((beta)^2)/N))
        names(par) <- c("log(sige^2)", "log(sigb^2)")   	
               out <- nlminb(start = par, objective = Qf, lower = c(-20, -20),  upper = c(20, 20))
       out$hessian <- optimHess(out$par, Qf)
        value.of.Q <- -out$objective
              shes <- try(solve(out$hessian))
                se <- if (any(class(shes)%in%"try-error")) rep(NA_real_,2) else sqrt(diag(shes)) 
        names(par) <- c("log(sig2e)", "log(sig2b)")  
              }
# end of fitting 
############################################################################################  
#------------ refitting ---------------------
             beta <- solve(exp(-out$par[1])*XWX + exp(-out$par[2])*G, (exp(-out$par[1])*XWy)) 
                    #  H <- solve(XWX + lambda*G, XWX)
      names(beta) <- levels(x)
               fv <- beta[x]
              tr1 <- sum(t(beta)%*%G%*%beta)/(exp(out$par[2]))    # PLEASE CHECK THIS 
  attributes(tr1) <- NULL
              tr2 <- nobs-(sum(weights*(y-fv)^2))/(exp(out$par[1]))#is the right trace if sig2e is fitted 
            sig2e <- exp(out$par[1])
            sig2b <- exp(out$par[2])
             sige <- sqrt(exp(out$par[1]))
             sigb <- sqrt(exp(out$par[2]))
attributes(sig2e) <- NULL
attributes(sig2b) <- NULL
 attributes(sige) <- NULL
 attributes(sigb) <- NULL
              var <- diag( solve(XWX +  (sig2e/sig2b)*G, XWX))
              var <- var[x] 
  #------------ If plot   ---------------------
#     if ((plot==TRUE) && (!is.null(polys))){
#       drawmap(fv,polys, col = gray(0.7), type = 'l')
#     }
  #------------ get the output ---------------------    # 
           fit <- list(fitted = fv, 
                           df = tr1,
      #                    df2 = tr2,
                   value.of.Q = value.of.Q,
                   deviance.Q = -2*value.of.Q,
                          par = list(par=(out$par), se=se),
                        sig2e = sig2e,
                        sig2b = sig2b,
                       lambda = sig2e/sig2b,
                         sige = sige,
                         sigb = sigb,
                            y = y,
                            x = x, 
                            G = G, 
                          var = var,
                       #    sumW = (sum(weights*(y-fv)^2)),
                           nobs =nobs,  
                       #    gaGga = sum(t(beta)%*%G%*%beta), 
                           beta = beta, # test                        
                      weights =  weights,
                            N = N,
                       method = "Q-function",
                         call = scall,
                          rss = sum(weights*(y-fv)^2),
                          aic = sum(weights*(-2*dNO(y, mu=fv, sigma=sqrt(exp(out$par[1])), log=TRUE)))+2*(tr1+1) , 
                          sbc = sum(weights*(-2*dNO(y, mu=fv, sqrt(exp(out$par[1])), log=TRUE)))+log(nobs)*(tr1+1),#
                     deviance = sum(weights*(-2*dNO(y, mu=fv, sigma=sqrt(exp(out$par[1])), log=TRUE))))
           class(fit) <- c( "MRF")
  fit
  }
# functions 4 disSmo methods
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# functions for MRF methods
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# methods for disSmo objects                    
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------                    
fitted.MRF<-function(object,...) 
{
  object$fitted
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
coef.MRF<-function(object,...) 
{
  coef <- object$par[["par"]]
  attr(coef, "se") <-  object$par[["se"]]
  coef
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
residuals.MRF<-function(object,type = c("z-scores", "simple"), ...) 
{
  type <- match.arg(type)
  res <- if (type == "z-scores") 
  {qNO(object$weights*pNO(object$y, mu=object$fitted, sigma=object$sige))}
  else
  {object$y-fitted(object)}
  res <- ifelse(res==-Inf, NA, res)
  res <- res[!is.na(res)]
  attr(res, "class") <- attr(object$y, "class")
  attr(res, "tsp") <- attr(object$y, "tsp")
  res
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
AIC.MRF <- function(object, ...,k=2,type = c("GD", "marginal"))
{
  type <- match.arg(type)
  if (type == "GD")
  {
    if (length(list(...))) {
      object <- list(object, ...)
      #isgamlss <- unlist(lapply(object, is.gamlss))
      #if (!any(isgamlss)) 
      #    stop("some of the objects are not gamlss")
      df <- as.numeric(lapply(object, function(x) x$df))
      AIC <- as.numeric(lapply(object, function(x) deviance(x) + 
        x$df * k))
      val <- as.data.frame(cbind(df, AIC))
      Call <- match.call()
      Call$k <- NULL
      row.names(val) <- as.character(Call[-1])
      val <- val[order(AIC), ]
      val
    }
    else 
    {
      val <- if (is(object, "MRF")) 
        object$deviance + (object$df) * k
      else stop(paste("this is not a disSmo object"))
    }   
  }  
  if (type == "marginal")
  {
    if (length(list(...))) {
      object <- list(object, ...)
      #isgamlss <- unlist(lapply(object, is.gamlss))
      #if (!any(isgamlss)) 
      #    stop("some of the objects are not gamlss")
      #md <- as.numeric(lapply(object, function(x) length(coef(x, type="marginal"))))
      df <- as.numeric(lapply(object, function(x) x$mdf))
      AIC <- as.numeric(lapply(object, function(x) deviance(x, type="marginal") + x$mdf * k))
      val <- as.data.frame(cbind(df, AIC))
      Call <- match.call()
      Call$k <- NULL
      Call$type <- NULL
      row.names(val) <- as.character(Call[-1])
      val <- val[order(AIC), ]
      val
    }
    else 
    {
      val <- if (is(object, "MRF")) 
        object$deviance.Q + 2 * k
      else stop(paste("this is not a disSmo object"))
    }   
  }
  val
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
deviance.MRF<-function(object,type = c("GD", "marginal"), ...) 
{
  type <- match.arg(type)
  if (type == "GD") return(object$deviance)
  else return(object$deviance.Q)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
plot.MRF <- function (x, xvar=NULL, parameters=NULL, ts=FALSE, summaries=TRUE, ...) 
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
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
print.MRF  <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{   
  cat("\nMarkov Random Fields fit", "\n")
  cat("Fitting method:", deparse(x$method), "\n")
  cat("\nCall: ", deparse(x$call), "\n", fill = TRUE)
  cat("Coefficients")
  co <- coef(x)
  cat("\n  ", names(co), "\n")
  cc <-simplify2array(co, higher=TRUE)
  cat(cc, " \n")
  cat("\n Degrees of Freedom for the fit:", x$df, "Residual Deg. of Freedom  ", 
      x$N-x$df, "\n")
  cat("Global Deviance:    ", format(signif(x$deviance)), 
      "\n            AIC:    ", format(signif(x$aic)), "\n            SBC:    ", 
      format(signif(x$sbc)), "\nMarginal Devia.:    ",
      format(signif(x$deviance.Q))  
  )
  invisible(x)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
summary.MRF  <- function (object, digits = max(3, getOption("digits") - 3), ...) 
{   
  cat("\nMarkov Random Fields fit\n")
  cat("Fitting method:", deparse(object$method), "\n")
  cat("\nCall: ", deparse(object$call), "\n", fill = TRUE)
  coef  <- coef(object)
  se.coef <- attributes(coef)$se
  tval <- coef/se.coef
  matcoef <- cbind(coef, se.coef, tval, 2*(1-pnorm(abs(tval))))
  dimnames(matcoef) <- list(names(tval), c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
  cat("\nCoefficient(s):\n")
  printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
  
  
  cat("\n Degrees of Freedom for the fit:", object$df, "Residual Deg. of Freedom  ", 
      object$N-object$df, "\n")
  cat("Global Deviance:    ", format(signif(object$deviance)), 
      "\n            AIC:    ", format(signif(object$aic)), "\n            SBC:    ", 
      format(signif(object$sbc)), "\nMarginal Devia.:    ",
      format(signif(object$deviance.Q))  
      )
  invisible(object)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
logLik.MRF <- function (object, ...) 
{
  val <- -object$deviance/2
  attr(val, "nall") <- object$N
  attr(val, "nobs") <- object$noObs
  attr(val, "df") <- object$df
  class(val) <- "logLik"
  val
}
# methods are finish here    
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
predict.MRF <- function(object, newdata=NULL, ...)
{
if (is.null(newdata)) pred  <- fitted(object)
         else
         {# check if levels are in x levels
           newdata <- as.character(newdata)
           if (!all(newdata%in%levels(object$x))) stop("newdata contains levels not existing in x")
          pred<-object$beta[newdata]
         }
pred
}

#----------------------------------------------------------------------------------------
