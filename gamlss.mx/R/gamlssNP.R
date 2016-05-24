#-------------------------------------------------------------------------------
# TO DO
# i) Ardo residuals ??
#-------------------------------------------------------------------------------
gamlssNP<-function(formula,
                    random =~1,
                    family = NO(),
                      data = NULL,
                         K = 4,
                   mixture = c("np","gq"),
                       tol = 0.5,
                    weights,# we have to decide how this will be in Randon eff
                    pluginz,
                   control = NP.control(...),
                 g.control = gamlss.control(trace=FALSE), 
                    ...)
{
##  Nonparametric ML/Gaussian quadrature macros for maximum likelihood
## in GAMLSS
## R code originally by Ross Darnell (2002), modifications and extensions
## by Jochen Einbeck / John Hinde (2005).
## modified for gamlss by Mikis Stasinopoulos Thursday, July 20, 2006 
#-------------------------------------------------------------------------------
## functions within gamlssNP
#-------------------------------------------------------------------------------
.gamlss.bi.list <- eval(quote(.gamlss.bi.list), envir = getNamespace("gamlss"))
## expanding the data ----------------------------------------------------------
expand.vc <- function(x,ni)
{
  if (length(ni)==1)
  {
      if (ni<1) stop("ni should be greater than 1")
      xx <- x
      if (ni==1) return(xx)
       for ( i in 2:ni) xx <- rbind(x,xx)
      xx
   } 
   else 
   {
      n <- dim(x)[[1]]
      c <- dim(x)[[2]]
     xx <- NULL
      for ( i in seq(1,n)){
          xx <- rbind(xx,matrix(rep(x[i,],ni[i]),ncol=c,byrow=TRUE))
      }
      xx
  }
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
weightslogl.calc.w <- function(p,fjk,weights){
  # p is a vector of length K containing the mixture proportions
  # fjk is a JXK matrix of log densities
  logpf <- t(apply(fjk,1,"+",log(p)))
     pf <- exp(logpf)# p*fjk
    Spf <- as.vector(apply(pf,1,sum))#weights denominator length N
      w <- pf/Spf
 ML.dev <- -2*sum(weights*log(Spf))
  #l <- sum(log(Spf)) # log likelihood
  list(w=w,ML.dev=ML.dev)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
get.log.likelihood.firstTime <- function(obj,mu, ...)
{
if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
   fname <- obj$family[[1]]
 DistPar <- obj$parameters
   nopar <- length(DistPar) 
    dfun <- paste("d",fname,sep="")
 # the first time the length of y is N not NxK so we recalculate y
    if (fname%in%.gamlss.bi.list) 
        {
        yy <- rep(obj$y, K)
        bd <- rep(obj$bd, K) 
        } 
    else 
        {
        yy <- rep(obj$y,K) 
        }   
        if (("sigma"%in%obj$parameters)) sigma <- rep(fitted(obj,"sigma"),K) 
        if (("nu"%in%obj$parameters)) nu <- rep(fitted(obj,"nu"),K) 
        if (("tau"%in%obj$parameters)) tau <- rep(fitted(obj,"tau"),K)  
 switch(nopar,  
    {lik <- if (fname%in%.gamlss.bi.list) eval(call(dfun, x=yy, bd=bd,  mu=mu, log=TRUE))
            else   eval(call(dfun, x=yy, mu=mu, log=TRUE))},
    {lik <- if (fname%in%.gamlss.bi.list) eval(call(dfun, x=yy, bd=bd,  mu=mu, sigma=sigma, log=TRUE ))
            else   eval(call(dfun, x=yy, mu=mu, sigma=sigma, log=TRUE)) },
    {lik <-  if (fname%in%.gamlss.bi.list) eval(call(dfun, x=yy, bd=bd,  mu=mu, sigma=sigma, nu=nu, log=TRUE ))
            else eval(call(dfun,x=yy, mu=mu, sigma=sigma, nu=nu ,log=TRUE))},
    {lik <- if (fname%in%.gamlss.bi.list) eval(call(dfun, x=yy, bd=bd,  mu=mu, sigma=sigma, nu=nu, tau=tau, log=TRUE ))
            else eval(call(dfun,x=yy, mu=mu, sigma=sigma, nu=nu, tau=tau,log=TRUE))})
lik
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
get.log.likelihood <- function(obj, ...)
{
if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
   fname <- obj$family[[1]]
 DistPar <- obj$parameters
   nopar <- length(DistPar) 
    dfun <- paste("d",fname,sep="")
 switch(nopar,  
    { lik <- if (fname%in%.gamlss.bi.list) eval(call(dfun,x=obj$y, bd=obj$bd,  mu=fitted(obj,"mu"), log=TRUE))
             else   eval(call(dfun, x=obj$y, mu=fitted(obj, "mu"), log=TRUE))},
    {lik <- if (fname%in%.gamlss.bi.list)  eval(call(dfun,x=obj$y, bd=obj$bd,  mu=fitted(obj,"mu"), sigma=fitted(obj,"sigma"), log=TRUE ))
             else   eval(call(dfun, x=obj$y, mu=fitted(obj, "mu"), sigma=fitted(obj,"sigma"), log=TRUE)) },
    {lik <- if (fname%in%.gamlss.bi.list)  eval(call(dfun,x=obj$y, bd=obj$bd,  mu=fitted(obj,"mu"), sigma=fitted(obj,"sigma"), nu=fitted(obj,"nu"),log=TRUE ))
             else   eval(call(dfun,x=obj$y, mu=fitted(obj,"mu"), sigma=fitted(obj,"sigma"), nu=fitted(obj,"nu"),log=TRUE))},
    {lik <- if (fname%in%.gamlss.bi.list)  eval(call(dfun,x=obj$y, bd=obj$bd,  mu=fitted(obj,"mu"), sigma=fitted(obj,"sigma"), nu=fitted(obj,"nu"), tau=fitted(obj,"tau"),log=TRUE))
             else eval(call(dfun,x=obj$y, mu=fitted(obj,"mu"), sigma=fitted(obj,"sigma"), nu=fitted(obj,"nu"), tau=fitted(obj,"tau"),log=TRUE))})
lik
}                          
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# for getting the cumulative function
get.the.p.function <- function(object, ...)
{
if (!is.gamlss(object))  stop(paste("This is not an gamlss object", "\n", ""))
   fname <- object$family[[1]]
 DistPar <- object$parameters
   nopar <- length(DistPar) 
    dfun <- paste("p",fname,sep="")
 # binomial denominators
 switch(nopar,  
    {pfun <- if (fname%in%.gamlss.bi.list)  eval(call(dfun,q=object$y, bd=object$bd,  mu=fitted(object,"mu")))
             else  eval(call(dfun,q=object$y, mu=fitted(object,"mu")))},
    {pfun <- if (fname%in%.gamlss.bi.list)  eval(call(dfun,q=object$y, bd=object$bd,  mu=fitted(object,"mu"),sigma=fitted(object,"sigma")))
             else  eval(call(dfun,q=object$y, mu=fitted(object,"mu"), sigma=fitted(object,"sigma"))) },
    {pfun <- if (fname%in%.gamlss.bi.list)  eval(call(dfun,q=object$y, bd=object$bd,  mu=fitted(object,"mu"),sigma=fitted(object,"sigma"),  nu=fitted(object,"nu")))
             else eval(call(dfun,q=object$y, mu=fitted(object,"mu"), sigma=fitted(object,"sigma"), nu=fitted(object,"nu")))},
    {pfun <-  if (fname%in%.gamlss.bi.list)  eval(call(dfun,q=object$y, bd=object$bd,  mu=fitted(object,"mu"),sigma=fitted(object,"sigma"),  nu=fitted(object,"nu"),  tau=fitted(object,"tau")))
             else eval(call(dfun,q=object$y, mu=fitted(object,"mu"), sigma=fitted(object,"sigma"), nu=fitted(object,"nu"), tau=fitted(object,"tau")))})
pfun 
}             
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## The gamlssNP starts here 
#library(gamlss)
  mixture <- match.arg(mixture)
## getting the contol papameters
  EMn.cyc <- control$EMn.cyc
     EMcc <- control$EMcc
  verbose <- control$trace
 plot.opt <- control$plot.opt
     damp <- control$damp
#----
#gamlss.bi.list <- c("BI", "Binomial", "BB", "Beta Binomial") # binomial denominators
     call <- match.call()
# we need the data so we can expand them ---------------------------------------
 if (missing(data)) stop("the data argument is needed")
 if(!missing(data) & any(is.na(data))) 
               stop("The data contains NA's, use data = na.omit(mydata)") 
   family <- as.gamlss.family(family)     
        Y <- model.extract(model.frame(formula,data=data),"response")
        N <- if(is.null(dim(Y))) length(Y) else dim(Y)[1]   # dimension for y     
## prior weights
 pweights <- if (missing(weights)) rep(1,N) else weights  
     data <- data.frame(data, pweights)   
## initial fit and simple gamlss for K=1 -----------------------------------------------
   fitout <- gamlss(formula, family=family, weights=pweights, data = data, 
                                                           control=g.control)#    
       l0 <- length(coef(fitout,"mu"))   
  if (K == 1) {  return(fitout)   } # finish if K =1 
## Expand the response
    datak <- expand.vc(data,K)# expand data.
   kindex <- rep(1:K,rep(N,K))# index for the mixtures
     tmp1 <- gqz(K, minweight = 1e-14)
       z0 <- -tmp1$location
        z <- rep(-tmp1$location, rep(N, K))
 #   p <- tmp$w
 #    tmp <- hermite(K)# grab weights and abscissas
 #     z0 <- tmp$z    # for GQ - masspoints
 #      z <- rep(tmp$z,rep(N,K))
  #pweights <- rep(pweights,K)
       p <- tmp1$weight
   rform <- random
##  Generate the random design matrix and append to fixed matrix
   mform <- strsplit(as.character(random)[2],'\\|')[[1]]
   mform <- gsub(' ', '',mform)
## the length(mform) determines whether VC or not 
#--------------for VC models-------------------- 
  if (length(mform)==2)
    { 
  mform1 <- mform[1]
  mform2 <- mform[2]
   group <- factor(levels(as.factor(datak[,mform2])))##R.E.D. 15/2/06 
      nr <- nlevels(group)
 ijindex <- rep(1:N,K)
 groupij <- data[,mform2]
groupijk <- rep(groupij,K)
    }
   else {mform1<-mform} 
#------------------------------------------------
  if (mixture=='np')
      {
      # Nonparametric random effect
      datak$MASS <- gl(K,N) 
      if (mform1=='1') random <- formula(~MASS) # random <- formula(~MASS-1) 
      else {  
          # Nonparametric random coefficient
          #random <- formula(paste('~ MASS + ',paste(mform1,'MASS',sep=":",collapse='+'), '-1',sep='')) 
          random <- formula(paste('~ MASS + ',paste(mform1,'MASS',sep=":",collapse='+'),sep=""))   
           }
      } 
   else 
      {
      datak$z <- z
      if (mform1=='1') 
        {
        # datak$z <- z
          random <- formula(~ z )
        } 
      else 
        { 
        ## this has to be checked ???????
            random <- formula(paste('~ z + ',paste(mform1,'z',sep=":",collapse='+'),sep=""))
        }
      }
## I do not have provision to use different starting values as in gamlssMX  
if (missing(pluginz))
   {
    sz <- if ("sigma"%in%fitout$parameters)  tol* fitted(fitout,"sigma")*z else  tol*z 
   }
else 
   {
    sz <- if ("sigma"%in%fitout$parameters)    rep(pluginz-fitout$mu.coef[[1]],rep(N,K)) #?? is this correct  
   } 
Eta <- fitout$mu.lp + sz
## I have to check what is doing  
 if (mixture=="np")
  {
      tol <- max(min(tol,1),1-damp)  #For tol >  1 or damp=F no Damping
      if(length(fitout$coef)==1){followmass<-matrix(Eta[(1:K)*N],1,K)} else {followmass<-matrix(fitout$mu.coef[1]+sz[(1:K)*N],1,K)}
  } else {
      followmass<-NULL; tol<-1
  }
      Mu <- family$mu.linkinv(Eta)
## expanded fitted values
    logf <- get.log.likelihood.firstTime(fitout,Mu)
    logf <- ifelse(logf>-700,logf,-700)
## Calculate the weights from initial model
   if (length(mform)==2) # if repeated measurments
    { 
   groupk <- interaction(groupijk,factor(kindex))
   logmik <- matrix(tapply(logf,groupk,sum),nrow=nr,ncol=K)
   #hoehle 09.04.2010: Problem if there are no instances of a specific
   #groupk factor level -> NA. But we don't want NAs
   logmik[is.na(logmik)] <- 0
      tmp <- weightslogl.calc.w(p,logmik,pweights[1:nr])# this maybe is rubish and the prior weighrs should be in the logf calulation
       ww <- as.vector(tmp$w[match(groupij,group),])
       ww <- ifelse(ww<1e-20,0,ww)  
    }
    else
    {  
      tmp <- weightslogl.calc.w(p,matrix(logf,ncol=K),rep(pweights,K))
       ww <- as.vector(tmp$w)
       ww <- ifelse(ww<1e-20,0,ww)
    } 
## add weights  into the datak    
    datak <- data.frame(datak, ww)
## add the MASS into the formula    
       ex <- parse(text=paste(deparse(formula[[3]]),deparse(random[[2]]),sep="+"))
  formula[[3]] <- ex[[1]] #
  # Initialize for EM loop
   ML.dev <- ML.dev0 <- deviance(fitout)
     iter <- ml <- 1
converged <- FALSE
                                                                                
 while (iter <= EMn.cyc && (!converged || (iter<=9 && mixture=='np'  )))
  {
  ##########Start of EM ##############################
  if (verbose){ if (iter%%17==16) cat(iter, "..\n") else cat(iter,"..")  }
# the gamlss fitting  
fitout <- if (iter==1)
          { gamlss(formula, family=family, weights=ww*pweights, data = datak, control=g.control, ...)}
        else
          { gamlss(formula, family=family, weights=ww*pweights, data = datak, control=g.control, start.from=fitout,  ...)}  
# Save the EM Trajectories
   if (mixture=="np")
      {   
      attcoef <- attr(coef(fitout),"names")
      if (mform1=='1')
       {
      masspoint <-  c(coef(fitout)[ grep("(Intercept)",attcoef)], coef(fitout)[grep("(Intercept)",attcoef)]+
                                    coef(fitout)[grep("MASS",attcoef)] )# 
       }
     else
       {
       sdif<-setdiff( grep("MASS",attcoef),grep(mform1,attcoef))
      masspoint <- c(coef(fitout)[ grep("(Intercept)",attcoef)], coef(fitout)[grep("(Intercept)",attcoef)]+
                                    coef(fitout)[sdif] )# 
       }
      followmass <- rbind(followmass, masspoint)
      }
   
## Calculate loglikelihood for expanded model for this iteration
     logf <- get.log.likelihood(fitout)  
     logf <- ifelse(logf>-740,logf,-740)  #avoid zero weights
     logf <- matrix(logf, ncol=K)
## Calculate the component proportions from the weights
    if (mixture=='np') p <- as.vector(tapply(ww,kindex,mean))
      # Calculate updated weights and loglikehood
    if (length(mform)==2) # if repeated measurments
    { 
       logmik <- matrix(tapply(logf,groupk,sum),nrow=nr,ncol=K)
       #hoehle 09.04.2010: Problem if there are no instances of a specific
       #groupk factor level -> NA. But we don't want NAs
       logmik[is.na(logmik)] <- 0
       tmp <- weightslogl.calc.w(p,logmik,pweights[1:nr])
        ww <- as.vector(tmp$w[match(groupij,group),])
        ww <- ifelse(ww<1e-20,0,ww)  
    }
    else
    {  
       tmp <- weightslogl.calc.w(p,matrix(logf,ncol=K),pweights[1:N])
        ww <- as.vector(tmp$w)
        ww <- ifelse(ww<1e-20,0,ww)
    } 
  datak$ww <- ww
      ML.dev[iter+1] <- tmp$ML.dev# -2 * log L max
      if (ML.dev[iter+1]>ML.dev0) {ml<-ml+1}  #only relevant for graphical output
      converged <- abs(ML.dev[iter+1] - ML.dev[iter])< EMcc
      iter <- iter + 1

   }###########################End of EM loop#############
  if (verbose)
  {
      cat("\n")
      if (converged)
      {
          cat("EM algorithm met convergence criteria at iteration  ", iter-1,"\n")
      } 
      else 
      {
        cat("EM algorithm failed to meet convergence criteria at iteration # ",
        iter-1,"\n")  
      }
  }
  mass.points <- masses <- NULL
    np <- length(fitout$mu.coef)
   ebp <- apply(ww*matrix(fitout$mu.lp,N,K,byrow=FALSE),1,sum)  #Emp. Bayes Pred. (Aitkin, 96)
  #m <- seq(1,np)[substr(attr(fitout$mu.coefficients,'names'),1,4)=='MASS']
  #mass.points <- masspoint c(coef(fitout)[ grep("(Intercept)",attcoef)], coef(fitout)[grep("(Intercept)",attcoef)]+
  #                                  coef(fitout)[grep("MASS",attcoef)] )# #fitout$mu.coefficients[m]
 # I am not sure what this does
 # if (is.na(fitout$mu.coefficients[np])){length(fitout$mu.coefficients)<-np-1}# if one variable is random and fixed
  if (plot.opt==3 && mixture=="np")
  {
      op<-par(mfrow=c(2,1),cex=0.5,cex.axis=1.5,cex.lab=1.5)
      on.exit(par(op))
  }
  if (plot.opt==1|| plot.opt==3)
  {
      if  ((family$type=="continuous" ) && damp  && mixture=='np' && iter>=max(8,ml+1))
      {
          #Linear interpolation for initial cycles
          ML.dev[2: max(7,ml)]<-ML.dev0+ 1:max(6,ml-1)/ max(7,ml)*(ML.dev[max(8,ml+1)]-ML.dev0) 
      }  
      plot(0:(iter-1),ML.dev, col=1,type="l",xlab='EM iterations',ylab='-2logL')
      if (verbose){ cat("Global deviance trend plotted.\n")}
  }
# if ( mform=='1'){
#            ylim<- c(min(na.omit(R)),max(na.omit(R)))
#            if (ylim[1]==-Inf){ylim[1]<-min(followmass[,])} ;if (ylim[2]==Inf){ylim[2]<-max(followmass[,])}
#      } else  { 
#          ylim<-c(min(followmass[,]),max(followmass[,]))
#      }  
    if (mixture=="np")
     {
      ylim<-c(min(followmass[,], na.rm = TRUE),max(followmass[,], na.rm = TRUE))
       if (plot.opt==2|| plot.opt==3)
         {
            plot(0:(iter-1),followmass[,1],col=1,type='l',ylim=ylim,ylab='mass points',xlab='EM iterations')
            for (i in 1:K)
            { lines(0:(iter-1), followmass[,i],col=i)
                        #if (mform=='1'){ points(rep(iter-1,length(R)),R)}
            }
            if (verbose){ cat("EM Trajectories plotted.\n")}
          }
     }
## 
## residuals ----
 # there is a problem with binomial also for "qp" prob is fixed 
              
               prob <- apply(matrix(ww,ncol=K),2,mean)
                 WF <- matrix(get.the.p.function(fitout), ncol=K, nrow=N) 
                 for (i in 1:K)     WF[,i] <- prob[i]*WF[,i] 
                res <- qnorm(rowSums(WF))
#               res2 <- apply(ww*matrix(fitout$mu.lp,N,K,byrow=FALSE),1,sum) 
#get.the.p.function(fitout)
#WF <- matrix(0, ncol=K, nrow=N)
# for (i in 1:K)     WF[,i] <- prob[i]*get.the.p.function(allModels[[i]]) 
#res <-qnorm(rowSums(WF))
#----------------
      fitout$df.fit <- if (mixture=="np") fitout$df.fit+K-1 else fitout$df.fit
 fitout$df.residual <- N-fitout$df.fit 
  fitout$G.deviance <- ML.dev[iter]
        fitout$call <- call 
     fitout$formula <- formula 
      fitout$random <- rform
 #fitout$mass.points <- mass.points
   fitout$post.prob <- list(matrix(ww,nrow=N,byrow=FALSE))
       fitout$data  <- data
      # fitout$dataK <- datak
         fitout$ebp <- ebp
      fitout$EMiter <- iter - 1
    fitout$pweights <- pweights
 fitout$EMconverged <- converged
         fitout$aic <- fitout$G.deviance+2*fitout$df.fit 
         fitout$sbc <- fitout$G.deviance+log(N)*fitout$df.fit 
fitout$allresiduals <- fitout$residuals
   fitout$residuals <- res
           fitout$N <- N  
           fitout$K <- K
       fitout$type  <- "Mixture"  
 if (mixture=="np") 
    {
      if (mform1=='1')
       {
 fitout$mass.points <- masspoint
       }
     else
       {
                 lf <- length(mform1)
      MASSPOINT     <- matrix(0, ncol=lf+1, nrow=K)
      MASSPOINT[,1] <- masspoint
        for (i in 1:lf)
           {  
           mm <- grep( mform1[i],attcoef)
      MASSPOINT[,i+1] <- c(coef(fitout)[mm[1]], coef(fitout)[mm[1]]+
                                    coef(fitout)[mm[2:K]] )#
           } 
 fitout$mass.points <- MASSPOINT   
       }
         masses     <- prob 
     names(masses)  <- paste('MASS',1:K,sep='')
    fitout$prob     <- masses
 fitout$orig.family <- fitout$family
     fitout$family <- paste( fitout$family, "Mixture with NP")             
   #   class(fitout) <- list("gamlssNP", "gamlss")
    } 
 else 
    {
  fitout$mass.points <- fitout$coef[1]+fitout$coef[np]*z0
  fitout$orig.family <- fitout$family
     fitout$family   <- paste( fitout$family[[1]], "Mixture with NO")           
         fitout$prob <- list(tmp1$weight)          
   #   class(fitout) <- list("gamlssNO", "gamlss")   
    }
      class(fitout) <- list("gamlssNP", "gamlss")
  fitout
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#get.log.likelihood <- function(obj,y,mu, ...)
#{
#if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
#   fname <- obj$family[[1]]
# DistPar <- obj$parameters
#   nopar <- length(DistPar) 
#    dfun <- paste("d",fname,sep="")
# switch(nopar,  
#    {lik <- eval(call(dfun,y=y, mu=mu, log=TRUE))},
#    {lik <- eval(call(dfun,y=y, mu=mu, sigma=fitted(obj,"sigma"),log=TRUE)) },
#    {lik <- eval(call(dfun,y=y, mu=mu, sigma=fitted(obj,"sigma"), nu=fitted(obj,"nu"),log=TRUE))},
#    {lik <- eval(call(dfun,y=y, mu=mu, sigma=fitted(obj,"sigma"), nu=fitted(obj,"nu"), tau=fitted(obj,"tau"),log=TRUE))})
#lik
#}             
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
NP.control <- function(EMcc = 0.001, EMn.cyc = 200, damp = TRUE, trace = TRUE, plot.opt = 3, ...)
{
        if(EMcc <= 0) {
warning("the value of cc supplied is zero or negative the default value of 0.001 was used instead")
                c.crit <- 0.001}
        if(EMn.cyc < 1) {
warning("the value of no cycles supplied is zero or negative the default value of 200 was used instead")
                n.cyc <- 200}
        if(is.logical(trace)) trace <- trace 
        else if (is.numeric(trace) & trace <= 0) 
              {warning("the value of trace supplied is less or equal t zero the default of 1 was used instead")
                trace <- 1
              }    
        list(EMcc = EMcc, EMn.cyc = EMn.cyc, trace = as.logical(trace)[1], damp=as.logical(damp)[1], plot.opt=plot.opt)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
################################################################################
#-------------------------------------------------------------------------------
 # Auxiliary functions:
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
print.gamlssNP <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    K <- x$K 
    cat("\nMixing Family: ", deparse(x$family), "\n", fill=TRUE)
    cat("Fitting method: EM algorithm \n")
    cat("\nCall: ", deparse(x$call),  "\n", fill=TRUE)
        if ("mu" %in% x$parameters) 
         {      
          cat("Mu Coefficients : \n")   
          print.default(format(coef(x, "mu"), digits = digits), print.gap = 2, quote = FALSE)        
         }
        if ("sigma" %in% x$parameters) 
         { 
          cat("Sigma Coefficients :\n")           
          print.default(format(coef(x, "sigma"), digits = digits), print.gap = 2, quote = FALSE)        
         }
        if ("nu" %in% x$parameters) 
         {       
         cat("Nu Coefficients : \n")    
         print.default(format(coef(x, "nu"), digits = digits), print.gap = 2, quote = FALSE)        
         }
        if ("tau" %in% x$parameters) 
         {
          cat("Tau Coefficients : \n")          
          print.default(format(coef(x, "tau"), digits = digits), print.gap = 2, quote = FALSE)        
         }
       if(!is.list(x$prob)) 
        {cat("\nEstimated probabilities:", x$prob, "\n" )}
    cat("\nDegrees of Freedom for the fit:", x$df.fit, "Residual Deg. of Freedom  ", 
        x$df.residual, "\n")
   cat("Global Deviance:    ", format(signif(x$G.deviance)), 
        "\n            AIC:    ", format(signif(x$aic)), "\n            SBC:    ", 
        format(signif(x$sbc)), "\n")
    invisible(x)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# This need attention MS
fitted.gamlssNP<-function (object, K=1, what = c("mu", "sigma", "nu", "tau"), ... ) 
{  
  what <- match.arg(what)
if (K%in%seq(1:object$K))
 {
 index <- seq(1+(K-1)*object$N, K*object$N) 
     x <- object[[paste(what,"fv",sep=".")]][index] #fitted(object$mu.fv, ...)[index]
 }
else 
  { # ?? THIS IS NOT WORKING  MS 
   x <-  object[[paste(what,"fv",sep=".")]]
  }
 x
}
#-------------------------------------------------------------------------------
family.gamlssNP <- function(object, ...) {
     object$family
 }
#-------------------------------------------------------------------------------
# this is the plot.gamlss function 
# created by PA  May 2002
# last change by MS Friday, Wednesday, December 17, 2003 at 09:11
# to incoorporate options in the plotying parameters
# the following options have been used for the BCT paper 
# par(mfrow=c(2,2), mar=par("mar")+c(0,1,0,0), col.axis="blue4", col="blue4", col.main="blue4",col.lab="blue4",pch="+",cex=.45, cex.lab=1.2, cex.axis=1, cex.main=1.2  )
#-------------------------------------------------------------------------------
plot.gamlssNP<- function (x, xvar=NULL, parameters=NULL, ts=FALSE, summaries=TRUE, ...) 
{
    if (!is.gamlss(x))  stop(paste("This is not an gamlss object", "\n", ""))
# chech for the residuals 
    if (is.null(x$residuals)) #    
         stop(paste("There are no randomised quantile residuals"))
# whether index or x-variable
    residx <- resid(x) # get the residuals 
         w <- x$weights
    xlabel <- if(!missing(xvar)) deparse(substitute(xvar)) else deparse(substitute(index))
    if(is.null(xvar))  xvar <- seq(1,length(resid(x)),1) # MS
# plot parameters
    if(is.null(parameters))
          op <- par(mfrow=c(2,2), mar=par("mar")+c(0,1,0,0), col.axis="blue4", col.main="blue4", col.lab="blue4",  col="darkgreen", bg="beige" )
    else  op <- parameters 
    # top  figures 
    if(identical(ts, TRUE))
     {
    # require(stats)
     acf.new<-acf(residx,plot=FALSE)
     plot(acf.new,xlim=c(2,length(acf.new$acf)),ylim=range(acf.new$acf[-1]))   # ms Tuesday, August 19, 2003 at 11:04
     pacf(resid(x))
     }
     else 
     {
     fittedvalues <- if(is.null(fitted(x))) fitted(x,"sigma") else fitted(x) # MS Wednesday, September 10, 2003 at 21:20
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
            lines(residx, residx, col="red" , lwd=.4, cex=.4 )
  
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
residuals.gamlssNP <- function(object,...)
{
object$residuals
}
#-------------------------------------------------------------------------------
#par(mfrow=c(2,2), mar=par("mar")+c(0,1,0,0), col.axis="blue4", col="blue4", col.main="blue4",col.lab="blue4",pch="+",cex=.45, cex.lab=1.2, cex.axis=1, cex.main=1.2  )
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# this is  Nick Sofroniou function from the package 
# npmlreg of  Jochen Einbeck, Ross Darnell and John Hinde (2006).
 gqz <- function (numnodes = 20, minweight = 1e-06) 
{
    out <- gauss.quad(numnodes, "hermite")
    h <- rbind(out$nodes * sqrt(2), out$weights/sum(out$weights))
    ord <- order(h[1, ], decreasing = TRUE)
    h <- h[, ord]
    h <- cbind(h[1, ], h[2, ])
    h <- subset(as.data.frame(h), (h[, 2] >= minweight))
    names(h) <- c("location", "weight")
    h
}
#-------------------------------------------------------------------------------
# this  Gordon Smyth function from package 
# package:statmod 
# see there for documantation
gauss.quad <- function (n, kind = "legendre", alpha = 0, beta = 0) 
{
    n <- as.integer(n)
    if (n < 0) 
        stop("need non-negative number of nodes")
    if (n == 0) 
        return(list(nodes = numeric(0), weights = numeric(0)))
    kind <- match.arg(kind, c("legendre", "chebyshev1", "chebyshev2", 
        "hermite", "jacobi", "laguerre"))
    i <- 1:n
    i1 <- i[-n]
    switch(kind, legendre = {
        muzero <- 2
        a <- rep(0, n)
        b <- i1/sqrt(4 * i1^2 - 1)
    }, chebyshev1 = {
        muzero <- pi
        a <- rep(0, n)
        b <- rep(0.5, n - 1)
        b[1] <- sqrt(0.5)
    }, chebyshev2 = {
        muzero <- pi/2
        a <- rep(0, n)
        b <- rep(0.5, n - 1)
    }, hermite = {
        muzero <- sqrt(pi)
        a <- rep(0, n)
        b <- sqrt(i1/2)
    }, jacobi = {
        ab <- alpha + beta
        muzero <- 2^(ab + 1) * gamma(alpha + 1) * gamma(beta + 
            1)/gamma(ab + 2)
        a <- i
        a[1] <- (beta - alpha)/(ab + 2)
        i2 <- 2:n
        abi <- ab + 2 * i2
        a[i2] <- (beta^2 - alpha^2)/(abi - 2)/abi
        b <- i1
        b[1] <- sqrt(4 * (alpha + 1) * (beta + 1)/(ab + 2)^2/(ab + 
            3))
        i2 <- i1[-1]
        abi <- ab + 2 * i2
        b[i2] <- sqrt(4 * i2 * (i2 + alpha) * (i2 + beta) * (i2 + 
            ab)/(abi^2 - 1)/abi^2)
    }, laguerre = {
        a <- 2 * i - 1 + alpha
        b <- sqrt(i1 * (i1 + alpha))
        muzero <- gamma(alpha + 1)
    })
    A <- rep(0, n * n)
    A[(n + 1) * (i - 1) + 1] <- a
    A[(n + 1) * (i1 - 1) + 2] <- b
    A[(n + 1) * i1] <- b
    dim(A) <- c(n, n)
    vd <- eigen(A, symmetric = TRUE)
    w <- rev(as.vector(vd$vectors[1, ]))
    w <- muzero * w^2
    x <- rev(vd$values)
    list(nodes = x, weights = w)
}
#-------------------------------------------------------------------------------
