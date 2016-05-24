# Authors Mikis Stasinopoulos Bob Rigby and Vlasios Voudouris
# created 11-04-12
# upadated  8--6-14
#-------------------------------------------------------------------------------
# what is new June 2014
# i)  the power transformation (x^p) function pt() now is defined to include 
#      zero as log(x)
# ii) the power transformation uses GAIC instead GD which was not reliable 
#      since different transformation will use differt effective df's 
# iii) a prediction function for an lms object is created
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#This function is design to help the user to construct centile estimation.
#It is only applicable when only "one" explanatory variable is available 
#   (usually age).
#It stats by fitting a normal error distribution and smooth function for mu and 
# then proceeds by fitting #several appropriate distributions.
# The set of gamlss.family distribution to fit are specified in the argument
#   families. 
#The default families arguments is LMS=c("BCCGo",  "BCPEo", "BCTo") that is the
# LMS class of distributions.
#Note that this class is only appropriate when y is positive (with no zeros). 
# If the response variable contains negative values and/or zeros then use 
# the argument theSHASH theSHASH <-  c("NO", "SHASHo") or add any other 
#  distribution which you think is appropriate
#-------------------------------------------------------------------------------
#  the LMS familily of distributions
      LMS <- c("BCCGo",  "BCPEo", "BCTo")
# the SHASH
 theSHASH <-  c("NO", "SHASHo")
#-------------------------------------------------------------------------------
#------------------------------------------------------------------------------- 
lms <- function(y, x,
        families = LMS,
            data = NULL, 
               k = 2, # for the AIC
            cent = 100*pnorm((-4:4)*2/3),
     calibration = TRUE,
         trans.x = FALSE,
       fix.power = NULL,
       lim.trans = c(0, 1.5),   
            prof = FALSE,
            step = 0.1, 
          legend = FALSE,
           mu.df = NULL,
        sigma.df = NULL,
           nu.df = NULL,
          tau.df = NULL,
          c.crit = 0.01,
       method.pb = c("ML", "GAIC"),
              ... 
                )  
{
#-------------------------------------------------------------------------------
# local function 
findPower <- function(y, x, data = NULL,  lim.trans = c(0, 1.5), prof=FALSE, k=2,  c.crit = 0.01, step=0.1)  
{
  cat("*** Checking for transformation for x ***", "\n") 
 ptrans<- function(x, p) if (abs(p)<=0.0001) log(x) else I(x^p)
  fn <- function(p) GAIC(gamlss(y~pb(pt(x,p)), c.crit = c.crit, trace=FALSE), k=k)
  if (prof) # profile dev
  {
    pp <- seq(lim.trans[1],lim.trans[2], step) 
    pdev <- rep(0, length(pp)) 
    for (i in 1:length(pp)) 
    {
      pdev[i] <- fn(pp[i])  
      #   cat(pp[i], pdev[i], "\n")
    }
    plot(pdev~pp, type="l")
    points(pdev~pp,col="blue")
    par <- pp[which.min(pdev)]
    cat('*** power parameters ', par,"***"," \n") 
  } else
  {
    fn <- function(p) GAIC(gamlss(y~pb(pt(x,p)), c.crit = c.crit, trace=FALSE), k=k)
    par <- optimise(fn, lower=lim.trans[1], upper=lim.trans[2])$minimum
    cat('*** power parameters ', par,"***"," \n") 
  }  
  par
}
# end of local function
#-------------------------------------------------------------------------------
## the families to fit
        FAM <- families      
## which method
  method.pb <- match.arg(method.pb)
## get the variables  
       ylab <- deparse(substitute(y))
       xlab <- deparse(substitute(x))
          y <- if (!is.null(data)) get(deparse(substitute(y)), envir=as.environment(data)) else y
          x <- if (!is.null(data)) get(deparse(substitute(x)), envir=as.environment(data)) else x
## ----------------------------------------------------------------------------- 
## if need to check for transformation in x    
  if (is.null(fix.power))
    {
    if (trans.x) # if x^p
    {
     ptrans<- function(x, p) if (p==0) log(x) else I(x^p)
      par <- findPower(y, x,   lim.trans = lim.trans, prof=prof, k=k,  c.crit = c.crit, step=0.1)    
      ox <- x
      x <-  pt(x,par)
    } 
    } else
    { par <- fix.power
      cat('*** power parameters fixed at ', par,"***"," \n")  
       ox <- x
        x <-  pt(x,par)
    }  
##  starting  model for fitted values for mu (we assuming that this will work).
##  Note no sigma is fitted here
##  fit the model -------------------------------------------------------------- 
    cat('*** Initial  fit***'," \n")       
    switch(method.pb, 
        "ML"= {m0 <- gamlss(y~pb(x), sigma.formula=~1, data=data, c.crit = 0.01)},
      "GAIC"= {m0 <- gamlss(y~pb(x, method="GAIC", k=k), sigma.formula=~1, data=data, c.crit = 0.01)}) ## initial fit  finish
## creating lists etc ----------------------------------------------------------
     failed <- list() 
       fits <- list()
        aic <- AIC(m0, k=k)
       fits <- c(fits, aic) 
  whichdist <- 0
## fitting the diferent models in FAM ------------------------------------------       
  for (i in 1:length(FAM)) 
  {
    cat('*** Fitting', FAM[i], "***","\n")  
     switch(method.pb, 
         "ML"= { m1 <- try(gamlss(y ~ pb(x, df=mu.df),
              sigma.formula = ~pb(x, df=sigma.df),
                 nu.formula = ~pb(x, df=nu.df), 
                tau.formula = ~pb(x, df=tau.df), 
                family = FAM[i], data = data,
               mu.start = fitted(m0), ...), silent=TRUE)},
        "GAIC"= { m1 <- try(gamlss(y ~ pb(x,  method="GAIC", k=k, df=mu.df),
              sigma.formula = ~pb(x,  method="GAIC", k=k, df=sigma.df),
                 nu.formula = ~pb(x,  method="GAIC", k=k, df=nu.df), 
                tau.formula = ~pb(x,  method="GAIC", k=k, df=tau.df), 
                family = FAM[i], data = data,
               mu.start = fitted(m0), ...), silent=TRUE)
         })      
    if (any(class(m1)%in%"try-error")) # if fitting failed
    {
      cat(FAM[i], " failed", "\n")
          failed <- c(failed, FAM[i]) 
    }
    else
    {
             aic <- AIC(m1, k=k)
      names(aic) <- FAM[i]
            fits <- c(fits, aic)
      if (AIC(m1, k=k) < AIC(m0, k=k)) 
      {
        m0 <-m1 
  whichdist <- i
      }
    }
  }
if(whichdist==0) 
{ # refitting the Normal with sigma  if not of the models is any good
  cat('*** Refitting', "NO", "***","\n")  
  m0 <-  switch(method.pb, 
                "ML"= {m0 <- gamlss(y~pb(x), sigma.formula=~pb(x), data=data, 
                                    c.crit = 0.01)},
                "GAIC"= {m0 <- gamlss(y~pb(x, method="GAIC", k=k), 
                              sigma.formula=~pb(x), 
                        data  =data, c.crit = 0.01)}) ## initial fit  finish
}
## changing the call t look better in the output -------------------------------
m0$call$mu.start <- NULL # this works OK
    m0$call$data <- substitute(data) # this is OK
  m0$call$family <- if(whichdist==0) "NO" else FAM[whichdist] # this is OK
# I am stack here
#m0$call$formula[3] 
#if (is.null(mu.df)) 
#  {m0$call$formula[[3]] <- sub("mu.df", "NULL", m0$call$formula[3])} else 
#  {m0$call$formula[[3]] <- sub("mu.df", mu.df, m0$call$formula[3])}
#
# FaM<-as.character(FAM[i])
# FaM
# m0$call$family <- substitute(Fam)
# m0$call
# 
# sub("FAM[i]", as.character(FAM[i]), m0$call,)
## transformation needed -------------------------------------------------------        
     if (trans.x||!is.null(fix.power))   
       { 
          x <- ox
   m0$power <- par 
       }
## save the rest information ---------------------------------------------------       
  m0$failed <- failed
       fits <- unlist(fits)
    m0$fits <- fits[order(fits)] 
    m0$xvar <- x#with(DaTa,x)
       m0$y <- y#with(DaTa,y)
    m0$ylab <- ylab
    m0$xlab <- xlab
  if (!is.null(data)) m0$call$data  <- substitute(data)
## calibration -----------------------------------------------------------------
  if (calibration)
  {
    calibration(m0, xvar=x, cent=cent, pch = 15, cex = 0.5, col = gray(0.7), ylab=ylab, xlab=xlab, legend=legend)	
  } 
  else 
  {
    centiles(m0, xvar=x, cent=cent, pch = 15, cex = 0.5, 
             col = gray(0.7), ylab=ylab, xlab=xlab, legend=legend)		
  }
 ## saving the fitted functions for mu sigma nu and tau  for prediction --------
if ("mu"%in%m0$par)
{
     muFun <- splinefun(x, fitted(m0,"mu"), method="natural")
  m0$muFun <- muFun
}
if ("sigma"%in%m0$par)
{
  sigmaFun <- splinefun(x, fitted(m0,"sigma"), method="natural")
  m0$sigmaFun <- sigmaFun
}
if ("nu"%in%m0$par)
{
  nuFun <- splinefun(x, fitted(m0,"nu"), method="natural")
  m0$nuFun <- nuFun
}
if ("tau"%in%m0$par)
{
  tauFun <- splinefun(x, fitted(m0,"tau"), method="natural")
  m0$tauFun <- tauFun
}
 #------------------------------------------------------------------------------
class(m0) <- c("lms", class(m0))
  m0  # save the last model
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# se ??
predict.lms <- function(object, 
                        what = c("all","mu", "sigma", "nu", "tau"),
                        parameter= NULL,
                     newdata = NULL,
                     ...
                     )
{
  what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
 # if newdata is not null check the class
 # if data.frame check for  object$xlab
 # other if vector use it 
  if (!is.null(newdata))
  {
  if (class(newdata)=="data.frame")
    {
    if (!object$xlab%in%names(newdata)) stop("the name in the data.frame do not much the x-variable in the model")
    x <- newdata[[object$xlab]]
    } else 
    {
     x <- newdata  
    }  
  }
  if (is.null(newdata))
  {
  out <-   switch(what,
          "all" = predictAll(object),
           "mu" = fitted(object, "mu"),
        "sigma" = fitted(object, "sigma"),
           "nu" = fitted(object, "nu"),
          "tau" = fitted(object, "tau"))
  } else
  {
    if(what=="all")
    {
      out <- list() 
      if ("mu"%in%object$par)
      {
        mu <- object$muFun(x)
      out$mu <- mu  
      }
      if ("sigma"%in%object$par)
      {
        sigma <- object$sigmaFun(x)
    out$sigma <- sigma  
      }
      if ("nu"%in%object$par)
      {
        nu <- object$nuFun(x)
    out$nu <- nu  
      }
      if ("tau"%in%object$par)
      {
        tau <- object$tauFun(x)
    out$tau <- tau  
      }
    }
    if(what=="mu")      out <- object$muFun(x)
    if(what=="sigma")   out <- object$sigmaFun(x)
    if(what=="nu")      out <- object$nuFun(x)
    if(what=="tau")     out <- object$tauFun(x)
  }  
out    
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# this function is appropriate to used when fitted model fails to c
calibration <- function(object, xvar, cent=100*pnorm((-4:4)*2/3), legend=FALSE, fan=FALSE,  ...)
{
 z   <-  quantile(resid(object), probs = cent/100)
 p   <-  pNO(z, mu=0, sigma=1)
 percent <- 100*p
 if (fan)
 {
   centiles.fan(object, xvar=xvar, cent=percent,   ...)  
 }
 else
 {
   centiles(object, xvar=xvar, cent=percent, legend=legend,  ...)
 }
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#------------------------------------------------------------------------------- 
# findPower <- function(y, x, data = NULL,  lim.trans = c(0, 1.5), prof=FALSE, k=2,  c.crit = 0.01, step=0.1)  
# {
#    ylab <- deparse(substitute(y))
#    xlab <- deparse(substitute(x))
#       y <- if (!is.null(data)) get(deparse(substitute(y)), envir=as.environment(data)) else y
#       x <- if (!is.null(data)) get(deparse(substitute(x)), envir=as.environment(data)) else x
#   ##  checking for transformation in x        
#     cat("*** Checking for transformation for x ***", "\n") 
#     ptrans<- function(x, p) if (p==0) log(x) else I(x^p)
#      fn <- function(p) GAIC(gamlss(y~pb(pt(x,p)),data=data, c.crit = c.crit, trace=FALSE), k=k)
#  if (prof) # profile dev
#  {
#    pp <- seq(lim.trans[1],lim.trans[2], step) 
#  pdev <- rep(0, length(pp)) 
#    for (i in 1:length(pp)) 
#      {
#      pdev[i] <- fn(pp[i])  
#   #   cat(pp[i], pdev[i], "\n")
#      }
#    plot(pdev~pp, type="l")
#    points(pdev~pp,col="blue")
#  par <- pp[which.min(pdev)]
#  cat('*** power parameters ', par,"***"," \n") 
#  } else
#  {
# #     fn <- function(p) GAIC(gamlss(y~pb(pt(x,p)),sigma.fo=~pb(pt(x,p)),nu.fo=~pb(pt(x,p)), data=data,tau.fo=~pb(pt(x,p)), c.crit = c.crit, trace=FALSE, family=BCT), k=k)
# #    fn <- function(p) GAIC(gamlss(y~pb(pt(x,p)),sigma.fo=~pb(pt(x,p)), data=data, c.crit = c.crit, trace=FALSE), k=k)
#        fn <- function(p) GAIC(gamlss(y~pb(pt(x,p)), data=data, c.crit = c.crit, trace=FALSE), k=k)
#    par <- optimise(fn, lower=lim.trans[1], upper=lim.trans[2])$minimum
#   # browser()
#  #  par <- optim(.5, fn, lower=lim.trans[1], upper=lim.trans[2], method="L-BFGS-B")$par
#    cat('*** power parameters ', par,"***"," \n") 
#  }  
#    par
#   }  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
z.scores <- function(object, y,x)
{
  if (!is(object,"lms"))  stop(paste("This is not an lms object", "\n", "")) 
  if (is.null(y)) stop("the y values should be set for z-scores")
  if (is.null(x)) stop("the x values should be set for z-scores")
  if (length(y)!= length(x)) stop("length of x and y is not the same")
     pred <- predict(object, newdata=x)
    fname <- object$family[1]
    lpar <- length(object$par)
     qfun <- paste("p",fname,sep="")
  if(lpar==1) 
   {newcall <-call(qfun,y,mu=pred$mu) }
  else if(lpar==2)
  {newcall <-call(qfun,y, mu=pred$mu,sigma=pred$sigma) }
  else if(lpar==3)
  {newcall <-call(qfun,y,mu=pred$mu,sigma=pred$sigma,nu=pred$nu) }
  else 
  {newcall <-call(qfun,y,mu=pred$mu,sigma=pred$sigma,nu=pred$nu,tau=pred$tau) }
  cdf <- eval(newcall)       
  rqres <- qnorm(cdf)
  rqres  
}

  
 
