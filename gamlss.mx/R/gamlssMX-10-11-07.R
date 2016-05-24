#-------------------------------------------------------------------------------
## this  a general function for fitting mixures within gamlss
## created by Mikis Stasinopoulos and Bob Rigby 
## LAST CHECKED 5-8-14 
## latest change Tuesday, April 10, 2007 at 08:54
## Some of the models which can be fitted 
## with the function gamlssNP can be identical 
## it needs a summary functions  
## In this version the prior probabilities
## pi can be modelled 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
gamlssMX  <- function ( formula = formula(data), 
                     pi.formula = ~1, 
                         family = "NO", # note is character
                         weights, 
                              K = 2, 
                           prob = NULL,
                           data = sys.parent(),
                        control = MX.control(), 
                      g.control = gamlss.control(trace=FALSE),
                 zero.component = FALSE,     
                      ...)
 {
#-------------------------------------------------------------------------------
.gamlss.bi.list <- eval(quote(.gamlss.bi.list), envir = getNamespace("gamlss"))
# extra functions within -------------------------------------------------------
# for getting the commulative function
get.the.p.function <- function(object, ...)
{
#gamlss.bi.list <- c("BI", "Binomial", "BB", "Beta Binomial")
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
    {pfun <- if (fname%in%.gamlss.bi.list)  eval(call(dfun,q=object$y, bd=object$bd,  mu=fitted(object,"mu"),sigma=fitted(object,"sigma"), nu=fitted(object,"nu")))
             else eval(call(dfun,q=object$y, mu=fitted(object,"mu"), sigma=fitted(object,"sigma"), nu=fitted(object,"nu")))},
    {pfun <- if (fname%in%.gamlss.bi.list)  eval(call(dfun,q=object$y, bd=object$bd,  mu=fitted(object,"mu"),sigma=fitted(object,"sigma"), nu=fitted(object,"nu"), tau=fitted(object,"tau")))
             else eval(call(dfun,q=object$y, mu=fitted(object,"mu"), sigma=fitted(object,"sigma"), nu=fitted(object,"nu"), tau=fitted(object,"tau")))})
pfun 
}
#-------------------------------------------------------------------------------
get.likelihood <- function(obj)
{
#gamlss.bi.list <- c("BI", "Binomial", "BB", "Beta Binomial")
if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
   fname <- obj$family[[1]]
 DistPar <- obj$parameters
   nopar <- length(DistPar) 
    dfun <- paste("d",fname,sep="")
 switch(nopar,  
    {lik <- if (fname%in%.gamlss.bi.list)  eval(call(dfun,x=obj$y, bd=obj$bd,  mu=fitted(obj,"mu")))
             else   eval(call(dfun, x=obj$y, mu=fitted(obj)))},
    {lik <- if (fname%in%.gamlss.bi.list)  eval(call(dfun,x=obj$y, bd=obj$bd,  mu=fitted(obj,"mu"), sigma=fitted(obj,"sigma") ))
             else   eval(call(dfun,x=obj$y, mu=fitted(obj), sigma=fitted(obj,"sigma"))) },
    {lik <- if (fname%in%.gamlss.bi.list)  eval(call(dfun,x=obj$y, bd=obj$bd,  mu=fitted(obj,"mu"), sigma=fitted(obj,"sigma"), nu=fitted(obj,"nu") ))
            else eval(call(dfun,x=obj$y, mu=fitted(obj), sigma=fitted(obj,"sigma"), nu=fitted(obj,"nu")))},
    {lik <- if (fname%in%.gamlss.bi.list)  eval(call(dfun,x=obj$y, bd=obj$bd,  mu=fitted(obj,"mu"), sigma=fitted(obj,"sigma"), nu=fitted(obj,"nu"), tau=fitted(obj,"tau") ))
            else eval(call(dfun,x=obj$y, mu=fitted(obj), sigma=fitted(obj,"sigma"), nu=fitted(obj,"nu"), tau=fitted(obj,"tau")))})
lik
}           
#-------------------------------------------------------------------------------
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
      for ( i in seq(1,n))
      {
          xx <- rbind(xx,matrix(rep(x[i,],ni[i]),ncol=c,byrow=TRUE))
      }
      xx
  }
}
#-------------------------------------------------------------------------------
# the proper function starts here
#-------------------------------------------------------------------------------
#library(gamlss)
#library(nnet)
 gamlssMXcall <- match.call()  #   the function call
## checking for NA in the data 
 if(!missing(data) & any(is.na(data)))   
      stop("The data contains NA's, use data = na.omit(mydata)") 
## prob are the starting values of the pi's   
 if(zero.component == TRUE)
 { KK <- K+1 
  if (is.null(prob)) 
   { 
     prob <- rep(1/KK, KK)
   }
  else
   {
    if (length(prob)!=KK) stop("the length of prob must be  K+1(for zero component) ")
   } 
 }
 else
 {
   if (is.null(prob)) 
     { 
     prob <- rep(1/K, K)
     }
 else
   {
    if (length(prob)!=K) stop("the length of prob must be K")
   } 
 }
## expand the formula if not list
 if (!is.list(formula))
  {
   allFormula <- vector("list",K)
   for (i in 1:K) allFormula[[i]] <- formula
  }
 else allFormula <- formula
## expand the family if not list
 if (!is.list(family))
  {
   allFamily <- vector("list",K)
   for (i in 1:K) allFamily[[i]] <- family
  }
  else allFamily <- family
## check that the length for both the family and formula is K
if (length(allFamily)!=K) stop("the length of the list for the family is not equal the number of componets in the mixture")
if (length(allFormula)!=K) stop("the length of the list for the formula is not equal the number of componets in the mixture")
## get response for his length
        Y <- model.extract(model.frame(allFormula[[1]],data=data),"response")
        N <-  if(is.null(dim(Y))) length(Y) else dim(Y)[1]# calculate the dimension for y  
pweights. <- pweights. <<- if (missing(weights)) rep(1,N) else weights  
## get things from control
  set.seed(control$seed) # set the seed 
prob.sample <- if (is.null(control$sample)) 10/N else control$sample  
## creating the matrix of posterior probabilities (weights)
        W <- matrix(1, nrow=N, ncol=K ) 
allModels <- vector("list",K)
## get starting values for the weights (posterior probabilities)
 for (i in 1:K)
   {         wSam <- sample(c(0,1), N, replace = TRUE, prob=c(1-min(.5,prob.sample),min(.5,prob.sample)))  
              ww. <- ww. <<- wSam # W[,i] # put in 
             #data1 <<- data.frame(data, ww, pweights)
    allModels[[i]] <- gamlss(allFormula[[i]], weights=ww.*pweights., data = data, family = allFamily[[i]], control=g.control, ...)
             W[,i] <- get.likelihood(allModels[[i]]) 
   } 
 if(zero.component == TRUE) W[,KK] <- prob[KK]*ifelse(Y==0,1,0)
      SumLik <- rowSums(W)
           W <- W/SumLik
           W <- ifelse(W<1e-20,0,W)
        prob <- colSums(W)/N
 # -2 * logLikelihood
       newdv <- -2*sum(log(SumLik))    # the global deviance
       olddv <- newdv + 1  
    iter.num <- 1 
 trace.print <- 0 
    dev.fits <- rep(0,control$n.cyc)
#-------------------------------------------------------------------------------
# if require model for prior propabilities pi's 
     modelPi <- if (pi.formula[[2]]==1) FALSE else TRUE
  if (modelPi)
   {
    if (missing(data)) stop("the data argument is needed if pi's are modelled")
     if(zero.component == TRUE) 
      {
      dataKK <- expand.vc(data,KK)# expand data.
     fac.fit <- gl(KK,N)
      }
      else
      { 
      dataKK. <<- expand.vc(data,K)# expand data
      fac.fit <- gl(K,N) # get the factor
      }
     res.var <- model.matrix(~fac.fit-1) # get the matrix of zeroes ans ones
   form.prob <-  update(pi.formula, res.var~.)# formula(paste("res.var~", pi.formula[[2]])) # get the formula 
        PROB <-  matrix(prob, ncol=K, nrow=N, byrow=TRUE) # expand prior prob
   }
## EM starts here--------------------------------------------------------------- 
 while ( abs(olddv-newdv) > control$cc && iter.num < control$n.cyc ) # MS Wednesday, June 26, 2002 
  {
  for (i in 1:K)
    {          ww. <- ww. <<- W[,i]
    allModels[[i]] <- gamlss(allFormula[[i]], weights=ww.*pweights., data = data, family = allFamily[[i]], control=g.control,...)
             W[,i] <- if (modelPi) PROB[,i]*get.likelihood(allModels[[i]]) else prob[i]*get.likelihood(allModels[[i]])
    } 
     if(zero.component == TRUE)
      {
      W[,KK] <-  if (modelPi) PROB[,KK]*ifelse(Y==0,1,0) else prob[KK]*ifelse(Y==0,1,0)  
      }
            SumLik <-rowSums(W)
                 W <- W/SumLik 
                 W <- ifelse(W<1e-20,0,W) 
                di <- -2*log(SumLik) 
             olddv <- newdv
             newdv <- sum(di)
       if (modelPi)
          {
          dataKK.$wWw <- wWw <- as.vector(W)
          dataKK.$res.var <- res.var
          mult.mod <- multinom(form.prob, weights=wWw, data=dataKK., trace=FALSE)
           fitted(mult.mod)[1:11,]
              PROB <- fitted(mult.mod)[1:N,] # take only N of them 
          }
       else { prob <- colSums(W)/N}
   if (control$trace == TRUE) cat("GAMLSS-MX iteration ",iter.num, "Global deviance =", newdv, "\n"  )
   else if (is.numeric(control$trace))
        { 
          trace.print <- trace.print+1
         if (trace.print==control$trace)
          {cat("GAMLSS-MX iteration ",iter.num, "Global deviance =", newdv, "\n"  )
          trace.print <- 0 
          }
        }
dev.fits[iter.num] <- newdv 
          iter.num <- iter.num+1    
  }
if (control$plot==TRUE) plot(dev.fits[dev.fits!=0], type="l", xlab="EM iterations", ylab="-2 logLik")
## EM finish here--------------------------------------------------------------- 
## the output starts here 
## residuals -------------------------------------------------------------------
                WF <- matrix(0, ncol=K, nrow=N)
   if (modelPi)
    {
     for (i in 1:K)     WF[,i] <- PROB[,i]*get.the.p.function(allModels[[i]]) 
    }
    else 
    {
    for (i in 1:K)     WF[,i] <- prob[i]*get.the.p.function(allModels[[i]]) 
    }
               res <-qnorm(rowSums(WF))

#------------------------------------------------------------------------------- 
## degrees of freedom
df.fit <- 0 
for (i in 1:K) df.fit <- df.fit + allModels[[i]]$df.fit
df.fit <- if(modelPi) df.fit+mult.mod$edf  else df.fit+K-1 # 
if (modelPi)  colnames(PROB) <- paste("pi",seq(1,K), sep="")
familyAll <- rep("NULL",K)
for (i in 1:K) familyAll[i] <- allModels[[i]]$family[1]
## all output
out <- list()
out <- list(models = allModels,
          model.pi = if (modelPi) mult.mod else NULL,
        G.deviance = newdv, 
            df.fit = df.fit, 
       df.residual = N-df.fit,
         post.prob = W, 
              prob = if (modelPi) PROB else prob,
            family = familyAll,
              call = gamlssMXcall,
               aic = newdv+2*df.fit,
               sbc = newdv+log(N)*df.fit,
                 K = K,
                 N = N,
           weights = pweights., 
         residuals = res,
              seed = control$seed
                 )
class(out) <- list("gamlssMX", "gamlss")
if (modelPi) on.exit(rm(ww.,pweights.,dataKK., envir=sys.frame(0)))
else on.exit(rm(ww.,pweights., envir=sys.frame(0)))
out
 }  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
MX.control <- function(cc=0.0001, n.cyc=200, trace=FALSE, seed=NULL, plot=TRUE, sample=NULL, ...)
{
  ##list(cc=0.0001, n.cyc=100, trace=TRUE),
  ##  Control iteration for GAMLSS
  ##  Mikis Stasinopoulos Monday, March 25, 2002 at 16:17
  ##
  if(cc <= 0) {
    warning("the value of cc supplied is zero or negative the default value of 0.0001 was used instead")
    c.crit <- 0.0001}
  if(n.cyc < 1) {
    warning("the value of no cycles supplied is zero or negative the default value of 20 was used instead")
    n.cyc <- 100}
  if(is.logical(trace)) trace <- trace 
  else if (is.numeric(trace) & trace <= 0) 
  {warning("the value of trace supplied is less or equal t zero the default of 1 was used instead")
   trace <- 1
  }
  seed <- if (is.null(seed) ) floor(runif(min=0, max=100000, n=1)) else seed 
  if (!is.null(sample)) 
  { if (sample>1||sample<0)
  {
    sample <-0.1 
    warning("the value of sample supplied is not probability the default of 0.10 was used instead")
  }
  else sample <- sample  
  } 
  plot <- if(is.logical(plot)) plot else TRUE     
  list(cc = cc, n.cyc = n.cyc, trace = trace, seed=seed, plot=plot, sample=sample)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# METHODS for "gamlssMX"
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
print.gamlssMX <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    K <- x$K 
    cat("\nMixing Family: ", deparse(x$family), "\n", fill=TRUE)
    cat("Fitting method: EM algorithm \n")
    cat("\nCall: ", deparse(x$call),  "\n", fill=TRUE)
     for (i in 1:K)
       {
        if ("mu" %in% x$models[[i]]$parameters) 
         {      
          cat("Mu Coefficients for model:", i,  "\n")   
          print.default(format(coef(x$models[[i]], "mu"), digits = digits), print.gap = 2, quote = FALSE)        
         }
        if ("sigma" %in% x$models[[i]]$parameters) 
         { 
          cat("Sigma Coefficients for model:", i,  "\n")           
          print.default(format(coef(x$models[[i]], "sigma"), digits = digits), print.gap = 2, quote = FALSE)        
         }
        if ("nu" %in% x$models[[i]]$parameters) 
         {       
         cat("Nu Coefficients for model:", i,  "\n")    
         print.default(format(coef(x$models[[i]], "nu"), digits = digits), print.gap = 2, quote = FALSE)        
         }
        if ("tau" %in% x$models[[i]]$parameters) 
         {
          cat("Tau Coefficients for model:", i,  "\n")          
          print.default(format(coef(x$models[[i]], "tau"), digits = digits), print.gap = 2, quote = FALSE)        
         }
       }
       if(is.matrix(x$prob)) 
         {
          cat("model for pi: \n")
          print(coef(x$model.pi))
          cat("\nEstimated probabilities: \n")
          print(x$prob[1:3,])
          cat("... \n" )
         }
       else   cat("\nEstimated probabilities:", x$prob, "\n" )
    cat("\nDegrees of Freedom for the fit:", x$df.fit, "Residual Deg. of Freedom  ", 
        x$df.residual, "\n")
   cat("Global Deviance:    ", format(signif(x$G.deviance)), 
        "\n            AIC:    ", format(signif(x$aic)), "\n            SBC:    ", 
        format(signif(x$sbc)), "\n")
    invisible(x)
}

################################################################################
#                         fitted.gamlssMX
################################################################################
# this gives the compoment fv if K is set otherwise average the componets 
# using the fitted probabilities
# MS is this justified if mu not the mean? 
#  also does it makes sense for other parameters??
fitted.gamlssMX<-function (object, K=1, ... ) 
{
if (K%in%seq(1:object$K))
 {
  x <- fitted(object$models[[K]], ...)
 }
else 
  {
  WF <- matrix(0, ncol=object$K, nrow=object$N)
  for (i in 1:object$K)     WF[,i] <- object$prob[i]*fitted(object$models[[i]]) 
  x <-rowSums(WF)
  }
 x
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# coefficints of the componets
coef.gamlssMX <- function(object, K=1, ...)
{
coef(object$models[[K]], ...)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
formula.gamlssMX <- function(x, K=1, ...)
{
formula(x$models[[K]], ...)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
model.matrix.gamlssMX <- function(object, K=1, ...)
{
model.matrix(object$models[[K]], ...)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
terms.gamlssMX <- function(x, K=1, ...)
{
terms(x$models[[K]], ...)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
predict.gamlssMX <- function(object, K=1, ...)
{
predict(object$models[[K]], ...)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
residuals.gamlssMX <- function(object,...)
{
object$residuals
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# this function fits n models with different starting values 
# and select the one with smallest deviance 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# this function fits n models with different starting values 
# and select the one with smallest deviance 
# last modification Thursday, April 30, 2009 
gamlssMXfits <- function(    n = 5,
                       formula = formula(data), 
                    pi.formula = ~1, 
                        family = "NO", # note is character
                        weights, 
                             K = 2, 
                          prob = NULL,
                          data = sys.parent(),
                       control = MX.control(), 
                     g.control = gamlss.control(trace=FALSE),
                     zero.component = FALSE,   
                     ...)
{
 gamlssMXfitscall <- match.call()  
        modellist <- list()
              dev <- rep(0,length=n) 
             seed <- floor(runif(min=0,max=10000, n=n))
for (i in (1:n))
 {
  m1 <- try(gamlssMX(formula=formula,  pi.formula = pi.formula, family = family ,K=K,  prob=prob, data = 
                    data, control=MX.control(seed=seed[i]), g.control = gamlss.control(trace=FALSE), 
                    zero.component=zero.component)) #
             if (any(class(m1)%in%"try-error")) 
               { 
                cat("model=", i, "failed", "\n") 
                dev[i] <- NA    
        modellist[[i]] <- NA
                 next
               }
 modellist[[i]] <- m1 
         dev[i] <- deviance(modellist[[i]])
 cat("model=", i,"\n")                    
 }
                   II <- which.min(dev)
                model <- modellist[[II]]
   gamlssMXfitscall$n <- NULL
gamlssMXfitscall[[1]] <- as.name("gamlssMX")
           model$call <- gamlssMXfitscall
          model$extra <- list(dev=dev, seed=seed, which=II)
 model
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
update.gamlssMX <- function (object, 
                          formula., 
                          ..., 
                          what = c("mu", "sigma", "nu", "tau"),
                          evaluate = TRUE) 
{
    call <- object$call
    if (is.null(call)) 
        stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.)) 
        {
        what <- match.arg(what) 
        if (what=="mu") 
         { call$formula <- update.formula(formula(object,what), formula.) }
        else  
         {
          call[[paste(what,"formula",sep=".")]] <- 
         if (length(update.formula(formula(object,what), formula.))==2)
          update.formula(formula(object,what), formula.)
         else
          update.formula(formula(object,what), formula.)[-2]
         }
        }
    if (length(extras) > 0) 
        {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) 
           {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
           }
        }
    if (evaluate) 
        eval(call, parent.frame())
    else call
}
################################################################################
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# this is a function to calculate the pdf of y fY(y)=p1*f1(y)+p2*f(y), of a MX model
# MS + BR 
dMX <- function(y, 
                 mu = list(mu1=1,mu2=5), 
              sigma = list(sigma1=1,sigma2=1), 
                 nu = list(nu1=1,nu2=1), 
                tau = list(tau1=1,tau2=1), 
                 pi = list(pi1=.2,pi2=.8),
                # K = 2, 
             family = list(fam1="NO", fam2="NO"),  
                log = FALSE,
                 ...)
  {
#  gamlss.bi.list <- c("BI", "Binomial", "BB", "Beta Binomial")
#  .gamlss.bi.list<-c("BI", "Binomial", "BB", "Beta Binomial")
 ## checking the probabilities 
     sump <- 0
        K <- length(pi) ## set the length using the length of propabilities
  for (i in 1:K) sump <-pi[[i]]+ sump 
  if (any(sump!=1)) stop(paste("the vector pi should sum to 1"))
 ## get the length of y and the length of the parameters
 if(is.null(dim(y))) N <- length(y) else N <- dim(y)[1]  
       Prob <- matrix(1, nrow=N, ncol=K ) # rep(prob,rep(N,K)) 
     for (i in 1:K) # go throught the components
   { 
        fam <-  as.gamlss.family(family[[i]])
      fname <- fam$family[[1]]
    DistPar <- fam$parameters
      nopar <- length(DistPar)  
       dfun <- paste("d",fname,sep="") 
       if (fname%in%.gamlss.bi.list) # for binomial data
    {
     if (NCOL(y) == 1)
        {
        y <- y 
       bd <- rep(1, length(y))
        }
     else if (NCOL(y) == 2)
        {
        y <- y[,1]
       bd <- y[,1]+y[,2]
        }
     else
        {
      stop("wrong response variable for binomial data")
        }
     }           
   switch(nopar,  
      {lik <- if (fname%in%.gamlss.bi.list){eval(call(dfun, x=y, bd=bd,  mu=mu[[i]]))}
            else  eval(call(dfun, x=y, mu=mu[[i]]))
      },
      {lik <- if (fname%in%.gamlss.bi.list) eval(call(dfun, x=y, bd=bd,  mu=mu[[i]], sigma=sigma[[i]]))
            else  eval(call(dfun, x=y, mu=mu[[i]], sigma=sigma[[i]])) 
      },
      {lik <- if (fname%in%.gamlss.bi.list) eval(call(dfun, x=y, bd=bd,  mu=mu[[i]], sigma=sigma[[i]], nu=nu[[i]])) 
            else eval(call(dfun,x=y, mu=mu[[i]], sigma=sigma[[i]], nu=nu[[i]]))},
      {lik <-  if (fname%in%.gamlss.bi.list) eval(call(dfun, x=y, bd=bd,  mu=mu[[i]], sigma=sigma[[i]], nu=nu[[i]], tau=tau[[i]])) 
            else eval(call(dfun,x=y, mu=mu[[i]], sigma=sigma[[i]], nu=nu[[i]], tau=tau[[i]]))})
    Prob[,i]<-pi[[i]]*lik
    }  # finish go throught the components
       fy <- rowSums(Prob)
       fy <- if(log == FALSE) fy else log(fy)
       fy 
  }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# this is a function to calculate the cdf y from a MX model
pMX <- function(q, 
                 mu = list(mu1=1,mu2=5), 
              sigma = list(sigma1=1,sigma2=1), 
                 nu = list(nu1=1,nu2=1), 
                tau = list(tau1=1,tau2=1), 
                  pi = list(pi1=.2,pi2=.8),
                #  K = 2, 
             family = list(fam1="NO", fam2="NO"),  
                 log = FALSE,
                 ...)
  {
 # gamlss.bi.list<-c("BI", "Binomial", "BB", "Beta Binomial")
## checking the probabilities 
       sump <- 0
          K <- length(pi) 
  for (i in 1:K) sump <-pi[[i]]+ sump 
  if (any(sump!=1)) stop(paste("the vector pi should sum to 1"))
 ## get the length of q and the length of the parameters
 if(is.null(dim(q))) N <- length(q) else N <- dim(q)[1]  
       Prob <- matrix(1, nrow=N, ncol=K ) # rep(prob,rep(N,K))  
     for (i in 1:K)
   { 
        fam <- as.gamlss.family(family[[i]])
      fname <- fam$family[[1]]
    DistPar <- fam$parameters
      nopar <- length(DistPar)  
       dfun <- paste("p",fname,sep="")
        if (fname%in%.gamlss.bi.list) # for binomial data
    {
     if (NCOL(q) == 1)
        {
        q <- q 
       bd <- rep(1, length(q))
        }
     else if (NCOL(q) == 2)
        {
        q <- q[,1]
       bd <- q[,1]+q[,2]
        }
     else
        {
      stop("wrong response variable for binomial data")
        }
     }            
   switch(nopar,  
     {lik <- if (fname%in%.gamlss.bi.list) eval(call(dfun, q=q, bd=bd[[i]],  mu=mu[[i]]))
            else   eval(call(dfun, q=q, mu=mu[[i]]))
     },
     {lik <- if (fname%in%.gamlss.bi.list) eval(call(dfun, q=q, bd=bd[[i]],  mu=mu[[i]], sigma=sigma[[i]]))
            else   eval(call(dfun, q=q, mu=mu[[i]], sigma=sigma[[i]])) 
     },
     {lik <- if (fname%in%.gamlss.bi.list) eval(call(dfun, q=q, bd=bd[[i]],  mu=mu[[i]], sigma=sigma[[i]],  nu=nu[[i]]))
            else eval(call(dfun,q=q, mu=mu[[i]], sigma=sigma[[i]], nu=nu[[i]]))},
     {lik <- if (fname%in%.gamlss.bi.list) eval(call(dfun, q=q, bd=bd[[i]],  mu=mu[[i]], sigma=sigma[[i]],  nu=nu[[i]], tau=tau[[i]]))
            else eval(call(dfun,q=q, mu=mu[[i]], sigma=sigma[[i]], nu=nu[[i]], tau=tau[[i]]))})
    Prob[,i]<-pi[[i]]*lik
    } 
            fy <-rowSums(Prob)
       fy <- if(log == FALSE) fy else log(fy)
       fy 
  }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------