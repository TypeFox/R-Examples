# ms Tuesday, July 7, 2009 
# fit smoothing terms using the  curfit.free.knot() function
# which is used in the backfitting 
# TO DO:
#----------------------------------------------------------------------------------------
fk <-function(x, start=NULL, control=fk.control(...), ...) 
{ 
#------------------------------------------
# function starts here
#------------------------------------------
    scall <- deparse(sys.call(), width.cutoff = 500L)
# get where "gamlss" is in system call
# it can be in gamlss() or predict.gamlss()       
    rexpr <- grepl("gamlss",sys.calls()) ## 
for (i in length(rexpr):1)
   { 
 position <- i # get the position, we are geting the fist from the last
 if (rexpr[i]==TRUE) break
   }
gamlss.env <- sys.frame(position) #gamlss or predict.gamlss
#--------
## get a random name to use it in the gamlss() environment
#--------
               sl <- sample(letters, 4)
      fourLetters <- paste(paste(paste(sl[1], sl[2], sep=""), sl[3], sep=""),sl[4], sep="")
  startLambdaName <- paste("start.Lambda",fourLetters, sep=".")
## put the starting values in the gamlss()environment
#--------
   assign(startLambdaName, list(start=start,iter=1), envir=gamlss.env)
      len <- length(x) # get the lenth of the data
## out
     xvar <- rep(0,  len) #
   attr(xvar, "x")             <- x
   attr(xvar,"control")        <- control
   attr(xvar, "gamlss.env")    <- gamlss.env
   attr(xvar, "NameForLambda") <- startLambdaName
   attr(xvar, "call")          <- substitute(gamlss.fk(data[[scall]], z, w, ...)) 
   attr(xvar, "class")         <- "smooth"
   xvar
}
##---------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
fk.control <-  function ( degree=1, all.fixed=FALSE, fixed = NULL, base=c("trun","Bbase"))
{
  list(knots=knots,  degree=1, all.fixed=all.fixed, fixed = NULL, base=c("trun","Bbase"))
}                      
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# the definition of the backfitting additive function
gamlss.fk <-function(x, y, w, xeval = NULL, ...)
{              
      xvar <-  if (is.null(xeval)) attr(x,"x")
             else  attr(x,"x")[seq(1,length(y))]
        control <- as.list(attr(x, "control")) 
         degree <- control$degree
     gamlss.env <- as.environment(attr(x, "gamlss.env"))
startLambdaName <- as.character(attr(x, "NameForLambda"))  
         lambda <- get(startLambdaName, envir=gamlss.env)$start
    ifFirstIter <- get(startLambdaName, envir=gamlss.env)$iter
      ## geting the starting knots 
   #   cat(lambda)
     if (control$all.fixed==TRUE||ifFirstIter)
     {
      fit <- fitFixedKnots(y=y, x=xvar,  weights=w, degree=degree, knots = lambda, fixed=control$fixed,  base=control$base)
     }
     else
     {
       fit <- fitFreeKnots(y=y, x=xvar,  weights=w, degree=degree, knots = lambda, fixed=control$fixed, base=control$base)     
     }
     #     browser()      
    #   cat("knot", knots(fit), "\n")
    #  plot(y~xvar)
    #  lines(fitted(fit)~xvar, col="red")
        assign(startLambdaName, list(start=fit$breakPoints,iter=0), envir=gamlss.env)
  if (is.null(xeval))
   {
   list(fitted.values=fitted(fit), residuals=y-fitted(fit),  nl.df = fit$df-1,# -1 if linear is not in
      lambda=knots(fit), ## we nead df's here 
     coefSmo = fit, var=NA)
   }    
else 
   {
      lenN <- length(attr(x,"x"))
     nxval <- attr(x,"x")[(length(y)+1):lenN]
      pred <- predict(fit,newdata=nxval)
      pred  
   }         

}
      
