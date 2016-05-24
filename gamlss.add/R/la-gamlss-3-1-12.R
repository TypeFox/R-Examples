# ms Tuesday, July 7, 2009 
# MS Thusday 3-1-12
# it uses penlags()
# which is used in the backfitting 
#----------------------------------------------------------------------------------------
la <-function(x, control=la.control(...), ...)
{ # how does it knows if x or y id fitted 
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
       LambdaName <- paste("LaMbdA",fourLetters, sep=".")
## put the starting values in the gamlss()environment
#--------
   assign(LambdaName, control$start.lambda, envir=gamlss.env)
      len <- length(x) # get the lenth of the data
## out
             xvar <- rep(0,  len) # model.matrix(as.formula(paste(paste("~",ff, sep=""), "-1", sep="")), data=Data) #
  # attr(xvar,"formula")     <- formula
      # xvar <- x  #rep(0,length(x)) # only the linear part in the design matrix the rest pass as artributes
      attr(xvar, "values")        <- x 
      attr(xvar, "control")       <- control 
      attr(xvar, "call")          <- substitute(gamlss.la(data[[scall]], z, w)) 
      attr(xvar, "gamlss.env")    <- gamlss.env
      attr(xvar, "NameForLambda") <- LambdaName
      attr(xvar, "class")         <- "smooth"
      xvar
}
##---------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
la.control <- function( lags = 10, from.lag=0, df = NULL, lambda = NULL,  start.lambda = 10, order = 1, plot = FALSE, method = c("ML", "GAIC"),  k = 2,  ...)
{ 
        if(lags<= 0) {
warning("the value of lags supplied is less than 0, the value of 10 was used instead")
                inter <- 10 }         
        if(order < 0) {
warning("the value of order supplied is zero or negative the default value of 2 was used instead")
                order <- 2}
        if(k <= 0) {
warning("the value of GAIC/GCV penalty supplied is less than zero the default value of 2 was used instead")
                k <- 2}   
method <- match.arg(method)                          
        list(lags=lags, from.lag=from.lag, df = df, lambda = lambda, order = order, start.lambda=start.lambda, plot=plot, method= method, k=k )
}
##---------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# the definition of the backfitting additive function
gamlss.la <-function(x, y, w, xeval = NULL, ...)
{                  
           xvar <-  if (is.null(xeval)) attr(x,"values")
                   else  attr(x,"values")[seq(1,length(y))]
     gamlss.env <- as.environment(attr(x, "gamlss.env"))
     LambdaName <- as.character(attr(x, "NameForLambda")) 
    lambdastart <- get(LambdaName, envir=gamlss.env) ## geting the starting knots 
        control <- as.list(attr(x, "control"))  
       # cat("df", df, "lambda"=lamnda, )
      # browser()
            fit <- penLags(y=y, x=xvar, lags=control$lags, from.lag=control$from.lag, weights=w, df=control$df, lambda=control$lambda, start.lambda=lambdastart,   order=control$order, plot=control$plot, method=control$method, k=control$k)     
          resid <- y-fitted(fit)
           #cat("knot", knots(fit), "\n")
        assign(LambdaName, fit$lambda, envir=gamlss.env)
       
  if (is.null(xeval))
    {
    list(fitted.values=fitted(fit), residuals=resid,  nl.df = fit$mu.df-1, 
         lambda=fit$lambda, coefSmo = fit, var=NA)
   }    
  else 
   {
      lenN <- length(attr(x,"x"))
     nxval <- attr(x,"x")[(length(y)+1):lenN]
      pred <- predict(fit,newdata=nxval)
      pred  
   }         

}
      
