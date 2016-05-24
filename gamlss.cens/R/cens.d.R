cens.d <- function(family = "NO", type = c( "right", "left", "interval"),...)
  {
    #if (!is.Surv(y)) stop(paste("the y variable is not a Surv object"))
     type <- match.arg(type)
    if (type=="counting") stop(paste("gamlss has no support for counting"))
   fname <- family
    if (mode(family) != "character" && mode(family) != "name")
   fname <- as.character(substitute(family))
 distype <- eval(call(family))$type
    dfun <- paste("d",fname,sep="")
    pfun <- paste("p",fname,sep="")
     pdf <- eval(parse(text=dfun))
     cdf <- eval(parse(text=pfun))
fun <- if (type=="left")  
       function(x, log = FALSE, ...)
        {
         if (!is.Surv(x)) stop(paste("the x variable is not a Surv object"))
        dfun <- ifelse(x[,"status"]==1, pdf(x[,1],log = TRUE,...),log(cdf(x[,1],...)))
        dfun <- if (log == TRUE) dfun else exp(dfun)
        dfun
       }
     else if (type=="right")
      function(x, log = FALSE, ...)
       {
        if (!is.Surv(x)) stop(paste("the x variable is not a Surv object"))
        dfun <- ifelse(x[,"status"]==1, pdf(x[,1],log = TRUE,...),log(1-cdf(x[,1],...)))
        dfun <- if (log == TRUE) dfun else exp(dfun)
        dfun
       } 
     else if (type=="interval")    
      function(x, log = FALSE, ...)
       {
        if (!is.Surv(x)) stop(paste("the x variable is not a Surv object"))
        dfun0 <-ifelse(x[,"status"]==0, cdf(x[,1], lower.tail=F, log.p=T, ...),0) # right  equivalent: log(1-cdf(y[,1],...)) cdf(y[,1], lower.tail=F, log.p=T, ...) 
        dfun1 <-ifelse(x[,"status"]==1, pdf(x[,1],log = TRUE,...),0)# death
        dfun2 <-ifelse(x[,"status"]==2, cdf(x[,1],  log.p=T, ...),0)# left  
        suppressWarnings(dfun3 <-ifelse(x[,"status"]==3, log(cdf(x[,2],...)-cdf(x[,1],...)),0))# interval
        dfun <- dfun0+dfun1+dfun2+dfun3
        dfun <- if (log == TRUE) dfun else exp(dfun)
        dfun
       }
  fun
  }
