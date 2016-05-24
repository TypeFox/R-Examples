cens.p <- function(family = "NO", type = c( "right", "left", "interval"),...)
  {
    #if (!is.Surv(y)) stop(paste("the y variable is not a Surv object"))
     type <- match.arg(type)
    if (type=="counting") stop(paste("gamlss has no not support counting"))
   fname <- family
   if (mode(family) != "character" && mode(family) != "name")
   fname <- as.character(substitute(family))
 distype <- eval(call(family))$type
    dfun <- paste("d",fname,sep="")
    pfun <- paste("p",fname,sep="")
     pdf <- eval(parse(text=dfun))
     cdf <- eval(parse(text=pfun))
fun <- if (type=="left")  
       function(q, log = FALSE, ...)
        {
         if (!is.Surv(q)) stop(paste("the q variable is not a Surv object"))
        pfun1 <- cdf(q[,1],...)
        pfun2 <- runif(length(q[,1]),0,pfun1) 
        pfun <- ifelse(q[,"status"]==1, pfun1,pfun2)
        pfun
       }
     else if (type=="right")
      function(q, log = FALSE, ...)
       {
        if (!is.Surv(q)) stop(paste("the q variable is not a Surv object"))
        pfun1 <- cdf(q[,1],...)
        pfun2 <- runif(length(q[,1]),pfun1,1) 
        pfun <- ifelse(q[,"status"]==1, pfun1,pfun2)
        pfun
       } 
     else if (type=="interval")    
      function(q, log = FALSE, ...)
       {# needs checking
        if (!is.Surv(q)) stop(paste("the q variable is not a Surv object"))
        pfun1 <- cdf(q[,1],...)
        pfun2 <- runif(length(q[,1]),0,pfun1) 
        pfun0 <- runif(length(q[,1]),pfun1,1) 
        suppressWarnings(pfun3<-runif( length(q[,1]), cdf(q[,1],...), cdf(q[,2],...) ) )# interval   )
        pfun0 <-ifelse(q[,"status"]==0, pfun0,0)  # right 
        pfun1 <-ifelse(q[,"status"]==1, pfun1,0)  # death
        pfun2 <-ifelse(q[,"status"]==2, pfun2,0)  # left
        pfun3 <-ifelse(q[,"status"]==3, pfun3,0)  # interval         
        dfun <- pfun0+pfun1+pfun2+pfun3
        dfun
       }
  fun
  }
