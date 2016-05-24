# ---------------------------------------------------------------------------------------
rqres <- function (pfun = "pNO", 
                   type = c("Continuous", "Discrete", "Mixed"),
               censored = NULL, 
                   ymin = NULL, 
                 mass.p = NULL, 
                prob.mp = NULL,
                      y = y,
                         ... )
{ 
 type <- match.arg(type)
  cdf <- eval(parse(text=pfun))
switch(type, 
  "Continuous"=                   # continous distributions
     {
    rqres <- qnorm(cdf(q=y,...))  # if censored cdf do the right thing
     },  
   "Discrete"=                    # discrete distributions
       {
       # randomize 
        if (is.null(censored)) # uncensored discrete
         {
       aval <- cdf(ifelse(y==ymin,y,y-1), ...) # lower quantile
       aval <- ifelse(y==ymin,0,aval)          # set to zero if y=0
       bval <- cdf(q=y,...)                      # upper quantile
       uval <- runif(length(y),aval,bval) #    gen rand. value
       uval <- ifelse(uval>0.99999,uval-0.1e-10,uval)# 
      rqres <- qnorm(uval)
         }
        else # censured discrete 
         {
        qq <- ifelse(y[,1]==ymin,y[,1], y[,1]-1)    
      aval <- cdf(Surv(qq,y[,2]), ...) # lower quantile
      aval <- ifelse(y[,1]==ymin,0,aval)          # set to zero if y=0
      bval <- cdf(q=y,...)                      # upper quantile
      uval <- runif(length(y[,1]),min=aval,max=bval) #    gen rand. value
      uval <- ifelse(uval>0.99999,uval-0.1e-10,uval)# 
     rqres <- qnorm(ifelse(y[,"status"]==1, uval, bval))
         }
       }, 
   "Mixed"=       # mixed distributions only up to two mass points are allowed
       {
      if (is.null(mass.p) && is.null(prob.mp)) 
              stop("For mixed distributions mass.p and prob.mp arguments have to be specified")
      length.mass.p <- length(mass.p)
      # At the moment we have only the following case
      # for y from 0 to Inf mass point at 0      case 1
      # for y from 0 to 1  mass point at  0      case 1
      # for y from 0 to 1 mass point at 1        case 2
      # for y from 0 to 1 mass point at 0 and 1  case 3
      switch(length.mass.p, 
             {  # if length  1 
                #  and the mass point is at 0    
                # case 1
                if(mass.p==0)
                { 
                 # browser()
                  uval <- ifelse(y==mass.p, runif(length(y),0,prob.mp),cdf(q=y,...))  
                } # if length  1 
                  #  and the mass point is at 2    
                  # case 2
                else if(mass.p==1)
                { 
                  uval <- ifelse(y==mass.p, runif(length(y),1-prob.mp,1),cdf(q=y,...))  
                }
                else 
                {stop("mass point is not at zero or one")}
             },
             {
              #if (mass.p) 
             uval <- ifelse(y==mass.p[1],runif(length(y),0,prob.mp[,1]),0)
             uval <- ifelse(y>mass.p[1] & y<mass.p[2], cdf(q=y,...),uval)
             uval <- ifelse(y==1,runif(length(y),1-prob.mp[,2],1),uval)
             }
             ) 
       rqres <- qnorm(uval)   
       }
       ) 
rqres
}
#----------------------------------------------------------------------------------------
# last change Tuesday, May 22, 2015 MS
rqres1 <- function (obj = NULL, setseed=NULL, save.resid=FALSE, ...)
{
  rqres <- function (pfun = "pNO", 
                     type = c("Continuous", "Discrete", "Mixed"),
                     censored = NULL,  
                     ymin = NULL, 
                     mass.p = NULL, 
                     prob.mp = NULL,
                     y = y,
                     ... )
  { }
  body(rqres) <-  eval(quote(body(rqres)), envir = getNamespace("gamlss"))
  #-------------------------------
  # local function
  #-------------------------------  
  
  newres <- function(obj = NULL, setseed=NULL)
  {  if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
     # if (obj$type!="Discrete" ) stop(paste("This is not discrete distribution ", "\n", ""))
     # if (howmany>10&plot=="all")  stop(paste("You can only have 10 or less plots" , "\n", ""))
     w <- obj$weights
     if (all(trunc(w)==w)) # if frequenies as weights
     { 
       y  <- rep(obj$y, w)
       mu <- rep(fitted(obj, "mu"),w)
       if(any(obj$family%in%.gamlss.bi.list)){ bd <- rep(obj$bd,w)} # MS Wednesday, July 23, 2003 at 12:03   
       if ("sigma"%in%obj$parameters)  sigma <- rep(fitted(obj,"sigma"),w)
       if ("nu"%in%obj$parameters)        nu <- rep(fitted(obj,"nu"),w)
       if ("tau"%in%obj$parameters)      tau <- rep(fitted(obj,"tau"),w)  
     }
     else   # note that weights=1 and weights not frequencies are treated equal here and this could create problems in the future
     {
       y  <- obj$y
       mu <- fitted(obj)
       if(any(obj$family%in%.gamlss.bi.list)){ bd <- obj$bd} # MS Wednesday, July 23, 2003 at 12:03   
       if ("sigma"%in%obj$parameters)  sigma <- fitted(obj,"sigma")
       if ("nu"%in%obj$parameters)        nu <- fitted(obj,"nu")
       if ("tau"%in%obj$parameters)      tau <- fitted(obj,"tau")
     }
     if (!is.null(setseed)) set.seed(setseed)
     res <- eval(obj$rqres)
  }
  
  lobj <- obj$noObs
  res <-  newres(obj, setseed)
  if (save.resid) 
  {
    obj$residuals <- res
    obj 
  } else
    res 
}
