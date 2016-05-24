gen.trun <-function(par = c(0), 
                     family = "NO", 
                     name = "tr", 
                     type = c("left", "right", "both"),
                  varying = FALSE,
                     ...)
 {
  type <- match.arg(type)
   #famAsname <- deparse(substitute(family))
   fam  <- as.gamlss.family(family) #ds Monday, March 10, 2008 at 10:07
  fname <- fam$family[[1]] 
 # fname <- family
 # if (mode(family) != "character" && mode(family) != "name")
 # fname <- as.character(substitute(family))
   dfun <- paste(paste("d",fname,sep=""), name, sep="")
   pfun <- paste(paste("p",fname,sep=""), name, sep="")
   qfun <- paste(paste("q",fname,sep=""), name, sep="")
   rfun <- paste(paste("r",fname,sep=""), name, sep="")
    fun <- paste(fname, name, sep="")
   alldislist <-c(dfun,pfun,qfun,rfun,fun)
   # generate d 
   eval(dummy <- trun.d(par, family = fname, type = type, varying = varying, ...))
   eval(call("<-",as.name(dfun),dummy), envir=parent.frame(n = 1))
   # generate p
   eval(dummy <- trun.p(par, family = fname, type = type, varying = varying, ...))
   eval(call("<-",as.name(pfun),dummy), envir=parent.frame(n = 1))
   # generate q
   eval(dummy <- trun.q(par, family = fname, type = type, varying = varying, ...))
   eval(call("<-",as.name(qfun),dummy), envir=parent.frame(n = 1))
   # generate r
   eval(dummy <- trun.r(par, family = fname, type = type, varying = varying, ...))
   eval(call("<-",as.name(rfun),dummy), envir=parent.frame(n = 1))
   # generate the fitting distribution
   eval(dummy <- trun(par, family = substitute(family), type = type, name=name, local=FALSE, varying = varying, ...))
   eval(call("<-",as.name(fun),dummy), envir=parent.frame(n = 1))
  cat("A truncated family of distributions from",  fname, "has been generated \n", 
  "and saved under the names: ", "\n",paste(alldislist,sep=","),"\n")#
  cat("The type of truncation is", type, "and the truncation parameter is", par, " \n") 
 }
