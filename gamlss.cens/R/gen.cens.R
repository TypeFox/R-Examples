gen.cens <- function(family = "NO", 
                      name = "cens", 
                      type = c( "right", "left", "interval"), 
                      ...)
 {
  type <- match.arg(type)
   fam  <- as.gamlss.family(family) # family 
    fname <- fam$family[[1]] # family name  
 # fname <- family
 # if (mode(family) != "character" && mode(family) != "name")
 # fname <- as.character(substitute(family))
   dfun <- paste(paste("d",fname,sep=""), substr(type,start=1,stop=1),substr(name,start=1,stop=1), sep="") # say dNOrc
   pfun <- paste(paste("p",fname,sep=""), substr(type,start=1,stop=1),substr(name,start=1,stop=1), sep="") # say pNOrc
   qfun <- paste(paste("q",fname,sep=""), substr(type,start=1,stop=1),substr(name,start=1,stop=1), sep="") # say qNOrc
#   dfun <- paste(paste("d",fname,sep=""), name, sep="")
#   pfun <- paste(paste("p",fname,sep=""), name, sep="")
#   qfun <- paste(paste("q",fname,sep=""), name, sep="")
#   rfun <- paste(paste("r",fname,sep=""), name, sep="")
    fun <- paste(fname, substr(type,start=1,stop=1),substr(name,start=1,stop=1), sep="")
   alldislist <-c(dfun,pfun,qfun,fun)
   # generate d 
   eval(dummy <- cens.d(family = fname, type = type, ...))
   eval(call("<-",as.name(dfun),dummy), envir=parent.frame(n = 1))
   # generate p
   eval(dummy <- cens.p(par, family = fname, type = type, ...))
   eval(call("<-",as.name(pfun),dummy), envir=parent.frame(n = 1))
   # generate q
   eval(dummy <- cens.q(par, family = fname, type = type, ...))
   eval(call("<-",as.name(qfun),dummy), envir=parent.frame(n = 1))   
   # generate the fitting distribution
   eval(dummy <- cens(family = fam, type = type, name=name, local=FALSE, ...))
   eval(call("<-",as.name(fun),dummy), envir=parent.frame(n = 1))
  cat("A censored family of distributions from",  fname, "has been generated \n", 
  "and saved under the names: ", "\n",paste(alldislist, sep=","),"\n")#
  cat("The type of censoring is", type,  " \n") 
 }
#----------------------------------------------------------------------------------------
# this is a dummy function to be able to generate a q function
# needed if centiles are to be form from a given fitted model 
cens.q <-function(family = "NO", ...)
  {
    fname <- family
   if (mode(family) != "character" && mode(family) != "name")
    fname <- as.character(substitute(family))
     qfun <- paste("q",fname,sep="")
   invcdf <- eval(parse(text=qfun))
   invcdf
  }
#----------------------------------------------------------------------------------------
