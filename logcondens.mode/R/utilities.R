####Charles Doss
##library(logcondens)

## ## Below we use the $CODEHOME bash environment variable.  I set it in bashrc
## ## to differing values on different systems.  NOTE: depending how you start
## ## emacs, the "ESS process" may start a shell which doesn't know about
## ## the bash environment.  If emacs is started from a bash shell (as opposed to
## ## a stand-alone icon) it will apparently know though.  If you start R directly
## ## from a bash shell, R will know the variables.  system("printenv") to check.
## codehome <- system2(command="echo", args="$CODEHOME", stdout=T,stderr=T); ## no "/" on end
## execPath <- paste(codehome,"/Documents/school/research/projects/NPE/logconcave_modeLR/code/", sep="")

## source(paste(execPath,"logcondens_new.R",sep=""))
## #####SHOULD UES dir.exists but because of cluster hvaing trouble
## ##### when i use pvm, i'm not right now ...
## ## {if (dir.exists(execPath)){  
## ##   ##setwd("/home/charles/Documents/school/research/NPE/logconcave_modeLR/")
## ##   source(paste(execPath,"logcondens_new.R",sep=""))
## ## }
## ## else 
## ##   stop("MLEcharacterization.R: Path Error. Perhaps home directory misbehaving or directory structure has changed.")
## ## }



## a should now be a numeric length 1 
x2z <- function(x,w,phi,aval){
  a <- geta(x=x,aval=aval)
  n1 <- length(x)
  if (!a$isx){
    if (a$val > x[n1])     m <- n1
    else if (a$val < x[1]) m <- 1
    else {
      modeNbhrs <- c(phi[a$idx-1],phi[a$idx])
      m <- a$idx-2 + which(modeNbhrs==max(modeNbhrs))
    }
    phi.a <- insert(phi[m],phi,a$idx)
    z=insert(a$val,x, a$idx)
    w.a=insert(0,w,a$idx)
    phi.a <- LocalNormalize(x=z,phi=phi.a)
    return(list(z=z,w.a=w.a,phi.a=phi.a,a=a));    
  }
  else return(list(z=x,w.a=w,phi.a=phi,a=a));
}


##get "a" for unconstrained data. useful for making some functions
##work for both constrained and unconstrained
## so this returns an 'a' object which is the correct mode for the given x and phi.
geta.UC <- function(x,phi){
  k <- which.max(phi)
  a <- list(idx=k, val=x[k], isx=TRUE)
  return(a)
}

##Note that k is 1 or length(x) if aval is outside of range.
## When converted to z, this will remain the same
getk <- function(x,aval){
  n1 <- length(x)
  if (aval<x[1])       k <- 1
  else if (aval>x[n1]) k <- n1+1
  else                 k <- min((1:n1)[x>=aval])
  return(k)
}  

##convert fixed mode location and data to an 'a'/'mode' object.
## 'a' should be numeric, length1.
geta <- function(aval, x){
  k <- getk(x=x,aval=aval)
  if (k<=length(x)) isx <- any((isx <- x[k]==aval) == TRUE)
  else isx <- FALSE
  return(list(val=aval, idx=k, isx=isx))
}


##distinct from rufibach's "intECDF"; returns a function
intECDFfn <- function(x){
  Y <- function(upper,lower=rep(x[1],length(upper))){
    return(intECDF(upper,x) - intECDF(lower,x))
  }
  return(Y)
}

## side is "left" or "right" note that the return value (function) of intFfn
##used to take arguments of the same name but in different order. R issued a
##warning.  Now have changed to comply ...

## But it's better if the "left" function has default of (upper,
## lower=default-arg) so that H(value) gives expected behavior. similarly for
## "right" function. So, changed it back.
intFfn <- function(x,phi,Fhat,prec=1e-10, side="left"){
  if (side=="left"){
    H <- function(upper, lower=rep(x[1],length(upper))){
      return(intF(upper,x,phi,Fhat,prec) - intF(lower,x,phi,Fhat,prec))
    }
  }
  #### This works fine
  else if (side=="right"){
    H <- function(lower, upper=rep(x[length(x)],length(lower))){
      return(intF(-lower,-rev(x),rev(phi),1-rev(Fhat),prec)
             - intF(-upper,-rev(x),rev(phi),1-rev(Fhat),prec)); ##switch, so -lower>-upper
    }
  }
  else
    print("Error: no good value for side passed in. Doing nothing.")
  return(H);
}


