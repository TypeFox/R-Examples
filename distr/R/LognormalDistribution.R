
################################
##
## Class: LnormParameter
##
################################

## Access Methods
setMethod("meanlog", "LnormParameter", function(object) object@meanlog)
setMethod("sdlog", "LnormParameter", function(object) object@sdlog)
## Replace Methods
setReplaceMethod("meanlog", "LnormParameter", 
                  function(object, value){ object@meanlog <- value; object})
setReplaceMethod("sdlog", "LnormParameter", 
                  function(object, value){ object@sdlog <- value; object})

setValidity("LnormParameter", function(object){
  if(length(sdlog(object)) != 1)
    stop("sdlog has to be a numeric of length 1")    
  if(sdlog(object) <= 0)
    stop("sdlog has to be positive")
  else return(TRUE)
}           
)

################################
##
## Class: lognormal distribution
##
################################

Lnorm <- function(meanlog = 0, sdlog = 1) 
             new("Lnorm", meanlog = meanlog, sdlog = sdlog)

## wrapped access methods
setMethod("meanlog", "Lnorm", function(object) meanlog(param(object)))
setMethod("sdlog", "Lnorm", function(object) sdlog(param(object)))

## wrapped replace methods
setMethod("meanlog<-", "Lnorm", 
           function(object, value) 
               new("Lnorm", meanlog = value, sdlog = sdlog(object)))
setMethod("sdlog<-", "Lnorm", 
           function(object, value) 
               new("Lnorm", meanlog = meanlog(object), sdlog = value))

setMethod("*", c("Lnorm","numeric"),
       function(e1, e2){
         if (length(e2)>1) stop("length of operator must be 1")
         if(isTRUE(all.equal(e2,0)))  
            return(new("Dirac", location = 0, .withArith = TRUE))
         nL <- new("Lnorm", meanlog = meanlog(e1) + log(abs(e2)), 
                        sdlog = sdlog(e1), .withArith = TRUE)
         if(e2 > 0) 
            return(nL)
         else 
            return(getMethod("*", c("AbscontDistribution","numeric"))(e1,e2))
            #return(-1 * as(nL, "AbscontDistribution"))
       })
