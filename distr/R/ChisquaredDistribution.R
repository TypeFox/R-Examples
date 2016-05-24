
################################
##
## Class: ChisqParameter
##
################################

## Access Methods
setMethod("df", "ChisqParameter", function(x, df1, df2, log = FALSE) x@df)
setMethod("ncp", "ChisqParameter", function(object) object@ncp)
## Replace Methods
setReplaceMethod("df", "ChisqParameter", function(object, value)
                 { object@df <- value; object})
setReplaceMethod("ncp", "ChisqParameter", function(object, value)
                 { object@ncp <- value; object})


setValidity("ChisqParameter", function(object){
  if(length(df(object)) != 1)
    stop("df has to be a numeric of length 1")    
  if(df(object) <= 0 + .Machine$double.eps^.5 )
    stop("df has to be positive")
  if(length(ncp(object)) != 1)
    stop("ncp has to be a numeric of length 1")    
  if(ncp(object) < 0)
    stop("ncp has to be not negative")
  else return(TRUE)  })


################################
##
## Class: Chi squared distribution
##
################################

Chisq <- function(df = 1, ncp = 0) new("Chisq", df = df, ncp = ncp)

## wrapped access methods
setMethod("df", "Chisq", function(x, df1, df2, log = FALSE) df(param(x)))
setMethod("ncp", "Chisq", function(object) ncp(param(object)))

## wrapped replace methods
setMethod("df<-", "Chisq", 
           function(object, value) 
               new("Chisq", df = value, ncp = ncp(object)))
setMethod("ncp<-", "Chisq", 
           function(object, value) new("Chisq", df = df(object), ncp = value))

setMethod("+", c("Chisq","Chisq"),
          function(e1,e2){
            newdf <- df(e1) + df(e2)
            newncp <- ncp(e1) + ncp(e2)
            return(new("Chisq", df = newdf, ncp = newncp, .withArith = TRUE))
          })

### new from version 1.9 on:

setMethod("shape", "Chisq", 
           function(object){ 
               if (isTRUE(all.equal(ncp(object),0))) 
                   df(object)/2
               else stop(gettextf("%s is not a Gamma distribution"),
                                   deparse(substitute(object))
                         )
                            }
           )

setMethod("scale", "Chisq", 
           function(x, center = TRUE, scale = TRUE){ 
               if (isTRUE(all.equal(ncp(x),0))) 
                                   2
               else stop(gettextf("%s is not a Gamma distribution"),
                                   deparse(substitute(object))
                         )
                            }
           )

###### end new version 1.9

