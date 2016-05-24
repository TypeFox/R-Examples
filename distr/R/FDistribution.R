
################################
##
## Class: FParameter
##
################################


## Access Methods
setMethod("df1", "FParameter", function(object) object@df1)
setMethod("df2", "FParameter", function(object) object@df2)
setMethod("ncp", "FParameter", function(object) object@ncp)
## Replace Methods
setReplaceMethod("df1", "FParameter", 
                  function(object, value){ object@df1 <- value; object})
setReplaceMethod("df2", "FParameter", 
                  function(object, value){ object@df2 <- value; object})
setReplaceMethod("ncp", "FParameter", 
                  function(object, value){ object@ncp <- value; object})

setValidity("FParameter", function(object){
  if(length(df1(object)) != 1)
    stop("df1 has to be a numeric of length 1")    
  if(df1(object) <= 0)
    stop("df1 has to be positive")
  if(length(df2(object)) != 1)
    stop("df2 has to be a numeric of length 1")    
  if(df2(object) <= 0)
    stop("df2 has to be positive")
  if(length(ncp(object)) != 1)
    stop("ncp has to be a numeric of length 1")      
  else return(TRUE)
})

################################
##
## Class: F distribution
##
################################

Fd <- function(df1 = 1, df2 = 1, ncp = 0) 
               new("Fd", df1 = df1, df2 = df2, ncp = ncp)

## wrapped access methods
setMethod("df1", "Fd", function(object) df1(param(object)))
setMethod("df2", "Fd", function(object) df2(param(object)))
setMethod("ncp", "Fd", function(object) ncp(param(object)))

## wrapped replace methods
setMethod("df1<-", "Fd", function(object, value) 
           new("Fd", df1 = value, df2 = df2(object), ncp=ncp(object)))
setMethod("df2<-", "Fd", function(object, value) 
           new("Fd", df1 = df1(object), df2 = value, ncp=ncp(object)))
setMethod("ncp<-", "Fd", function(object, value) 
           new("Fd", df1 = df1(object), df2 = df2(object), ncp=value))
