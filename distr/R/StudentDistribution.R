################################
##
## Class: TParameter
##
################################

## Access Methods
setMethod("df", "TParameter", function(x, ...) x@df)
setMethod("ncp", "TParameter", function(object) object@ncp)
## Replace Methods
setReplaceMethod("df", "TParameter", 
                  function(object, value){ object@df <- value; object})
setReplaceMethod("ncp", "TParameter", 
                  function(object, value){ object@ncp <- value; object})


setValidity("TParameter", function(object){
  if(length(df(object)) != 1)
    stop("df has to be a numeric of length 1")    
  if(df(object) <= 0)
    stop("df has to be positive")
  if(length(ncp(object)) != 1)
    stop("ncp has to be a numeric of length 1")      
  else return(TRUE)
})

################################
##
## Class: Student distribution
##
################################

Td <- function(df = 1, ncp = 0) new("Td", df = df, ncp = ncp)

## wrapped access methods
setMethod("df", "Td", function(x, ...) df(param(x)))
setMethod("ncp", "Td", function(object) ncp(param(object)))

## wrapped replace methods
setMethod("df<-", "Td", 
           function(object, value) new("Td", df = value, ncp = ncp(object)))
setMethod("ncp<-", "Td", 
           function(object, value) new("Td", df = df(object), ncp = value))
