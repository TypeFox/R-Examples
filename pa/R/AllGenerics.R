################################################################################
##
## Generics for the pa package.
##
################################################################################

## portfolioMatch class methods.

if(!isGeneric("exposure"))
  setGeneric("exposure", function(object,
                                  var, 
                                  ...) standardGeneric("exposure"))

if(!isGeneric("returns"))
  setGeneric("returns", function(object,
                                 ...) standardGeneric("returns"))

