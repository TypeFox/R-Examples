# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March 2016
# Version 1.1
# Licence GPL v3

# 
# if (!isGeneric("names")) {
#   setGeneric("names", function(x)
#     standardGeneric("names"))
# }

setMethod('names', signature(x='sdmdata'), 
          function(x) {
            x@species.names
          }
)
# 
# if (!isGeneric("names<-")) {
#   setGeneric("names<-", function(x,value)
#     standardGeneric("names<-"))
# }


setReplaceMethod('names', signature(x='sdmdata'), 
                 function(x,value) {
                   if (!is.character(value)) stop('value should be a character vector')
                   if (length(x@species.names) != length(value)) stop('the length of names is not equal to the number of species...!')
                   x@species.names <- value
                   names(x@species) <- value
                   if (length(x@sdmFormula@formula) == 3) {
                     x@sdmFormula@formula <- as.formula(paste(paste(value,collapse='+'),'~',deparse(x@sdmFormula@formula[[3]]),sep=''),env = parent.frame())
                   } else if (length(x@sdmFormula@formula) == 2 && deparse(x@sdmFormula@formula[[1]]) == '~') {
                     x@sdmFormula@formula <- as.formula(paste(paste(value,collapse='+'),'~',deparse(x@sdmFormula@formula[[2]]),sep=''),env = parent.frame())
                   } else {
                     x@sdmFormula@formula <- as.formula(paste(paste(value,collapse='+'),'~',paste(x@sdmFormula@vars,collapse='+'),sep=''),env = parent.frame())
                     warning('formula in the @sdmFormula may have problem!')
                   }
                   x@sdmFormula@species <- value
                   x
                 }
)