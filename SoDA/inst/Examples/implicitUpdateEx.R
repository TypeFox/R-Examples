updFormula <- function(model, var) {
       var <- as.character(substitute(var))
       tnames <- colnames(attr(terms(model), "factors"))
       hasVar <- grep(paste(":",var,":", sep = ""),
                      paste(":", tnames, ":", sep=""), fixed = TRUE)
       fText <- paste("~ .", paste("-", tnames[hasVar], collapse = " "))
       eval(parse(text = fText)[[1]])
   }

setGeneric("updFormula", signature = "model")

setMethod(updFormula, "formula", function(model, var) {
       eval(parse(text = paste("~ . -", as.character(var))))
   })

setGenericImplicit("updFormula")
