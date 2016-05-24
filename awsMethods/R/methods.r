risk <- function(y, ...) cat("Data extraction not defined for this class:",class(y),"\n")

setGeneric("risk", function(y, ...) standardGeneric("risk"))

extract <- function(x, ...) cat("Data extraction not defined for this class:",class(x),"\n")

setGeneric("extract", function(x, ...) standardGeneric("extract"))

