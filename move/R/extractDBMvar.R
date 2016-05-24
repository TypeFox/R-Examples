setGeneric(".extractDBMvar", function(object) {standardGeneric(".extractDBMvar")})
setMethod(f = ".extractDBMvar", 
          signature = "list", 
          definition = function(object) {
            new()
            means <- unlist(lapply(DBMvarLST, slot, "means"))
            in.windows <- c(in.windows, object[[i]]@DBMvar@in.windows)
            interest <- c(interest, object[[i]]@DBMvar@interest)
            
            return(dBMvarianceStack)
          })


setGeneric(".extractDBMvar", function(object) {
  standardGeneric(".extractDBMvar")
})
setMethod(f = ".extractDBMvar", signature = "DBBMM", definition = function(object) {
  return(object@DBMvar)
})
