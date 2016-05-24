##
## Methods for updating an existing object
## =======================================
##
## Method definition for objects of class "NeosOff"
##
setMethod("update", signature(object = "NeosOff"), function(object, formula.,..., evaluate = TRUE){
  call <- object@CALL
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate)
    eval(call, parent.frame())
  else call
})
