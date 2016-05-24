##
## Methods for objects of class PortSol
##
## Weights-method
##
setMethod(f = "Weights", signature(object = "PortSol"), definition = function(object){
  ans <- slot(object, "weights")
  return(ans)
})
##
## show-method
##
setMethod(f = "show", signature(object = "PortSol"), definition = function(object){
  cat("\n")
  cat("\nOptimal weights for porfolio of type:\n")
  cat(paste(object@type))
  cat("\n")
  cat("\n")
  print(round(Weights(object), 4))
  cat("\n")
})
## 
## Solution-method
##
setMethod(f = "Solution", signature(object = "PortSol"), definition = function(object){
  ans <- slot(object, "opt")
  return(ans)
})
##
## update-method
##
setMethod("update", signature(object = "PortSol"), function(object, formula.,..., evaluate = TRUE){
  call <- object@call
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

