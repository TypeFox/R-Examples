##
## Methods for updating an existing object
## =======================================
##
## Method definition for objects of class "GoGARCH"
##
setMethod("update", signature(object = "GoGARCH"), function(object, formula.,..., evaluate = TRUE){
  call <- object@CALL
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.))
    call$formula <- update.formula(formula(object), formula.)
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
##
## Method definition for objects of class "Goestica"
## "Goestica" extends directly "GoGARCH"
##
setMethod("update", signature(object = "Goestica"), function(object, formula., ..., evaluate = TRUE){
  call <- object@CALL
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.))
    call$formula <- update.formula(formula(object), formula.)
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
##
## Method definition for objects of class "Goestmm"
## "Goestmm" extends directly "GoGARCH"
##
setMethod("update", signature(object = "Goestmm"), function(object, formula., ..., evaluate = TRUE){
  call <- object@CALL
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.))
    call$formula <- update.formula(formula(object), formula.)
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
##
## Method definition for objects of class "Goestnls"
## "Goestnls" extends directly "GoGARCH"
##
setMethod("update", signature(object = "Goestnls"), function(object, formula., ..., evaluate = TRUE){
  call <- object@CALL
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.))
    call$formula <- update.formula(formula(object), formula.)
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
##
## Method definition for objects of class "Goestml"
## "Goestml" extends directly "GoGARCH"
##
setMethod("update", signature(object = "Goestml"), function(object, formula., ..., evaluate = TRUE){
  call <- object@CALL
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.))
    call$formula <- update.formula(formula(object), formula.)
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
