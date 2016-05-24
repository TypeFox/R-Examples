model.frame.gls <- function(object, ...) model.frame(formula(object), nlme::getData(object),...)
	
model.frame.lme <- function(object, ...) object$data

model.matrix.aovlist <- function(object, ...)
    stop(sQuote("predicted.means"), " does not support objects of class ",
         sQuote("aovlist"))
		 
model.matrix.gls <- function(object, ...)
    model.matrix(terms(object), data = nlme::getData(object), ...)
	
model.matrix.lme <- function(object, ...)
    model.matrix(terms(object), data = model.frame(object), ...)

terms.gls <- function(object, ...) terms(model.frame(object),...) 

terms.lme <- function (object, ...) {
  v <- object$terms
  if (is.null(v)) stop("no terms component")
  return(v)
}

#terms.merMod <- function (object, ...) {
#  v <- terms(object)
#  attr(v, "dataClasses") <- sapply(all.vars(formula(object, fixed.only=TRUE)),function(x) class(model.frame(object)[[x]])[1])
#  return(v)
#}
#	