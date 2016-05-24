## S3 generic for re-fitting a model object with new weights
reweight <- function(object, weights, ...) UseMethod("reweight")

reweight.linearModel <- function(object, weights, ...) {
    fit <- linearModel@fit    
    do.call("fit", c(list(object = object$ModelEnv, weights = weights), object$addargs))
}

reweight.glinearModel <- function(object, weights, ...) {
    fit <- glinearModel@fit    
    do.call("fit", c(list(object = object$ModelEnv, weights = weights), object$addargs))
}

reweight.survReg <- function(object, weights, ...) {
     fit <- survReg@fit
     do.call("fit", c(list(object = object$ModelEnv, weights = weights), object$addargs))
}
