## S4 StatModel object
btReg <- function(type = "loglin", ref = NULL, undecided = NULL, position = NULL) {
  new("StatModel",
    capabilities = new("StatModelCapabilities"),
    name = "Bradley-Terry regression model",
    dpp = ModelEnvFormula,
    fit = function(object, weights = NULL, ...){

        ## extract response (there are no regressors)
        y <- object@get("response")[[1]]
	
        ## call btReg.fit()
        z <- btReg.fit(y = y, weights = weights, type = type, ref = ref,
	  undecided = undecided, position = position)
	z$ModelEnv <- object
	z$addargs <- list(...)	
        z
    }
  )
}

reweight.btReg <- function(object, weights, ...) {
     fit <- btReg(type = object$type, ref = object$ref,
       undecided = object$undecided, position = object$position)@fit
     do.call("fit", c(list(object = object$ModelEnv, weights = weights), object$addargs))
}

estfun.btReg <- function(x, ...) x$estfun

bread.btReg <- function(x, ...) x$vcov * x$n

