#  Calculate Effects for term(s) in a Multivariate Linear Model 
#  2014-03-12: Introduced allEffects.mlm(). J. Fox


Effect.mlm <- function(focal.predictors, mod, response, ...) {
	if (missing(response)) {
		mod.frame <- model.frame(mod)
    response <- colnames(model.response(mod.frame))
	}
	else if (is.numeric(response)) {
		mod.frame <- model.frame(mod)
    response.names <- colnames(model.response(mod.frame))
    response <- response.names[response]
	}
	
	if (length(response)==1) {
			mod.1 <- update(mod, as.formula(paste(response, " ~ .")))
			result <- Effect(focal.predictors, mod.1,  ...)
	}
	else {
		result <- as.list(NULL)
		for (resp in response) {
			mod.1 <- update(mod, as.formula(paste(resp, " ~ .")))
			lab <- resp
			result[[lab]] <- Effect(focal.predictors, mod.1,  ...)
		}
		class(result) <- "efflist"
	}
	result
}

allEffects.mlm <- function(mod, ...){
    result <- NextMethod()
    class(result) <- "mlm.efflist"
    result
}

plot.mlm.efflist <- function(x, ...){
    x <- do.call(c, x)
    class(x) <- "efflist"
    plot(x, ...)
}

summary.mlm.efflist <- function(object, ...){
    object <- do.call(c, object)
    for (effect in names(object)){
        cat("\n\nResponse:", object[[effect]]$response, "\n")
        print(summary(object[[effect]], ...))
    }
}

print.mlm.efflist <- function(x, ...){
    x <- do.call(c, x)
    for (effect in names(x)){
        cat("\n\nResponse:", x[[effect]]$response, "\n")
        print(x[[effect]], ...)
    }
    invisible(x) 
}

