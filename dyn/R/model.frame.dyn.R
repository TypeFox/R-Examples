model.frame.dyn <- function (formula, data = NULL, ...)
{
	if (!inherits(formula, "formula")) {
		if (nargs() == 1) return(NextMethod("model.frame", formula))
		return(NextMethod("model.frame", formula, data = data, ...))
	}
	class(formula) <- "formula"
	if (is.null(data)) data <- parent.frame()
	env <- environment(formula)
	.Class <- if (is.environment(data) || is.list(data))
			class(eval(as.list(formula)[[2]], data, env))
		  else class(data)
	mf <- NextMethod("model.frame", formula, data = data, ...)
	attr(mf, ".Class") <- .Class
	predvars <- attr(attr(mf, "terms"), "predvars")
	attr(mf, "series") <- eval(predvars, as.list(data), env)
	for(i in seq(mf)) 
		if (!inherits(mf[[i]], "factor")) 
			# mf[[i]] <- as.vector(unclass(mf[[i]]))
			mf[[i]] <- coredata(mf[[i]])
	mf
}

model.frame.zooreg <- model.frame.ts <- model.frame.irts <- model.frame.its <-
model.frame.zoo <- function (formula, data = NULL, ...)
{
	if (is.environment(data))
		formula <- terms(formula, data = data)
	else {
		data <- as.list(data)
		data1 <- as.data.frame(lapply(data, "[", 1))
		names(data1) <- names(data)
		formula <- terms(formula, data = data1)
	}
	cl <- attr(formula, "variables")
	cl[[1]] <- as.name("merge.zoo")
	cl$retclass <- "data.frame"
	cl$retclass <- "list"
	cl$all <- TRUE
	attr(formula, "predvars") <- as.call(cl)
	NextMethod("model.frame", formula, data = data, ...)
}

model.matrix.dyn <- function(object, ...) {
	model.matrix(terms(object), model.frame(object), ...)
}

