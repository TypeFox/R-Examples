###
# Abstract Model
# 
# used to represent underlying model/predict calls in a generic way
###
setClass("AbstractModel",
	representation(
		model="character",
		model_args="list",
		predict="function",
		predict_args="list",
		xtrans="function",
		ytrans="function",
		formula="logical"
	),
	prototype = prototype(xtrans=function(X){X}, ytrans=function(X){X})
)

ab.create <- function(model.call, model.args=list(), predict.args=list(), formula=NA,...)
{
	#see if the first arg of the model function is formula, and if so assume we use formula
	if (is.na(formula)) {
		formula <- names(formals(model.call))[1] == 'formula'
	}
	
	b <- new("AbstractModel", model=model.call, model_args=model.args, predict_args=predict.args, formula=formula, ...) 
	return(b)
}

ab.arglist <- function(ab, X, Y)
{
	x <- ab@xtrans(X)
	y <- ab@ytrans(Y)
	if (ab@formula) {
		arglist <- list(Y~., data=cbind(x, Y=y))
	} else {
		arglist <- list(x=x, y=y)
	}
	arglist <- c(arglist, ab@model_args)
}

ab.predict <- function(ab, mod, X)
{
	x = ab@xtrans(X)
	alist <- list(mod, X)
	alist <- c(alist, ab@predict_args)
	r <- do.call(predict, alist)
}

ab.model <- function(ab, X, Y)
{
	k <- do.call(ab@model, ab.arglist(ab, X, Y))	
}

ab.prargs <- function(x, ...) {
	if (!is.list(x)) {
		print(paste("x =", x))
	} else {
		print(unlist(x))
	}
}

###
# print for the abstract model class
# move this to AbstractModel.R and also fill out details
###
print.AbstractModel <- function(x, ...) { cat('ab model print')}
