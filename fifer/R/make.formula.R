##' \code{make.formula} is a function that easily converts a set of strings into a formula. It requires two arguments: a single response variable,
##' and a vector of strings. See examples. 
##'
##' @title Convert strings to a formula
##' @param response a single string used on the left side of a formula
##' @param predictors a string (or a vector of strings) representing the predictors. Each will be separated by a plus sign.
##' @return a formula object
##' @export
##' @author Dustin Fife
##' @examples
##' k = data.frame(matrix(rnorm(100), ncol=5))
##' names(k) = LETTERS[1:5]
##' formula = make.formula("A", LETTERS[2:5])
##' formula
##' lm(formula, data=k)
make.formula = function(response, predictors){
	formula(paste(response, "~", 
				paste(predictors, collapse="+"),
			sep=""))
}

