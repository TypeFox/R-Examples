##' Generic function for obtaining scaled coefficients
##'
##' Given an object of class lm, glm, or lda, this function will first standardize the variables, then run the model again. The resulting
##' coefficients will be standardized betas. 
##' @title Standardize coefficients
##' @param object an object resulting from glm, lm, or lda
##' @param scale.response should the response variable be scaled as well? (Usually not for glm or lda).
##' @return an object of the same class as the one outputted
##' @export
##' @author Dustin Fife
##' @examples
##' 		#### create random data with different means and variances
##' d = data.frame(matrix(rnorm(5*50, c(10,5,14,100, 33), c(3,5,4,3,5)), nrow=50, byrow=TRUE))
##' names(d) = LETTERS[1:5]
##' g = lm(C~B + A + D + E, data=d)
##' scaleB(g, TRUE)
##' ##### make a logistic
##' d$A = as.factor(as.numeric(cut(d$A, 2, labels=c(1,0))	))
##' object = glm(A~B + C + D + E, data=d, family=binomial)
##' scaleB(object)
##' ##### LDA
##' object = lda(A~B + C + D + E, data=d, family=binomial)
##' scaleB(object)
scaleB = function(object, scale.response=F){
	vars = row.names(attr(terms(object), "factors"))
	if (scale.response){
		part1 = paste("scale(", vars[1], ")~", sep="")
	} else {
		part1 = paste("(", vars[1], ")~", sep="")
	}
	new.form = formula(paste(part1,
					paste("scale(", vars[-1], ")", collapse="+", sep=""),
					sep=""))
	n.object = update(object, new.form, evaluate=FALSE)	
	n.object = eval.parent(n.object)		
	return(n.object)
}


