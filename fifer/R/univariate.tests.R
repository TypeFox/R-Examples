comp.p = function(variable, group){
	#### if it's a factor, compute chi square test
	if (is.factor(variable) | is.character(variable) | length(unique(variable))<4){
		variable = as.factor(variable)
		test = chisq.test(variable,group)
		return(test$p.value)
	} else if (length(unique(group))<3){
	##### do a t test
		test = t.test(variable~group)
		return(test$p.value)
	} else {
		k = lm(variable~group)
		return(as.numeric(unlist(anova(k)["Pr(>F)"]))[1])
	}
}

##' Extract the p value from a univariate significance test
##'
##' @description Given a dataframe, this function predicts the specified categorical variable using each column in the dataset, one at a time.
##' The function will automatically select whether to do a chi-square test, a t-test, or an ANOVA. See details. 
##' @details \code{univariate.tests} will look at each column in the dataframe, then perform a t-test, ANOVA, or chi-square test where
##' the grouping variable serves as the independent variable. The computer will chose a chi-square test of one of the following three conditions is
##' met: (1) the variable is a factor, (2) the variable is a character variable, or (3) the variable has less than four unique values. An ANOVA will be
##' used if the number of levels of the grouping variable is greater than two. In all other cases, a t-test will be used. 
##' @title Extract p values for a data frame
##' @param dataframe a data frame containing both the variables and the grouping variable
##' @param exclude.cols a vector indicating which columns should not have a significance test
##' @param group a string with the name of the grouping variable
##' @return a vector of p values
##' @author Dustin Fife
##' @export
##' @examples
##' k = data.frame(cbind(ID=1:100,
##'				A = rnorm(100),
##'				B = rnorm(100),
##'				C = rnorm(100),
##'				Group = rep(1:2, times=50)))
##' univariate.tests(dataframe = k, exclude.cols=1, group="Group")
univariate.tests = function(dataframe, exclude.cols=NULL, group){
	cl = which(names(dataframe)==group)
	if (is.null(exclude.cols)){
		exclude.cols = cl
	} else {
		exclude.cols = c(exclude.cols, cl)
	}
	p.vals = unlist(lapply(dataframe[,-c(exclude.cols)], comp.p, group=dataframe[,group]))
	return(p.vals)
}


#' Fictitious Medical Dataset
#'
#' @name fakeMedicalData
#' @docType data
#' @author Dustin Fife \email{blahblah@@roxygen.org}
#' @keywords data
NULL