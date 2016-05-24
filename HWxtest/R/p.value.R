# p.value   Generic function to return just the P value from a Hardy-Weinberg test
# (c) William R. Engels, 2014
# Accepts hwtest, list, matrix, vector, table, genotype, etc. 


#' Extract just the P value(s) from a Hardy-Weinberg test.
#' 
#' Use the \code{p.value} function to return just the P value(s) from the results of a call to \code{\link{hwx.test}}. If applied to a list of results, it will return a vector or matrix of P values. You can specify the \code{statName} as \dQuote{LLR}, \dQuote{Prob}, \dQuote{U} or \dQuote{Chisq}. You can also apply \code{p.value} to a matrix or vector and it will attempt to use \code{\link{hwx.test}} to return a P value. However, it's usually preferable to use \code{\link{hwx.test}} directly.
#' 
#' 
#' 
#' @references The methods are described by \href{http://dx.doi.org/10.1534/genetics.109.108977}{Engels, 2009. \bold{Genetics} 183:1431}.
#' 
#' @param x The result of a call to \code{\link{hwx.test}} or a list of such results. It can also be the genotype counts in any of the same formats as accepted by \code{\link{hwx.test}}

#' @param statName can be \dQuote{LLR}, \dQuote{Prob}, \dQuote{U}, or \dQuote{Chisq}
#' 
#' @return The P value
#' 
#' @examples
#' data(HWcases)
#' testResults <- hwx.test(HWcases)
#' p.value(testResults)
#' p.value(testResults, statName="U")

#' 
#' @export
p.value <- function(x, statName=NA){
	 UseMethod("p.value")}

#' @export
p.value.hwtest <- function(x, statName=NA){
	if(is.na(x$method))return(NA)
	if(is.na(statName)) statName <- x$statName;
	x$Pvalues[statName]
}
	
	
#' @export
p.value.list <- 
function(x,statName= NA) {sapply(x, p.value, statName=statName)}

#' @export
p.value.matrix <- 
function(x, statName="LLR") {hwx.test(x, statName=statName)$Pvalues[statName]} 

#' @export
p.value.integer <- 
function(x, statName="LLR"){ hwx.test(x, statName=statName)$Pvalues[statName]}

#' @export
p.value.table <- 
function(x, statName="LLR") {hwx.test(x, statName=statName)$Pvalues[statName]}

#' @export
p.value.numeric <- 
function(x, statName="LLR") {hwx.test(x, statName=statName)$Pvalues[statName]}

#' @export
p.value.genotype <- 
function(x, statName="LLR") {hwx.test(x, statName=statName)$Pvalues[statName]}

#' @export
p.value.logical <- 
function(x, statName="LLR") {hwx.test(x, statName=statName)$Pvalues[statName]}