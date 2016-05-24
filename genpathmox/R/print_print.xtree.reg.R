#' @title Print function for the Pathmox Segmentation Trees: linear regression and LAD
#' 
#' @description
#' The function \code{print.xtree.reg} print the \code{reg.pathmox} tree
#' 
#' @param x An object of class \code{"xtree.reg"}.
#' @param \dots Further arguments are ignored.
#'
#' @author Giuseppe Lamberti
#'  
#' @references Lamberti, G. (2014) \emph{Modeling with Heterogeneity.} PhD Dissertation. 
#'
#'
#' \code{\link{summary.xtree.pls}}.
#' @method print xtree.reg
#' @S3method print xtree.reg
#' @examples
#'
#'  \dontrun{
#' #example of LM in alumni satisfaction
#'  
#' data(fibtelereg)
#'
#  #Identify the segmentation variables  
#' segvar= fibtelereg[,2:11]
#'
#  #Select the variables
#' data.fib=fibtelereg[,12:18]          
#'
#  #re-ordering those segmentation variables with ordinal scale
#' segvar$Age 		= factor(segvar$Age, ordered=T)
#' segvar$Salary 	= factor(segvar$Salary, 
#'		levels=c("<18k","25k","35k","45k",">45k"), ordered=T)
#' segvar$Accgrade = factor(segvar$Accgrade, 
#'		levels=c("accnote<7","7-8accnote","accnote>8"), ordered=T)
#' segvar$Grade 	= factor(segvar$Grade, 
#'		levels=c("<6.5note","6.5-7note","7-7.5note",">7.5note"), ordered=T)
#'
#  #Regression PATHMOX
#' fib.reg.pathmox=reg.pathmox(Satisfact~.,data=data.fib,segvar,
#'		signif=0.05,deep=2,method="lm",size=0.15)
#'
#'  print(fib.reg.pathmox)
#'
#'}

print.xtree.reg <- function(x, ...)
{
	cat("\n") 
	print(x$MOX)
	 	
}
