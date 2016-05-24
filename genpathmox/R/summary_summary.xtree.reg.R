#' @title Summary function for the Pathmox Segmentation Trees: linaer regression and LAD
#' 
#' @description
#' The function \code{summary.xtree.reg} returns the most important results obtained 
#' by the function \code{reg.pathmox}. In order, it provides the parameters algorithm (
#' threshold significance,node size limit,tree depth level and the method used for the 
#' split partition), the basic characteristics of the tree (deep and number of terminal 
#' nodes), the basic characteristics of the nodes and the F-global and F-coefficients 
#' results. For the test results the significance level is indicated.
#' 
#' @param object An object of class \code{"xtree.reg"}.
#' @param \dots Further arguments are ignored.
#'
#' @author Giuseppe Lamberti
#'  
#' @references Lamberti, G. (2014) \emph{Modeling with Heterogeneity.} PhD Dissertation. 
#'
#'
#' \code{\link{summary.xtree.pls}}, \code{\link{reg.pathmox}}.
#'
#' @method summary xtree.reg
#' @S3method summary xtree.reg
#' 
#' 
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
#' data.fib= fibtelereg[,12:18]          
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
#'  summary(fib.reg.pathmox)
#'
#'}


summary.xtree.reg <- function(object, ...)
{
	
	info= object$model
	
	cat("\n")
	cat("REGRESSION SEGMENTATION TREE","\n")
	cat("\n")
	cat("---------------------------------------------")
	cat("\n")
	cat("Info Parameters Algorithm:","\n")
	info.value=rbind(info[[1]],info[[2]],info[[3]],info[[4]])
	dimnames(info.value)=NULL
	info.name= c("Threshold signif","Node size limit(%)","Tree depth level","Method")
	info.tree=data.frame(info.name,info.value)
	names(info.tree)= c("Parameters Algorithm", "value")
	print(info.tree)
	cat("---------------------------------------------")
	cat("\n")
	cat("Info Tree:","\n")
	tree= object$MOX
	info.value=rbind(max(tree[,3]),sum(length(which(tree[,5]=="yes"))))
	dimnames(info.value)=NULL
	info.name= c("Deep tree","Number terminal nodes")
	info.tree=data.frame(info.name,info.value)
	names(info.tree)= c("Parameters Tree", "value")
	print(info.tree)
	cat("---------------------------------------------")
	cat("\n")
	cat("Info nodes:","\n")
	print(tree)
	cat("---------------------------------------------")
	cat("\n")
	cat("Info Splits:","\n")
	cat("\n")
	cat("Variable:","\n")
	print(object$Fg.r[,c(1,4,5,6)])
	cat("\n")
	cat("F.statistic global:","\n")
	printCoefmat(object$Fg.r[,c(1,2,3)], P.values=TRUE, has.Pvalue=TRUE)
	cat("\n")
	cat("F.statistic coefficient:","\n")
	
	for(i in 1:length(object$Fc.r)) {
		cat("\n")
		cat(paste("Node",i,":"))
		cat("\n")
		printCoefmat(object$Fc.r[[i]], P.values=TRUE, has.Pvalue=TRUE)}
}
