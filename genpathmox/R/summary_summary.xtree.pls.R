#' @title Summary function for the Pathmox Segmentation Trees: PLS-PM
#' 
#' @description
#' The function \code{summary.xtree.pls} returns the most important results obtained 
#' by the function \code{pls.pathmox}. In order, it provides the parameters algorithm (
#' threshold significance,node size limit",tree depth level and the method used for the 
#' split partition), the basic characteristics of the tree (deep and number of terminal 
#' nodes), the basic characteristics of the nodes and the F-global the F-block and F-coefficients 
#' results. For the test results the significance level is indicated.
#' 
#' @param object An object of class \code{"xtree.pls"}.
#' @param \dots Further arguments are ignored.
#'
#' @author Giuseppe Lamberti
#'  
#' @references Lamberti, G. (2014) \emph{Modeling with Heterogeneity.} PhD Dissertation. 
#'
#'
#' \code{\link{pls.pathmox}}
#' @method summary xtree.pls
#' @S3method summary xtree.pls
#' @examples
#'  
#'  \dontrun{
#'  ## example of PLS-PM in alumni satisfaction
#'  
#'  data(fibtele)
#'  
#'  # select manifest variables
#'  data.fib <-fibtele[,12:35]
#'  
#'  # define inner model matrix
#'  Image 			= rep(0,5)
#'	 Qual.spec	= rep(0,5)
#'	 Qual.gen		= rep(0,5)
#'	 Value			= c(1,1,1,0,0)
#'	 Satis			= c(1,1,1,1,0)
#'  inner.fib <- rbind(Image,Qual.spec, Qual.gen, Value, Satis)
#'  colnames(inner.fib) <- rownames(inner.fib)
#'  
#'  # blocks of indicators (outer model)
#'  outer.fib <- list(1:8,9:11,12:16,17:20,21:24)
#'  modes.fib  = rep("A", 5)
#'  
#'  # apply plspm
#'  pls.fib <- plspm(data.fib, inner.fib, outer.fib, modes.fib)
#'                  
#'
#'  # re-ordering those segmentation variables with ordinal scale 
#'   seg.fib= fibtele[,2:11]
#'  
#'	 seg.fib$Age = factor(seg.fib$Age, ordered=T)
#'	 seg.fib$Salary = factor(seg.fib$Salary, 
#'			levels=c("<18k","25k","35k","45k",">45k"), ordered=T)
#'	 seg.fib$Accgrade = factor(seg.fib$Accgrade, 
#'			levels=c("accnote<7","7-8accnote","accnote>8"), ordered=T)
#'	 seg.fib$Grade = factor(seg.fib$Grade, 
#'	levels=c("<6.5note","6.5-7note","7-7.5note",">7.5note"), ordered=T)
#'
#'  # Pathmox Analysis
#'  fib.pathmox=pls.pathmox(pls.fib,seg.fib,signif=0.05,
#'					deep=2,size=0.2,n.node=20)
#'  
#'
#'  summary(fib.pathmox)
#'  }

summary.xtree.pls <- function(object, ...) 
{
	info=object$model
	
	cat("\n")
	cat("PLS-PM SEGMENTATION TREE","\n")
	cat("\n")
	cat("---------------------------------------------")
	cat("\n")
	cat("Info Parameters Algorithm:","\n")
	info.value = rbind(info[[1]],info[[2]],info[[3]])
	dimnames(info.value) = NULL
	info.name = c("Threshold signif","Node size limit(%)","Tree depth level")
	info.tree = data.frame(info.name,info.value)
	names(info.tree) = c("Parameters Algorithm", "value")
	print(info.tree)
	cat("---------------------------------------------")
	cat("\n")
	cat("Info Tree:","\n")
	tree = object$MOX
	info.value = rbind(max(tree[,3]),sum(length(which(tree[,5]=="yes"))))
	dimnames(info.value) = NULL
	info.name = c("Deep tree","Number terminal nodes")
	info.tree = data.frame(info.name,info.value)
	names(info.tree) = c("Parameters Tree", "value")
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
	cat("F.statistic block:","\n")
	print(object$Fb.r)
	cat("F.statistic coefficient:","\n")
	for (i in 1:length(object$Fc.r))
	{
		cat("\n")
		cat(paste("Node",substr(names(object$Fc.r)[i],5,5),":"))
		cat("\n")
		printCoefmat(object$Fc.r[[i]], P.values=TRUE, has.Pvalue=TRUE)
	}
}
