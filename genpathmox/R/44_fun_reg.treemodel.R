	
#' @title Regression results of terminal nodes from the Pathmox Segmentation Trees
#' 
#' @description
#' Calculates basic regression results for the terminal nodes of Pathmox Segmentation Trees: liner regression and LAD
#' trees
#' 
#' @details
#' 
#' The argument \code{xtree.reg} is an object of class \code{"xtree.reg"} returned by 
#'\code{\link{reg.pathmox}}.
#'
#' @param xtree.reg An object of class \code{"xtree.reg"} returned by
#' \code{\link{reg.pathmox}}.
#' @param terminal is string, if equal to \code{TRUE}, just the terminal nodes are considered 
#' for the output reults. when it is equal to \code{FALSE},the regression results are generated 
#' for all nodes of the tree
#' @param intercept if equal to \code{TRUE} also the intercept is considered in the estimation
#' @param label is a boolean. tI is false for defect. If it is \code{TRUE}, label.nodes has to be fix. 
#' @param label.nodes is a vector with the name of the nodes. It is null for defect. 
#' @param \dots Further arguments passed on to \code{\link{reg.treemodel}}. 
#' @return An object of class \code{"regtreemodel"}. Basically a list with the
#' following results:
#' @return \item{inner}{Matrix of the inner relationship between latent variables of the PLS-PM model}
#' @return \item{method}{A string containing the used method ("lm" or "lad"}
#' @return \item{coefficients}{Matrix coefficients for each terminal node}
#' @return \item{Std.}{Matrix of estandard deviation of coefficients for each terminal node}
#' @return \item{pval.coef}{Matrix of p-value significance for each terminal node}
#' @return \item{r2}{Matrix of r-squared coefficients for each terminal node}
#'
#' @author Giuseppe Lamberti
#' 
#' @references Sanchez, G. (2009) \emph{PATHMOX Approach: Segmentation Trees in
#' Partial Least Squares Path Modeling.} PhD Dissertation. 
#' 
#' @references Lamberti, G. (2014) \emph{Modeling with Heterogeneity.} PhD Dissertation. 
#'
#' @seealso \code{\link{pls.pathmox}}, \code{\link{plot.xtree.pls}}
#' @export
#' @examples
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
#' 			signif=0.05,deep=2,method="lm",size=0.15)
#'
#' #terminal nodes comparison
#' fib.node.comp=reg.treemodel(fib.reg.pathmox) 
#'
#'}
		
reg.treemodel 	<- function (xtree.reg,terminal=TRUE,intercept=FALSE,label=FALSE, label.nodes=NULL, ...)
{
    if (class(xtree.reg) != "xtree.reg") 
        stop("Argument 'xtree.reg' must be an object of class 'xtree.reg'")
 		
	x			=	xtree.reg$model$data
	method		=	xtree.reg$model$method
	MOX			=	xtree.reg$MOX
	q			=	ncol(x)
	
	if(terminal==TRUE)	{elementos	=	xtree.reg$terminal}
	if(terminal==FALSE) {elementos	=	xtree.reg$nodes}
	
	coefficients	=	NULL
	d.standard		=	NULL
	r2				=	NULL
	pval.coef		=	NULL
	
	#Calculo de los valores
	
	for ( i in 1 : length(elementos))
	{
		x.node	= 	x[elementos[[i]],]
		resp	=	as.matrix(x.node[,1])
		pred	=	as.matrix(x.node[,-1])
			
		if(method=="lm")
		{
			if(intercept==FALSE)	{reg.node	=	lm(resp~pred-1)}
			if(intercept==TRUE)		{reg.node	=	lm(resp~pred)}
			
			coefficients	=	round(cbind(coefficients,reg.node$coefficients),3)
			pval.coef		=	round(cbind(pval.coef,coef(summary(reg.node))[,4]),3)
			d.standard		=	round(cbind(d.standard,coef(summary(reg.node))[,2]),3)
    		r2[i]			=	1-((sum((resp-reg.node$fitted.values)^2))/(sum((resp-mean(resp))^2)))
		}
		if(method=="lad")
		{	
			if(intercept==FALSE)		{reg.node	=	rq(resp~pred-1,method="fn")}
			if(intercept==TRUE)			{reg.node	=	rq(resp~pred,method="fn")}

			coefficients	=	round(cbind(coefficients,reg.node$coefficients),3)
    	}
	}
	
	#Calculo las etiquetas
	
	if(method=="lm")
	{	
		if(terminal==TRUE) 
		{
			term.nodes <- which(MOX$Terminal == "yes") - 1
			tn.labs <- paste("Node", MOX$Node[term.nodes + 1], sep = "_")
			name.node=c("Root_Node", tn.labs)
			if(label==TRUE){name.node= label.nodes} 
		}
		if(terminal==FALSE) 
		{
			term.nodes <- which(MOX$Depth> 0) - 1
			tn.labs <- paste("Node", MOX$Node[term.nodes + 1], sep = "_")
			name.node=c("Root_Node", tn.labs)
			if(label==TRUE){name.node= label.nodes} 
		}
		colnames(coefficients)=colnames(d.standard)=colnames(pval.coef)=names(r2)=name.node
		
		res.p=list(method=method,coefficients=coefficients,d.standard=d.standard,pval.coef=pval.coef,r2=r2)
		
		
		res=list(method=method,coefficients=coefficients, Std. =d.standard,pval.coef=pval.coef,r2=r2)
	}
	
	if(method=="lad")
	{
		if(terminal==TRUE) 
		{
			term.nodes <- which(MOX$Terminal == "yes") - 1
			tn.labs <- paste("Node", MOX$Node[term.nodes + 1], sep = "_")
			name.node=c("Root_Node", tn.labs)
			if(label==TRUE){name.node= label.nodes} 
		}
				if(terminal==FALSE) 
		{
			term.nodes <- which(MOX$Depth  > 0) - 1
			tn.labs <- paste("Node", MOX$Node[term.nodes + 1], sep = "_")
			name.node=c("Root_Node", tn.labs)
			if(label==TRUE){name.node= label.nodes} 
		}
	
		colnames(coefficients)=name.node

		res=list(method=method,coefficients=coefficients)
	}
	class(res)="treemodelreg"
	res
}
