
#'@title Comparative plot between nodes from or the Pathmox Segmentation Trees: linear and LAD regression
#'
#'@description
#'Plot method for objects of class \code{"treemodelreg"}. Barplots of path
#'coefficients of terminal nodes with respect to those of the global (root)
#'model
#'
#'@details
#'This function aims to visualize the comparison between coefficients of
#'the terminal nodes against the coefficients coefficients of the global model in the
#'root node. 
#'
#'@param x An object of class \code{"treemodelreg"} returned by
#'\code{\link{reg.treemodel}}.
#'@param main.node It is string. If iequl to TRUE you have to inidcate the main of each barplot 
#'in \code{"names.nodes"}.
#'
#'@param names.nodes Optional vector of names for each the terminal node (must be a
#'vector of length equal to the number of terminal nodes).
#'
#' @param  eti  is string. If it is TRUE the label of each coefficients for all the terminal nodes 
#' must be  specify in  \code{"lab.vec"}. If it is false the labels are defined by the programe. 
#'
#' @param lab.vec Optional vector of names for each coefficient of the terminal nodes (must be a
#' vector of length equal to the number of coefficients).
#' @param short.min Integer number indicating the minimum length of the.
#' @param cex.names Allows to fix the size of coefficient labels. Equal to 1 to default.
#' @param cex.axis Allows to fix the size of axes. Equal to 1.2 to default.
#' @param cex.main Allows to fix the size of the main. Equal to 1 to default. 
#' @param lim Allows to fix the axes interval. Equal to (-0.5,0.5) to default.
#' @param short.labs Logical value indicating if the labels of the barplots.
#' @param \dots Further arguments passed on to \code{\link{plot.treemodelreg}}. 
#'
#' \code{\link{reg.treemodel}}, \code{\link{reg.pathmox}}
#'
#' @author Giuseppe Lamberti
#' 
#' @references Sanchez, G. (2009) \emph{PATHMOX Approach: Segmentation Trees in
#' Partial Least Squares Path Modeling.} PhD Dissertation. 
#' 
#' @references Lamberti, G. (2014) \emph{Modeling with Heterogeneity.} PhD Dissertation. 
#'
#' @method plot treemodelreg
#' @S3method plot treemodelreg
#'@examples
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
#'	levels=c("<18k","25k","35k","45k",">45k"), ordered=T)
#' segvar$Accgrade = factor(segvar$Accgrade, 
#'	levels=c("accnote<7","7-8accnote","accnote>8"), ordered=T)
#' segvar$Grade 	= factor(segvar$Grade, 
#'	levels=c("<6.5note","6.5-7note","7-7.5note",">7.5note"), ordered=T)
#'
#  #Regression PATHMOX
#' fib.reg.pathmox=reg.pathmox(Satisfact~.,data=data.fib,segvar,
#'	signif=0.05,deep=2,method="lm",size=0.15)
#'
#' #terminal nodes comparison
#' fib.node.comp=reg.treemodel(fib.reg.pathmox) 
#' 
#  #Drawing the bar-plots
#' plot(fib.node.comp)
#'
#'
#'}


plot.treemodelreg	<- function (x,main.node=FALSE,names.nodes = NULL,eti=FALSE,lab.vec=NULL,short.min = NULL,cex.names=1,cex.axis=1.2,cex.main=1,lim=c(-0.5,0.5),short.labs = TRUE,...) 
{	
	if (class(x) != "treemodelreg") 
        stop("Argument 'x' must be an object of class 'treemodelreg'")

	if (short.labs) 
        if (mode(short.min) != "numeric" || length(short.min) != 1 || (short.min%%1) != 0) 
        short.min <- 5
	
	x=x$coefficients
	
	x.new=NULL
	for(i in 1:(ncol(x)-1))
	{
		x.dif=x[,i+1]-x[,1]
		x.new=cbind(x.new,x.dif)
	}
	y=data.frame(x[,1],x.new)
	
	#Labals terminal nodes

	if(main.node==FALSE){nodes.name=colnames(x)}
	if(main.node==TRUE){nodes.name=names.nodes}
	
	#Labals bares
	
	if(eti==FALSE){lab.name=abbreviate(substring(rownames(x)[1:length(rownames(x))], first=5), minlength = short.min)}
	if(eti==TRUE){lab.name=lab.vec}
	
	rs <- c(1, 1, 1, 2, 2, 2, 2, 2, 3, 2, 3, 3, 4, 4)
	cs <- c(1, 2, 3, 2, 3, 3, 4, 4, 3, 5, 4, 4, 4, 4)
	index.mat <- cbind(1:14, rs, cs)
	lvs <- nrow(x)
	colors = rainbow(lvs, s = 0.5, v = 0.7)
	nn=ncol(x)
	par(mfrow = index.mat[nn + 1, 2:3])
	par(mar = c(3, 3, 3, 3))
	barplot(x[,1], main = paste("Global"), col = colors,cex.names=cex.names,cex.axis=cex.axis,cex.main=cex.main,names.arg=lab.name,ylim=lim)
	abline(h = 0)
	ylim <- 1.15 * c(min(x[1, ]), max(x[1, ]))
	for (n in 2:nn) 
	{
		barplot(x[,n], main = nodes.name[n],col = colors,cex.names=cex.names, cex.axis=cex.axis,cex.main=cex.main,names.arg=lab.name,ylim=lim)
	}
}
