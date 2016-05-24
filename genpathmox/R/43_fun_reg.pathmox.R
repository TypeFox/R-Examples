
#' @title PATHMOX-REG: Segmentation Trees in 
#'  linaer and LAD regression model
#' 
#' @description
#' The function \code{reg.pathmox} calculates a binary segmentation tree in the context 
#' of linear regression following the PATHMOX algorithm. This function generalizes the Pathmox 
#' algorithm introduced by Sanchez in 2009 to the context of linear and LAD regression.
#' 
#' @details
#' The argument \code{formula} is an object of class \code{"formula"} (or one that can be coerced to that class): 
#' a symbolic description of the model to be fitted. 
#'
#'
#' The argument \code{SVAR} must be a data frame containing segmentation
#' variables as factors. The number of rows in
#' \code{SVAR} must be the same as the number of rows in the data 
#'  
#' The argument \code{signif} represent the p-value level takes as reference  
#' to stop the tree partitions. 
#'
#' The argument \code{deep} represent the p-value level takes as reference  
#' to stop the tree partitions. 
#'
#' The argument \code{method} is a string contaning the criterion used to calculate the   
#' the test; if \code{method="lm"} the classic least square approach is used to perform the test; 
#' if \code{method="lad"} the lad (least absolute deviation) is used.
#'
#' The argument \code{size} has defined as a decimal value (i.e. proportion
#' of elements inside a node). 
#'
#' 
#' @param formula An object of class \code{"formula"}.
#' @param SVAR A data frame of factors contaning the segmentation variables.
#' @param signif A numeric value indicating the significance threshold of the
#' F-statistic. Must be a decimal number between 0 and 1.
#' @param deep An integer indicating the depth level of the tree. Must be an
#' integer greater than 1.
#' @param method A string indicating the criterion used to calculate the   
#' the test can be equal to \code{"lm"} or \code{"lad"}
#' node.
#' @param size A numeric value indicating the minimum size of elements inside a
#' node.
#'
#' @param tree A logical value indicating if the tree should be displayed
#' (\code{TRUE} by default).
#' @param data an optional data frame.
#'
#' @param \dots Further arguments passed on to \code{\link{reg.pathmox}}. 
#'
#' @return An object of class \code{"xtree.reg"}. Basically a list with the
#' following results:
#' @return \item{MOX}{Data frame with the results of the segmentation tree}
#' @return \item{root}{element of contaning in the root node}
#' @return \item{terminal}{element of contaning in the terminal nodes}
#' @return \item{nodes}{element of contaning in all nodes terminal and intermediate}
#' @return \item{candidates}{List of data frames containing the candidate 
#' splits of each node partition}
#' @return \item{Fg.r}{Data frame containing the results of the F-global test 
#' for each node partition}
#' @return \item{Fc.r}{A list of Data frames containing the results of the F-coefficients test 
#' for each node partition}
#' @return \item{model}{Information about the internal paramenters} 
#' 
#' @author Giuseppe Lamberti
#'  
#' @references Lamberti, G. (2014) \emph{Modeling with Heterogeneity.} PhD Dissertation. 
#'
#'
#' @references Sanchez, G. (2009) \emph{PATHMOX Approach: Segmentation Trees in
#' Partial Least Squares Path Modeling.} PhD Dissertation. 
#'
#'
#' @export
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
#'  
#'
#'}


reg.pathmox	<- function(formula,SVAR,signif,deep,method,size,tree=TRUE,data=NULL,...)
{
	x <- as.matrix(model.frame(formula,data=data,na.action=NULL))
	
	if (any(is.na(x)))
	{
        warning("Data contains NA: missing values are imputed by the Chained Equations method")
		x=as.matrix(complete(mice(x,m=1,seed=30,printFlag=FALSE)))
    }
    if (!is.data.frame(SVAR))
    {
        stop("Argument 'SVAR' must be a data frame containing factors")
    }
    if (any(is.na(SVAR)))
    {
        warning("data factors contains NA: missing values are imputed by the Chained Equations method")
		SVAR=complete(mice(SVAR,m=1,seed=30,printFlag=FALSE))
    }
    for (j in 1:ncol(SVAR)) if (!is.factor(SVAR[, j]))
    { 
        stop("One or more columns in 'SVAR' are not factors")
	}
	if (nrow(x) != nrow(SVAR))
	{
		stop("Arguments 'x' and 'SVAR' are incompatible. Different number of rows")
    }
	if (mode(signif) != "numeric" || length(signif) != 1 || signif <= 0 || signif >= 1 || is.na(mode(signif)))
	{
        warning("NOTICE: Invalid argument 'signif'. Default value 0.05 was used","\n")
        signif <- 0.05
    }
	if (mode(method) != "character" || method != "lm" && method != "lad" || is.na(mode(method)))
	{
		warning("NOTICE: Invalid argument 'method'. Default method lm was used","\n")
        method <- "lm"  
    }
    if (mode(size) != "numeric" || length(size) != 1 || size <=0 || is.na(mode(size)))
    {
        warning("NOTICE: Invalid argument 'size'. Default value 0.10 was used","\n")
        size <- 0.1
    }
    if (mode(deep) != "numeric" ||  deep < 1 || (deep%%1) != 0 || is.na(mode(deep)))
    {
        warning("NOTICE: Invalid argument 'deep'. Default value 1 was used","\n")
        deep <- 1
    }

	min.ind.node=percent.node(x,size)
	
	info.mox.reg(signif,size,deep,SVAR,method)
	
	model=list(signif=signif,size=size,deep=deep,method=method,data=x,SVAR=SVAR)
	
	#init: create the root
	
	id=0
	t = new ("tree",id=id)
	dim_row = nrow(x)
	elements = seq(1:dim_row)
	id=id+1
	root = new("node.reg",id=id,elements=elements,father=0,childs=0)
	new_nodes = list(root)
	while(length(new_nodes)>0){
		n = new_nodes[[1]]
		
		#init
		
		if(length(n@elements)>=min.ind.node$min.n.ind && showDeepth(n) < deep){
				
		d = x[n@elements,]
		s = SVAR[n@elements,]
		
		cat		=	sapply(s, is.factor)
		s[cat] 	=	lapply(s[cat], factor)
		
		tmp = partopt.reg(d,s,method)
		
		if(tmp$pvl.opt<= signif && any(!is.null(tmp$pvl.opt))==TRUE) {
				
				#before create child nodes
				
				variable		= tmp$variable.opt
				modalidad		= tmp$modalidad.opt
				candidates		= tmp$candidates
				modtwo 			= tmp$modtwo.opt
				
				mod=test.particion.reg(d,modtwo,signif,method)
				
				fglobal = mod$Fg
				fcoef   = mod$Fc
				pvg     = mod$pvg
				pvc     = mod$pvc
				
				#create child nodes
				
				for(i in 1:2){
						elements = n@elements[which(modtwo==i)]
						child_id = (n@id*2)+i-1
						child = new("node.reg",id=child_id,elements=elements,father=n@id)
						n@childs[i] = child_id
						new_nodes[[length(new_nodes)+1]]=child
					}
				n@info=new("info.reg",variable=variable,modalidad=modalidad, fgstatistic=fglobal,fpvalg=pvg,fcstatistic=fcoef,fpvalc=pvc, candidates=candidates )
			
		}}
		t@nodes[[length(t@nodes)+1]]=n
		new_nodes[1]=NULL
	}
	
	if(length(t@nodes)==1)
	{
		root=root.tree(t)
		warning("No sognificative partition faunded")
		res=list(root=root, data=x, method=method)
	}
	else
	{
		MOX=mox.tree(t)
		terminal=terminal.tree(t)
		nodes=nodes.tree(t)
		candidates=candidates.tree(t)
		Fg.r=fglobal.tree(t)
		Fc.r=fcoef.tree.reg(t)
		res=list(MOX=MOX,terminal=terminal,nodes=nodes,candidates=candidates,Fg.r=Fg.r,Fc.r=Fc.r,model=model)

	class(res)="xtree.reg"
	if (tree) 
        plot(res)
    }
res
}
