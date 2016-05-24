#' @title PATHMOX-PLS: Extended Segmentation Trees in 
#' Partial Least Squares Path Modeling 
#' 
#' @description
#' The function \code{pathmox.pls} calculates a binary segmentation tree in the context 
#' PLS-PM following the PATHMOX algorithm. This function extends the pathmox algorithm 
#' introduced by Sanchez in 2009 including the two new test: the F-block test (to detect
#' the responsible latent endogenous equations of the difference), the F-coefficient test
#' (to detect the path coefficients responsible of the difference).The F-tests used in the 
#' split process are implemented following the classic lest square estimation. An implementation 
#' of the tests following the LAD regression also are proposed to overcome the parametric 
#' hypothesis of the F-test.
#' 
#' @details
#' The argument \code{xpls} is object of class \code{"plspm"} returned by \code{\link{plspm}}.
#'
#' The argument \code{SVAR} must be a data frame containing segmentation
#' variables as factors. The number of rows in
#' \code{SVAR} must be the same as the number of rows in the data used in
#' \code{pls}. 
#'  
#' The argument \code{signif} represent the p-value level takes as reference  
#' to stop the tree partitions. 
#'
#' The argument \code{deep} represent the depth level of the tree takes as reference  
#' to stop the tree partitions. 
#'
#' The argument \code{method} is a string contaning the criterion used to calculate the   
#' tests; if \code{method="lm"} the classic least square approach is used to perform the tests; 
#' if \code{method="lad"} the LAD (least absolute deviation regression) is used.
#'
#' The argument \code{size} is defined as a decimal value (i.e. proportion
#' of elements inside a node). 
#'
#' The argument \code{n.node} is the minimum number of individuals to consider a candidate 
#' partition. If the candidate split produces a partition where the number of individuals is less  
#' then \code{n.node}, the partition is not considered.
#'
#' When the object \code{pls} does not contain a data matrix (i.e.
#' \code{pls$data=NULL}), the user must provide the data matrix or data frame in
#' \code{X}.
#' 
#' @param xpls An object of class \code{"plspm"} returned by \code{\link{plspm}}.
#' @param SVAR A data frame of factors contaning the segmentation variables.
#' @param signif A numeric value indicating the significance threshold of the
#' F-statistic. Must be a decimal number between 0 and 1.
#' @param deep An integer indicating the depth level of the tree. Must be an
#' integer greater than 1.
#' @param method A string indicating the criterion used to calculate the   
#' the test can be equal to \code{"lm"} or \code{"lad"}.
#' @param size A numeric value indicating the minimum size of elements inside a
#' node.
#' @param X Optional dataset (matrix or data frame) used when argument
#' \code{dataset=NULL} inside \code{pls}.
#' @param tree A logical value indicating if the tree should be displayed
#' (\code{TRUE} by default).
#' @param n.node It is the minimum number of individuals to consider a candidate 
#' partition (\code{30} by default).
#' 
#'
#' @param \dots Further arguments passed on to \code{\link{pls.pathmox}}. 
#'
#' @return An object of class \code{"xtree.pls"}. Basically a list with the
#' following results:
#' @return \item{MOX}{Data frame with the results of the segmentation tree}
#' @return \item{root}{List of elements contanined in the root node}
#' @return \item{terminal}{List of elements contanined in terminal nodes}
#' @return \item{nodes}{List of elements contanined in all nodes: terminal and intermediate}
#' @return \item{candidates}{List of data frames containing the candidate 
#' splits of each node partition}
#' @return \item{Fg.r}{Data frame containing the results of the F-global test 
#' for each node partition}
#' @return \item{Fb.r}{List of data frames containing the results of the F-block test 
#' for each node partition}
#' @return \item{Fc.r}{A list of data frames containing the results of the F-coefficients test 
#' for each node partition}
#' @return \item{model}{Informations about the internal paramenters} 
#' 
#' @author Giuseppe Lamberti
#' 
#'
#' @references Lamberti, G. (2014) \emph{Modeling with Heterogeneity.} PhD Dissertation. 
#'
#' @references Sanchez, G. (2009) \emph{PATHMOX Approach: Segmentation Trees in
#' Partial Least Squares Path Modeling.} PhD Dissertation. 
#' 
#'
#' @export
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
#'  }
#'
pls.pathmox = function(xpls,SVAR,signif,deep,method="lm",size,X = NULL,tree = TRUE,n.node=30,...)
{
    if (class(xpls) != "plspm") 
        stop("Argument 'pls' must be an object of class 'plspm'")
    
    if (!is.null(X)) {
        if (is.null(xpls$data)) {
            if (!is.matrix(X) && !is.data.frame(X)) 
                stop("Invalid object 'X'. Must be a numeric matrix or data frame.")
            if (nrow(X) != nrow(xpls$latents)) 
                stop("Argument 'pls' and 'X' are incompatible. Different number of rows.")
            if (nrow(X) != nrow(SVAR)) 
                stop("Arguments 'X' and 'SVAR' are incompatible. Different number of rows")
        }
    }
    else {
        if (is.null(xpls$data)) {
            stop("Argument 'X' is missing. No dataset available.")
        }
        else {
            if (nrow(xpls$data) != nrow(SVAR)) 
                stop("Arguments 'pls' and 'SVAR' are incompatible. Different number of rows")
        }
    }
    if (!is.data.frame(SVAR)) 
        stop("Argument 'SVAR' must be a data frame containing factors")
    for (j in 1:ncol(SVAR)) if (!is.factor(SVAR[, j])) 
        stop("One or more columns in 'SVAR' are not factors")
	
	x		=	xpls$data
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
    if (mode(signif) != "numeric" || length(signif) != 1 || signif <= 
        0 || signif >= 1) {
        warning("NOTICE: Invalid argument 'signif'. Default value 0.05 was used", 
            "\n")
        signif <- 0.05
    }
	if (mode(method) != "character" || method != "lm" && method != "lad" ){
        warning("NOTICE: Invalid argument 'method'. Default method lm was used", 
            "\n")
        method <- "lm"  
    }
    if (mode(size) != "numeric" || length(size) != 1 || size <= 
        0) {
        warning("NOTICE: Invalid argument 'size'. Default value 0.10 was used", 
            "\n")
        size <- 0.1
    }
    if (mode(deep) != "numeric" ||  deep < 
        1 || (deep%%1) != 0) {
        warning("NOTICE: Invalid argument 'deep'. Default value 1 was used", 
            "\n")
        deep <- 1
    }
	
	info.mox.pls(signif,size,deep,SVAR)

	model	= list(signif=signif,size=size,deep=deep,method=method,data=x,SVAR=SVAR)
	
	inner   	=	xpls$model$IDM
	outer	  	=	xpls$model$blocks
	mode		=	xpls$model$specs$modes
	scaling 	= 	xpls$model$specs$scaling
	scaled		=	xpls$model$specs$scaled
	scheme		=  	xpls$model$specs$scheme	
	
	min.ind.node = percent.node(x,size)
	
	#init: create the root
	id=0
	t = new ("tree",id=id)
	dim_row = xpls$model$gens$obs
	elements = seq(1:dim_row)
	id = id+1
	root = new("node",id=id,elements=elements,father=0,childs=0)
	new_nodes = list(root)
	while(length(new_nodes) > 0){
		n = new_nodes[[1]]
		
		#init
		
		if(length(n@elements)>=min.ind.node$min.n.ind && showDeepth(n) < deep){
				
		d = x[n@elements,]
		s = SVAR[n@elements,]
		cat		=	sapply(s, is.factor)
		s[cat] 	=	lapply(s[cat],factor)
		
		tmp = partopt.pls(d,s,inner,outer,mode,scheme,scaling,scaled,method,n.node)

		if(tmp$pvl.opt <= signif && any(!is.null(tmp$pvl.opt)) == TRUE) {
				
				#before create child nodes
				
				variable		= tmp$variable.opt
				level		= tmp$level.opt
				candidates	= tmp$candidates
				modtwo		= tmp$modtwo.opt
				
				mod=test.partition.pls(d,inner,outer,mode,scheme,scaling,scaled,modtwo,signif,method)
				
				fglobal = mod$Fg
				fblock  = mod$Fb
				fcoef   = mod$Fc
				pvg     = mod$pvg
				pvb     = mod$pvb
				pvc     = mod$pvc
				
				#create child nodes
				
				for(i in 1:2){
						elements = n@elements[which(modtwo==i)]
						child_id = (n@id*2)+i-1
						child = new("node",id=child_id,elements=elements,father=n@id)
						n@childs[i] = child_id
						new_nodes[[length(new_nodes)+1]]=child
					}
				n@info=new("info.pls",
							variable=variable,
							level=level,
							fgstatistic=fglobal,
							fpvalg=pvg,
							fbstatistic=fblock,
							fpvalb=pvb,
							fcstatistic=fcoef,
							fpvalc=pvc, 
							candidates=candidates)
			
		}}
		t@nodes[[length(t@nodes)+1]] = n
		new_nodes[1] = NULL
	}
	if	(length(t@nodes)==1)
	{
		root = root.tree(t)
		cat("No sognificative partition faunded")
		res = list(root=root,model=model)
	}
	else
	{
		MOX = mox.tree.pls(t)
		root = root.tree(t)
		terminal = terminal.tree(t)
		nodes = nodes.tree(t)
		candidates = candidates.tree(t)
		Fg.r = fglobal.tree.pls(t)
		Fb.r = fblock.tree.pls(t)
		Fc.r = fcoef.tree.pls(t)
		res = list(MOX=MOX,
				root=root,
				terminal=terminal,
				nodes=nodes,
				candidates=candidates,
				Fg.r=Fg.r,
				Fb.r=Fb.r,
				Fc.r=Fc.r,
				model=model)
	class(res) = "xtree.pls"
	if (tree) 
        plot(res)
	}
 res
}
