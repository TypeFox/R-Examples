
#' @title PLS-PM results of terminal nodes from the Pathmox Segmentation Trees
#' 
#' @description
#' Calculates basic PLS-PM results for the terminal nodes of PATHMOX
#' trees
#' 
#' @details
#' The argument \code{xpls} must be the same used for calculating the
#' \code{xtree} object.  When the object \code{xpls} does not contain a data
#' matrix (i.e. \code{pls$data=NULL}), the user must provide the data matrix or
#' data frame in \code{X}.
#' 
#' The argument \code{xtree} is an object of class \code{"xtree.pls"} returned by 
#'\code{\link{pls.pathmox}}.
#'
#' @param xpls An object of class \code{"plspm"} returned by \code{\link{plspm}}.
#' @param xtree An object of class \code{"xtree.pls"} returned by
#' \code{\link{pls.pathmox}}.
#' @param X Optional dataset (matrix or data frame) used when argument
#' \code{dataset=NULL} inside \code{xpls}.
#' @param alpha is numeric value indicating the significance threshold of the invariance test
#' @param terminal is string, if equal to \code{TRUE}, just the terminal nodes are considered 
#' for the output reults. when it is equal to \code{FALSE},the PLS-PM results are generated 
#' for all nodes of the tree
#' @param scaled to standardize the latent variables or not
#' @param label is a string. It is false for defect. If it is \code{TRUE}, label.nodes has to be fix. 
#' @param label.nodes is a vector with the name of the nodes. It is null for defect. 
#' @param \dots Further arguments passed on to \code{\link{pls.treemodel}}. 
#' @return An object of class \code{"treemodel.pls"}. Basically a list with the
#' following results:
#' @return \item{inner}{Matrix of the inner relationship between latent variables of the PLS-PM model}
#' @return \item{invariance.test}{A data frame containing the results of the invariance test}
#' @return \item{weights}{Matrix of outer weights for each terminal node}
#' @return \item{loadings}{Matrix of loadings for each terminal node}
#' @return \item{paths}{Matrix of path coefficients for each terminal node}
#' @return \item{r2}{Matrix of r-squared coefficients for each terminal node}
#' @return \item{sign}{list of matrix with the significance for each terminal node}
#'
#' @author Giuseppe Lamberti
#' 
#' @references Sanchez, G. (2009) \emph{PATHMOX Approach: Segmentation Trees in
#' Partial Least Squares Path Modeling.} PhD Dissertation. 
#' 
#' @references Lamberti, G. (2014) \emph{Modeling with Heterogeneity.} PhD Dissertation. 
#'
#' @seealso \code{\link{pls.pathmox}},\code{\link{plot.xtree.pls}}
#' @export
#' @examples
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
#'	 Qual.spec		= rep(0,5)
#'	 Qual.gen			= rep(0,5)
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
#'			levels=c("<6.5note","6.5-7note","7-7.5note",">7.5note"), ordered=T)
#'
#'  # Pathmox Analysis
#'  fib.pathmox=pls.pathmox(pls.fib,seg.fib,signif=0.05,
#'					deep=2,size=0.2,n.node=20)
#' 
#'  fib.comp=pls.treemodel(pls.fib,fib.pathmox)
#'  
#'  }
#'

pls.treemodel <- function (xpls,xtree,X=NULL,alpha=0.05,terminal=TRUE,scaled=FALSE, label=FALSE, label.nodes=NULL, ...)
{
  if (class(xpls) != "plspm") 
    stop("Argument 'pls' must be an object of class 'plspm'")
  if (class(xtree) != "xtree.pls") 
    stop("Argument 'xtree' must be an object of class 'xtree.pls'")
  if (nrow(xpls$scores) != xtree$MOX$Size[1]) 
    stop("Arguments 'xpls' and 'xtree' are incompatible. Different number of observations")
  if (!is.null(X)) {
    if (is.null(xpls$data)) {
      if (!is.matrix(X) && !is.data.frame(X)) 
        stop("Invalid object 'X'. Must be a numeric matrix or data frame.")
      if (nrow(X) != nrow(xpls$latents)) 
        stop("Argument 'pls' and 'X' are incompatible. Different number of rows.")
    }
  }
  else {
    if (is.null(xpls$data)) 
      stop("Argument 'X' is missing. No dataset available.")
  }  
    
  x			=	xtree$model$data
  inner   	=	xpls$model$IDM
  outer	  	=	xpls$model$blocks
  mode		=	xpls$model$specs$modes
  scaling 	= 	xpls$model$specs$scaling
  scaled		=	xpls$model$specs$scaled
  scheme		=  	xpls$model$specs$scheme	
  MOX			=   xtree$MOX
  
  
  if(terminal==TRUE)	{nodes	=	xtree$terminal}
  if(terminal==FALSE) {nodes	=	xtree$nodes}
  
  latent		=  	xpls$scores
  
  inv = invariance(x,xtree$terminal,inner,outer,mode,scheme,scaling,scaled)
  
  
  lvs			=	ncol(inner)
  lvs.names	=	row.names(inner)
  path.labs	=	NULL
  for (j in 1:lvs) for (i in j:lvs) if (inner[i, j] == 1) path.labs <- c(path.labs, paste(lvs.names[j], "->", lvs.names[i],sep = ""))
  
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
  
  
  weights=NULL
  loadings=NULL
  paths=NULL
  effect=NULL
  r2=NULL
  significo=list()

  for (i in 1 : length(nodes))
  {
    x.node	= 	x[nodes[[i]],]
    
    pls.node=plspm(x.node,inner,outer,mode,scaling,scheme,scaled=scaled)	
    
    weights		=	round(cbind(weights,pls.node$outer_model[,3]),3)
    loadings	=	round(cbind(loadings,pls.node$outer_model[,4]),3)
    paths		=	round(cbind(paths,pls.node$path_coefs[pls.node$path_coefs!=0]),3)
    effect		=	round(cbind(effect,pls.node$effects[,4]),3)
    rownames(effect)=pls.node$effects[,1]
    r2			=	round(cbind(r2,pls.node$inner_summary$R2),3)
    
    signific=NULL
    signific=list()
    for(k in 1:length(pls.node$inner_mode))	{signific[[length(signific)+1]]=round(cbind(pls.node$inner_model[[k]][,4]),3)}
    significo[[length(significo)+1]]=signific
    
  }
  
  p2=list()
  for(j in 1:length(pls.node$inner_mode)){
    p1=NULL
    for(i in 1:length(nodes)){
      p1=cbind(p1,significo[[i]][[j]])
    }
    colnames(p1)=name.node
    p2[[length(p2)+1]]=p1
  }
  names(p2)=names(pls.node$inner_model)
  
  colnames(weights)=colnames(loadings)=colnames(paths)=colnames(r2)=colnames(effect)=name.node
  rownames(paths)=path.labs
  rownames(weights)=pls.node$outer_model[,1]
  rownames(loadings)=pls.node$outer_model[,1]
  rownames(r2)=rownames(pls.node$inner_summary)
  
  if(inv$p.value>alpha){
    weights=cbind(weights,inv$avg.weights)
  }
  
  res=list(IDM=inner,invariance.test=inv$test,weights=weights,loadings=loadings,paths=paths,r2=r2,sign=p2)
  class(res)="treemodel"
  res
}
