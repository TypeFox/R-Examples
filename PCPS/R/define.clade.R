#' Define clade
#' 
#' Function to define groups (clades) in a phylogenetic tree. 
#' 
#' In the method threshold the total length of phylogenetic tree is used as cutting factor. 
#' If threshold is near to zero the cutting is near the root, if threshold near to one 
#' cutting is near the tips.
#' 
#' The phylogenetic tree must contain the node labels for the function work. Use the 
#' \code{\link{makeNodeLabel}} for defining node labels in a flexible way.
#' 
#' @encoding UTF-8
#' @importFrom phylobase phylo4 descendants
#' @importFrom ape makeNodeLabel node.depth.edgelength
#' @aliases define.clade
#' @param tree Phylogenetic tree.
#' @param threshold A threshold value to form the groups.
#' @param time A cutting height (age) to form the groups.
#' @param method Method to define the clades, "threshold" or "time".
#' @return \item{clades}{Tips and their clades.} \item{height}{The cutting height.}
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{makeNodeLabel}}
#' @keywords PCPS
#' @examples
#' 
#' require(ape)	
#' tree<-makeNodeLabel(rcoal(20))
#' clades<-define.clade(tree, threshold = 0.8, method = "threshold")
#' clades
#' plot.phylo(tree, show.node.label = TRUE)
#' abline(v = clades$height)
#' @export
define.clade<-function(tree,threshold,time,method=c("threshold","time")){
	if(is.null(tree$node.label)){
		stop("\n Node labels not found. Use the function makeNodeLabel\n")
	}
	if(tree$Nnode!=length(tree$node.label)){
		stop("\n Some nodes not found among all nodes in tree\n")
	}
	tree1<-phylobase::phylo4(tree)
	NoDes<-ape::node.depth.edgelength(tree) 
	names(NoDes)=c(tree$tip.label,tree$node.label)
	NoDes<-NoDes[-match(tree$tip.label,names(NoDes))]
	N1<-length(NoDes)
	if(method=="threshold"){
		NoDes2<-NoDes[NoDes/node.depth.edgelength(tree)[1]>=(threshold)]
		height<-node.depth.edgelength(tree)[1]*(threshold)
	}
	else{
		NoDes2<-NoDes[NoDes>=time]
		height<-time
	}
	clades<-vector(length=length(tree$tip.label))	
	names(clades)=tree$tip.label
	clades[]=tree$tip.label	
	N2<-length(NoDes2)
	if(N1==N2){
		N2=N2-1
	}
	n=1
	if(!N2==0){
		for (n in 1:N2){
			RM=n
			Descendants<-phylobase::descendants(tree1,names(sort(NoDes,decreasing=T))[RM])
			NoDe<-names(sort(NoDes,decreasing=T))[RM]
			clades[Descendants]=NoDe
		}
	}
return(list(clades=clades,height=height))
}