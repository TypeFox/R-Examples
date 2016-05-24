#' Paint Tree Branch Depth by Color
#' 
#' Paints the edges of a phylogeny with colors relative to their depth.
#' 
#' @details The only purpose of this function is to make an aesthetically-pleasing
#' graphic of one's tree, where branches are color-coded with a rainbow
#' palette, relative to their depth. Depth is defined relative to the number of
#' branching nodes between the basal node of a branch and the root, not the
#' absolute distance (i.e. branch length) to the root or the distance from the
#' tips.
#' 
#' @param tree A phylo object
#' @return No value returned, just plots a colorful phylogeny.
#' @examples
#' 
#' set.seed(444)
#' tree <- rtree(500)
#' depthRainbow(tree)
#' 
#' @export depthRainbow
depthRainbow<-function(tree){
	#plots a tree with edges color-coded to depth
	if(!inherits(tree, "phylo")){
		stop("tree isn't of class phylo")
		}
	tree<-ladderize(tree)
	ndepth<-node.depth.edgelength(tree)
	#nodelabels(ceiling(ndepth[(Ntip(tree):Nedge(tree))+1]),node=(Ntip(tree):Nedge(tree))+1)
	edepth<-ceiling((ndepth[(Ntip(tree):Nedge(tree))+1])[tree$edge[,1]-Ntip(tree)])+1
	col_edge<-rainbow(max(edepth))[edepth]
	plot(ladderize(tree),show.tip.label=FALSE,edge.color=col_edge);axisPhylo()
	}
