#' Scales Edge Lengths of a Phylogeny to a Minimum Branch Length
#' 
#' Rescales a tree with edge lengths so that all edge lengths are at least some minimum branch length (mbl),
#' without changing the relative distance of the tips from the root node. Edge lengths are transformed so they are
#' greater than or equal to the input minimum branch length, by subtracting edge length from more rootward edges
#' and added to later branches. 

#' @details
#' This function was formally an internal segment in \code{\link{timePaleoPhy}}, and now is called by \code{timePaleoPhy}
#' instead, allowing users to apply \code{minBranchLength} to trees that already have edge lengths.

#' @param tree A phylogeny with edge lengths of class 'phylo'.

#' @param mbl The minimum branch length 

#' @return
#' A phylogeny with edge lengths of class 'phylo'.

#' @seealso
#' This function was originally an internal piece of \code{\link{timePaleoPhy}}, which implements the minimum branch
#' length time-scaling method along with others, which may be what you're looking for
#' (instead of this miscellaneous function).

#' @author 
#' David W. Bapst

#' @examples
#' 
#' #simulation with an example non-ultrametric tree
#' 
#' tree<-rtree(20)
#' # randomly replace edges with ZLBs, similar to multi2di output
#' tree<-degradeTree(tree,0.3,leave.zlb=TRUE) 	
#' 
#' tree2<-minBranchLength(tree,0.1)
#' 
#' layout(1:2)
#' plot(tree);axisPhylo()
#' plot(tree2);axisPhylo()
#' 
#' layout(1)
#' 
#' #now let's try it with an ultrametric case
#' 
#' # get a random tree
#' tree<-rtree(30)
#' # randomly replace edges with ZLBs, similar to multi2di output
#' tree<-degradeTree(tree,0.5,leave.zlb=TRUE) 
#' # now randomly resolve	
#' tree<-di2multi(tree)
#' # give branch lengths so its ultrametric
#' tree<-compute.brlen(tree)
#' 
#' plot(tree) #and we have an ultrametric tree with polytomies, yay!
#' 
#' #now randomly resolve, get new branch lengths as would with real data
#' tree2<-multi2di(tree)
#' tree2<-minBranchLength(tree2,0.1)
#' 
#' layout(1:2)
#' plot(tree,show.tip.label=FALSE);axisPhylo()
#' plot(tree2,show.tip.label=FALSE);axisPhylo()
#' 
#' layout(1)

#' @name minBranchLength
#' @rdname minBranchLength
#' @export
minBranchLength<-function(tree, mbl){	
	#require(phangorn)
	#test arguments
	#tree - a tree with edge lengths
	if(!inherits(tree, "phylo")){
		stop("tree is not of class phylo")
		}
	if(is.null(tree$edge.length)){stop("Tree has no edge lengths")}
	timetree<-tree
	#mbl - a single numeric value
	if(!is.numeric(mbl) | length(mbl)!=1){
		stop("mbl is a not a single numeric value")}
	#
	#root_node<-Ntip(timetree)+1
	while(any(timetree$edge.length<mbl)){
		# get every edge that is short
		shortEdge<-(1:Nedge(timetree))[timetree$edge.length<mbl]
		# pick the shortest one, if multiple of that length, pick first one
		shortLength<-timetree$edge.length[shortEdge]
		shortestLength<-shortEdge[shortLength==min(shortLength)]
		mom<-timetree$edge[shortestLength[1],1]
		#make vector of every mom node that is ancestral
		mom<-c(mom,Ancestors(timetree,mom))
		debt<-mbl-min(timetree$edge.length[timetree$edge[,1]==mom[1]])
		timetree$edge.length[mom[1]==timetree$edge[,1]]<-timetree$edge.length[mom[1]==timetree$edge[,1]] + debt[1]
		#make vector of smallest brlen with each mom node as anc
		#calculate, simulatenously, the changes in debt and branch lengthening required as go down tree
		#change branch lengths; hypothetically, debt should then equal zero...
		if(length(mom)>1){for(i in 2:length(mom)){
			small<-min(timetree$edge.length[timetree$edge[,1]==mom[i]])
			mom_blen<-timetree$edge.length[timetree$edge[,1]==mom[i] & timetree$edge[,2]==mom[i-1]]
			debt[i]<-max(debt[i-1] - max(mom_blen-mbl,0),0) + max(mbl-small,0) 
			timetree$edge.length[timetree$edge[,1]==mom[i] & timetree$edge[,2]==mom[i-1]] <- 
			mom_blen - max(min(max(mom_blen-mbl,0),debt[i-1]),0) + max(mbl-small,0)
			timetree$edge.length[timetree$edge[,1]==mom[i] & timetree$edge[,2]!=mom[i-1]] <-  
			timetree$edge.length[timetree$edge[,1]==mom[i] & timetree$edge[,2]!=mom[i-1]] + debt[i]
			}}
		}
	return(timetree)
	}

	