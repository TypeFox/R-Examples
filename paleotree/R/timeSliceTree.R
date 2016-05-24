#' Time-Slicing a Phylogeny
#' 
#' Removes the portion of a tree after a set point in time, as if the tree
#' after that moment had been sliced away.
#' 
#' @details The function assumes that ttree will generally have an element called
#' $root.time, which is the time before present that the root divergence
#' occurred. If $root.time is not present as an element of ttree, then it is
#' assumed the tip furthest from the root is at time 0 (present-day) and a new
#' $root.time is calculated (a warning will be issued in this case).
#' 
#' The sliceTime is always calculated as on the same scale as ttree$root.time.
#' In other words, if root.time=100, then timeSlice=80 will slice the tree 20
#' time units after the root.
#' 
#' If drop.extinct=T, then extinct tips are dropped and (if present) the
#' $root.time of ttree is adjusted. This is done using the function
#' dropExtinct.
#' 
#' @param ttree A time-scaled phylogeny of class phylo
#' @param sliceTime Time to 'slice' the tree at. See details below.
#' @param drop.extinct If true, drops tips that go extinct before timeSlice.
#' @param plot If true, plots input and output trees for comparison.
#' @return Returns the modified phylogeny as an object of class phylo.
#' 
#' Tip labels for cut branches which held multiple descendant tips will be the
#' label for the earliest appearing tip descendant of that branch. This is
#' somewhat arbitrary; the actual morphotaxon present at that time might have
#' been a different taxon. For simulated datasets, use taxa2phylo.

#' @author David W. Bapst, with modification of code by Klaus Schliep to avoid use of
#' function \code{dist.nodes} which has difficulty with large trees and greatly
#' benefiting function runtime.

#' @seealso \code{\link{phyloDiv}}, \code{\link{dropExtinct}},
#' \code{\link{dropExtant}}
#' 
#' Also see the function treeSlice in the library phytools, which will slice a
#' tree at some point in and return all the subtrees which remain after the
#' slicing time. (Effectively the reverse of timeSliceTree.)
#' @examples
#' 
#' #a neat example of using phyloDiv with timeSliceTree 
#'    #to simulate doing extant-only phylogeny studies 
#'    #of diversification...in the past!
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' taxicDivCont(taxa)
#' #that's the whole diversity curve
#' 
#' #with timeSliceTree we could look at the lineage accumulation curve 
#'    #we'd get of species sampled at a point in time
#' tree <- taxa2phylo(taxa)
#' #use timeSliceTree to make tree of relationships up until time=950 
#' tree950 <- timeSliceTree(tree,sliceTime=950,plot=TRUE,drop.extinct=FALSE)
#' 
#' #use drop.extinct=T to only get the tree of lineages extant at time=950
#' tree950 <- timeSliceTree(tree,sliceTime=950,plot=FALSE,drop.extinct=TRUE)
#' #now its an ultrametric tree with many fewer tips...
#' #lets plot the lineage accumulation plot on a log scale
#' phyloDiv(tree950,plotLogRich=TRUE)
#' 
#' @export timeSliceTree
timeSliceTree<-function(ttree,sliceTime,drop.extinct=FALSE,plot=TRUE){
	#take a phylogeny and produce a phylogenetic 'slice' at time X (w/respect to root.time)
		#lineages extant at sliceTime sliced to end at that point
		#if no root.time, then it is assumed the tip furthest from the root is at 0 (present-day)
			#a warning will be issued if this is assumed
		#extinct tips will be dropped if drop.extinct=TRUE
	#require(ape)
	if(!inherits(ttree, "phylo")){
		stop("ttree is not of class phylo")
		}
	if(is.null(ttree$root.time)){
		ttree$root.time<-max(node.depth.edgelength(ttree)[1:Ntip(ttree)])
		message("Warning: no ttree$root.time! Assuming latest tip is at present (time=0)")
		}
	tslice<-ttree$root.time-sliceTime	#time from root to slice time
	#first let's drop all edges that branch later than the slice
	#make them all single lineages by dropping all but one taxon
	dnode<-node.depth.edgelength(ttree)
	#identify the ancestor nodes of edges which cross the tslice
	cedge<-which((dnode[ ttree$edge[, 1] ] < tslice) & (dnode[ttree$edge[, 2] ]  >= tslice))
	droppers<-numeric()
	propPartTree<-prop.part(ttree)
	for(i in 1:length(cedge)){
		desc<-ttree$edge[cedge[i],2]
		if(desc>Ntip(ttree)){	#if an internal edge that goes past the tslice
			desctip<-propPartTree[[desc-Ntip(ttree)]]	#drop all but one tip
			droppers<-c(droppers,desctip[-1])
		}}
	stree<-drop.tip(ttree,droppers)
	#which edges cross over tslice?
	dnode<-node.depth.edgelength(stree)
	cedge<-(dnode[stree$edge[, 2] ]  >= tslice)
	cnode_depth<-dnode[stree$edge[cedge,1]]
	stree$edge.length[cedge]<-tslice-cnode_depth
	stree$root.time<-ttree$root.time		#root.time shouldn't (can't?) change!
	if(drop.extinct){
		stree1<-dropExtinct(stree,ignore.root.time=TRUE)
	}else{
		stree1<-stree
		}
	if(plot){
		layout(1:2)
		plot(ladderize(ttree),show.tip.label=FALSE)
			axisPhylo()
		plot(ladderize(stree1),show.tip.label=FALSE)
			axisPhylo()
		layout(1)
		}
	return(stree1)
	}
