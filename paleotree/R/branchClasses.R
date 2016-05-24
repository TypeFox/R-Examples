#' Partitions the branch lengths of a tree into several classes based on their placement.
#'
#' @details This function will partition the internode (node to node, including internal node to terminal tip) branch lengths of a tree into 
#' four separate classes: all (all the internode branches of a tree), int (internal branches
#' which run from one internode to another), live (terminal branches which run from an internal node to 
#' a terminal tip representing an extinction event before the present) and dead (terminal branches 
#' which run from an internal node to a terminal tip at the modern day, reflecting a still-living taxon).
#'
#' The depths of the internal 'mother' node (i.e. time of origin, before the modern day 
#' of each branch length are included as the /code{names} of the branch length vectors.
#'
#' This function is mainly of use for modeling internode branch lengths in a phylogeny including fossil taxa.

#' @param tree A time-scaled phylogeny to be analysed, as an object of class phylo, ideally with a $root.time element as is
#' typical for paleotree output phylogenies. If $root.time is not present, the most recent tips will be interpreted as 
#' being at the modern day (i.e. 0 time-units before present).

#' @param whichExtant A logical vector with length equal to number of tips in the tree. TRUE indicates a taxon that is extant
#' at the moden day. If present 

#' @param tol Tolerance used to distinguish extant taxa, if whichExtant is not provided, to avoid issues with number rounding. Taxa within tol of 
#' the modern day will be considered extant.

#' @return The output is a list consisting of four vectors, with the /code{names} of the vectors being their
#' corresponding time of origin. See details.

#' @examples
#' #simulated example
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=c(10,20))
#' taxa<-fossilRecord2fossilTaxa(record)
#' tree <- taxa2phylo(taxa)
#' brlenRes <- branchClasses(tree)
#'
#'#see frequency histograms of branch lengths
#' layout(1:4)
#' for(x in 1:length(brlenRes)){ 
#' 	hist(brlenRes[[x]],main="Branch Lengths",xlab=names(brlenRes)[x])
#' 	}
#'
#'#see frequency histograms of branch depths
#' layout(1:4)
#' for(x in 1:length(brlenRes)){ 
#' 	hist(as.numeric(names(brlenRes[[x]])),main="Branch Depths",xlab=names(brlenRes)[x])
#' 	}
#'
#' layout(1)

#' @export
branchClasses<-function(tree,whichExtant=NULL,tol=0.01){
	#names will be time of origin for each branch
	#require(ape)
	if(!inherits(tree,"phylo")){stop("tree is not of class phylo")}
	dists<-node.depth.edgelength(tree)
	if(is.null(whichExtant)){
		dists1<-dists[1:Ntip(tree)]
		if(is.null(tree$root.time)){modern<-max(dists1)}else{modern<-tree$root.time}
		modern.tips<-which(dists1>(modern-tol))
		dead.tips<-which(dists1<=(modern-tol))
	}else{
		if(length(whichExtant)!=Ntip(tree)){stop("Length of whichExtant is not equal to number of tip taxa on tree")}
		modern.tips<-which(sapply(tree$tip.label,function(x) any(x==names(which(whichExtant==1)))))
		dead.tips<-which(sapply(tree$tip.label,function(x) any(x==names(which(whichExtant==0)))))
		}
	int.edges<-tree$edge[,2]>Ntip(tree)
	live.edges<-sapply(tree$edge[,2],function(x) any(x==modern.tips))
	dead.edges<-sapply(tree$edge[,2],function(x) any(x==dead.tips))
	if(is.null(tree$root.time)){
		depths<-max(dists)-dists[tree$edge[,1]]
		}else{
		depths<-tree$root.time-dists[tree$edge[,1]]
		}
	brlen.all<-tree$edge.length
		names(brlen.all)<-depths
	brlen.int<-tree$edge.length[int.edges]
		names(brlen.int)<-depths[int.edges]
	brlen.live<-tree$edge.length[live.edges]
		names(brlen.live)<-depths[live.edges]		
	brlen.dead<-tree$edge.length[dead.edges]
		names(brlen.dead)<-depths[dead.edges]
	list(brlen.all=brlen.all,brlen.int=brlen.int,brlen.live=brlen.live,brlen.dead=brlen.dead)
	}