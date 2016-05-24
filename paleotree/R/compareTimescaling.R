#' Comparing the Time-Scaling of Trees
#' 
#' These functions take two trees and calculate the changes in node ages (for
#' compareNodeAges) for shared clades or terminal branch lengths leading to
#' shared tip taxa (for compareTermBranches).
#' 
#' @details For their most basic usage, these functions compare the time-scaling of two
#' trees. Any taxa not-shared on both trees are dropped before analysis, based
#' on tip labels.
#' 
#' As with many paleotree functions, calculations relating to time on trees are
#' done with respect to any included $root.time elements. If these are not
#' present, the latest tip is assumed to be at the present day (time=0).
#' 
#' compareNodeAges calculates the changes in the clade ages among those clades
#' shared by the two trees, relative to the first tree in absolute time. For
#' example, a shift of +5 means the clade originates 5 time-units later in
#' absolute time on the second tree, while a shift of -5 means the clade
#' originated 5 time-units prior on the second tree.
#' 
#' For compareNodeAges, if tree2 is actually a multiPhylo object composed of
#' multiple phylogenies, the output will be a matrix, with each row
#' representing a different tree and each column a different clade shared
#' between at least some subset of the trees in tree2 and the tree in tree1.
#' values in the matrix are the changes in clade ages between from tree1 (as
#' baseline) to tree2, with NA values representing a clade that isn't contained
#' in the tree represented by that row (but is contained in tree1 and at least
#' one other tree in tree2). The matrix can be reduced to only those clades
#' shared by all trees input via the argument dropUnshared. Note that this
#' function distinguishes clades based on their shared taxa, and cannot infere
#' that two clades might be identical if it were not for single taxon within
#' the crown of one considered clade, despite that such a difference should
#' probably have no effect on compare a node divergence date. Users should
#' consider their dataset for such scenarios prior to application of
#' compareNodeAges, perhaps by dropping all taxa not included in all other
#' trees to be considered (this is NOT done by this function).
#' 
#' compareTermBranches calculates the changes in the terminal branch lengths
#' attached to tip taxa shared by the two trees, relative to the first tree.
#' Thus, a shift of +5 means that this particular terminal taxon is connected
#' to a terminal branch which is five time-units longer.


#' @aliases compareTimescaling compareNodeAges compareTermBranches

#' @param tree1 A time-scaled phylogeny of class 'phylo'

#' @param tree2 A time-scaled phylogeny of class 'phylo'; for compareNodeAges,
#' tree2 can also be an object of class 'multiPhylo' composed of multiple
#' phylogenies. See below.

#' @param dropUnshared If TRUE, nodes not shared across all input trees are
#' dropped from the final output for compareNodeAge. This argument has no
#' effect if tree2 is a single phylogeny (a 'phylo'-class object).

#' @return compareTermBranches returns a vector of temporal shifts for terminal
#' branches with the shared tip names as labels.
#' 
#' compareNodeAges, if both tree1 and tree2 are single trees, outputs a vector
#' of temporal shifts for nodes on tree2 with respect to tree1. If tree2 is
#' multiple trees, then a matrix is output, with each row representing each
#' tree in tree2 (and carrying the name of each tree, if any is given). The
#' values are temporal shifts for each tree in tree2 with respect to tree1. For
#' either case, the column names or element names (for a vector) are the sorted
#' taxon names of the particular clade, the dates of which are given in that
#' column. See above for more details. These names can be very long when large
#' trees are considered.


#' @seealso \code{\link{dateNodes}}, \code{\link{taxa2phylo}}, \code{\link{phyloDiv}}

#' @examples
#' 
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' #get the true tree
#' tree1 <- taxa2phylo(taxa)
#' #simulate a fossil record with imperfect sampling with sampleRanges()
#' rangesCont <- sampleRanges(taxa,r=0.5)
#' #let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
#' cladogram <- taxa2cladogram(taxa,plot=TRUE)
#' #Now let's try timePaleoPhy using the continuous range data
#' tree2 <- timePaleoPhy(cladogram,rangesCont,type="basic")
#' #let's look at the distribution of node shifts
#' hist(compareNodeAges(tree1,tree2))
#' #let's look at the distribution of terminal branch lengths
#' hist(compareTermBranches(tree1,tree2))
#' 
#' #testing ability to compare multiple trees with compareNodeAges
#' trees <- cal3TimePaleoPhy(cladogram,rangesCont,brRate=0.1,extRate=0.1,
#'     sampRate=0.1,ntrees=10)
#' nodeComparison <- compareNodeAges(tree1,trees)
#' #plot it as boxplots for each node
#' boxplot(nodeComparison,names=NULL);abline(h=0)
#' #plot mean shift in node dates
#' abline(h=mean(apply(nodeComparison,2,mean,na.rm=TRUE)),lty=2)
#' 
#' #just shifting a tree back in time
#' set.seed(444)
#' tree1 <- rtree(10)
#' tree2 <- tree1
#' tree1$root.time <- 10
#' compareNodeAges(tree1,tree2)
#' compareTermBranches(tree1,tree2)
#' 
#' @name compareTimescaling
#' @rdname compareTimescaling
#' @export
compareNodeAges<-function(tree1,tree2,dropUnshared=FALSE){
	#output vector of shifts in node dates
		#08-02-12: Allows multiple trees2 to be multiple trees
			#will produce a matrix, each row is a tree in tree2, each column a different but commonly shared clade
	#require(ape)
	if(!inherits(tree1, "phylo")){
		stop("tree1 is not of class phylo")
		}
	tree1orig<-tree1
	if(!inherits(tree2, "phylo")){
		if(!inherits(tree2, "multiPhylo")){
			stop("tree2 is not of class phylo or multiphylo")
			}
		trees2<-tree2
	}else{		#if it isn't multiphylo, make it into one!
		trees2<-list(tree2)
		class(trees2)<-"multiPhylo"
		}
	#okay, need to find all matches common to tree1 and tree2
		#we'll make a MATRIX of all clades held in common between each tree in trees2 and to tree1
		#each row will be a different tree, each column a different clade
		#node dates for each tree will get added as a new column or to an old column
	matchMat<-NULL
	for(i in 1:length(trees2)){
		tree2<-trees2[[i]]
		tree1<-tree1orig
		#incredibly, all of the following is necessary to properly adjust node dates (??!!)
		matches1<-which(!is.na(match(tree1$tip.label,tree2$tip.label)))[1]
		if(length(matches1)<1){stop(paste("No shared taxa between tree2 and tree1[[",i,"]]!"))}
		tipmatch<-tree1$tip.label[matches1]
		mtimeA<-node.depth.edgelength(tree1)[matches1]
		mtimeB<-node.depth.edgelength(tree2)[match(tipmatch,tree2$tip.label)]
		tree1<-drop.tip(tree1,tree1$tip.label[is.na(match(tree1$tip.label,tree2$tip.label))])
		tree2<-drop.tip(tree2,tree2$tip.label[is.na(match(tree2$tip.label,tree1$tip.label))])
		ntime1<-node.depth.edgelength(tree1)
		ntime2<-node.depth.edgelength(tree2)
		mtime1<-ntime1[match(tipmatch,tree1$tip.label)]
		mtime2<-ntime2[match(tipmatch,tree2$tip.label)]
		if(!is.null(tree1$root.time)){
			tree1$root.time<-tree1$root.time-(mtimeA-mtime1)
			ntime1<-tree1$root.time-ntime1
			ntime1<-round(ntime1,6)
			if(min(ntime1)<0){stop(paste("tree1$root.time is less than total depth of tree1!"))}
		}else{
			ntime1<-max(ntime1)-ntime1
			}
		if(!is.null(tree2$root.time)){
			tree2$root.time<-tree2$root.time-(mtimeB-mtime2)
			ntime2<-tree2$root.time-ntime2
			ntime2<-round(ntime2,6)
			if(min(ntime2)<0){stop("tree2[",i,"]$root.time is less than total depth of that tree!")}
		}else{
			ntime2<-max(ntime2)-ntime2
			}
		clades1<-lapply(prop.part(tree1),function(x) sort(tree1$tip.label[x]))
		clades2<-lapply(prop.part(tree2),function(x) sort(tree2$tip.label[x]))
		matches<-match(clades1,clades2)
		cladesMatches<-clades2[matches[!is.na(matches)]]
		if(length(matches[!is.na(matches)])==1){cladesMatches<-list(cladesMatches)}
		ages1<-ntime1[Ntip(tree1)+which(!is.na(matches))]
		ages2<-ntime2[Ntip(tree2)+matches[!is.na(matches)]]
		age_diff<-ages1-ages2
		names(age_diff)<-NULL
		#okay, need to find all matches common to tree1 and tree2
			#we'll make a MATRIX of all clades held in common between each tree in trees2 and to tree1
			#each row will be a different tree, each column a different clade
			#node dates for each tree will get added as a new column or to an old column
		cladesMatches<-sapply(cladesMatches,function(x) paste(x,collapse=","))
		if(is.null(matchMat)){	#if the first tree examined...
			matchMat<-matrix(age_diff,1,)
			colnames(matchMat)<-cladesMatches
		}else{
			currMatches<-match(cladesMatches,colnames(matchMat))
			matchMat<-rbind(matchMat,rep(NA,ncol(matchMat)))
			for(j in 1:length(currMatches)){
				if(!is.na(currMatches[j])){
					matchMat[i,currMatches[j]]<-age_diff[j]
				}else{
					matchMat<-cbind(matchMat,c(rep(NA,i-1),age_diff[j]))
					colnames(matchMat)[ncol(matchMat)]<-cladesMatches[j]
					}
				}
			}
		}
	if(dropUnshared){
		matchMat<-matchMat[,apply(matchMat,2,function(x) all(!is.na(x))),drop=FALSE]
		}
	rownames(matchMat)<-names(trees2)
	if(length(trees2)==1){
		matchMat<-matchMat[1,]
		names(matchMat)<-NULL
		}
	return(matchMat)
	}
	
#' @rdname compareTimescaling
#' @export
compareTermBranches<-function(tree1,tree2){
	#output vector of shifts in terminal branch lengths
	#require(ape)
	if(!inherits(tree1, "phylo")){stop("tree1 is not of class phylo")}
	if(!inherits(tree2, "phylo")){stop("tree2 is not of class phylo")}
	tree1<-drop.tip(tree1,tree1$tip.label[is.na(match(tree1$tip.label,tree2$tip.label))])
	tree2<-drop.tip(tree2,tree2$tip.label[is.na(match(tree2$tip.label,tree1$tip.label))])
	term1<-tree1$edge.length[tree1$edge[,2]<=Ntip(tree1)]
	term1<-term1[order(tree1$edge[tree1$edge[,2]<=Ntip(tree1),2])]
	term2<-tree2$edge.length[tree2$edge[,2]<=Ntip(tree2)]
	term2<-term2[order(tree2$edge[tree2$edge[,2]<=Ntip(tree2),2])]
	term2<-term2[match(tree1$tip.label,tree2$tip.label)]
	term_diff<-term2-term1
	names(term_diff)<-tree1$tip.label
	return(term_diff)
	}