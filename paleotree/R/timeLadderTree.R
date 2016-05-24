#' Resolve Polytomies by Order of First Appearance
#' 
#' Resolves polytomies in trees with lineages arranged in a pectinate pattern
#' (i.e. a ladder-like subtree), ordered by the time of first appearance (FAD)
#' for each lineage.
#' 
#' @details
#' This method of resolving polytomies assumes that the order of stratigraphic
#' appearance perfectly depicts the order of branching. This may not be a good
#' assumption for poorly sampled fossil records.
#' 
#' This function is for resolving trees when a continuous time-scale is known.
#' For discrete time-scales, see the function bin_timePaleoPhy.
#' 
#' Taxa with the same identical first appearance date will be ordered randomly.
#' Thus, the output is slightly stochastic, but only when ties exist. This is
#' probably uncommon with real data on continuous time-scales.
#' 
#' Taxa not shared between the tree and the timeData matrix, or listed as
#' having a FAD or LAD of NA in timeData will be dropped and will not be
#' included in the output tree.
#' 
#' See this blog post for more information:
#'
#' http://nemagraptus.blogspot.com/2012/07/resolving-polytomies-according-to.html
#' 

#' @param tree A phylo object

#' @param timeData Two-column matrix of per-taxon first and last occurrances in
#' absolute continous time

#' @return Returns the modified tree as an object of class phylo, with no edge
#' lengths.

#' @author David W. Bapst

#' @seealso \code{\link{di2multi}}

#' @examples
#' 
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(100,200))
#' taxa<-fossilRecord2fossilTaxa(record)
#' tree<-taxa2cladogram(taxa)
#' ranges<-sampleRanges(taxa,r=0.5)
#' tree1<-timeLadderTree(tree,ranges)
#' layout(1:2)
#' plot(ladderize(tree),show.tip.label=FALSE)
#' plot(ladderize(tree1),show.tip.label=FALSE)
#' 
#' #an example with applying timeLadderTree to discrete time data
#' rangeData<-binTimeData(ranges,int.len=5)	#sim discrete range data
#' tree2<-bin_timePaleoPhy(tree,timeList=rangeData,timeres=TRUE)
#' plot(ladderize(tree),show.tip.label=FALSE)
#' plot(ladderize(tree2),show.tip.label=FALSE)
#' axisPhylo() 
#'
#' layout(1)
#' 
#' @export timeLadderTree
timeLadderTree<-function(tree,timeData){
	#resolves all polytomies in a tree as ladders to match FADs in timeData
	#only applicable to continuous time data
	#require(ape)	
	#first sanitize data
	if(!inherits(tree,"phylo")){
		stop("tree must be of class 'phylo'")
		}
	if(!inherits(timeData,"matrix")){
		if(inherits(timeData,"data.frame")){
			timeData<-as.matrix(timeData)
		}else{
			stop("timeData not of matrix or data.frame format")
			}
		}
	if(ncol(timeData)==6){	#also allow it to accept taxad objects
		timeData<-timeData[,3:4,drop=FALSE]
		}	
	#first clean out all taxa which are NA or missing in timeData
	#remove taxa that are NA or missing in timeData
	tree<-drop.tip(tree,tree$tip.label[is.na(match(tree$tip.label,names(which(!is.na(timeData[,1])))))])
	if(Ntip(tree)<2){
		stop("Less than two valid taxa shared between the tree and temporal data")
		}
	timeData<-timeData[!is.na(timeData[,1]),]
	if(any(is.na(timeData))){
		stop("Weird NAs in Data??")
		}
	if(any(timeData[,1]<timeData[,2])){
		stop("timeData is not in time relative to modern (decreasing to present)")
		}
	#node IDs of polytomies
	nodes<-Ntip(tree)+1:Nnode(tree)
	polys<-nodes[sapply(nodes,function(x) sum(tree$edge[,1]==x)>2)]
	#error if no polytomies
	if(length(polys)==0){stop("No Polytomies in Cropped Tree?")}
	while(length(polys)>0){
		node<-polys[1]
		#make a list of the first FAD of each descendant
		desc<-tree$edge[tree$edge[,1]==node,2]
		dTips<-lapply(desc,function(x) if(x>Ntip(tree)){
			tree$tip.label[prop.part(tree)[[x-Ntip(tree)]]]
			}else{tree$tip.label[x]})
		#get max FAD
		dFADs<-sapply(dTips,function(x) max(timeData[x,1]))
		#build a ladderized subtree of the node
		subtree<-stree(length(dFADs),type="left")
		#figure out the order of tips, use a second vector of random numbers for ties
		subtree$tip.label<-desc[order(-dFADs,sample(1:length(dFADs)))]
		#stick on tip labels
		for(i in desc){
			dtip<-which(subtree$tip.label==i)
			if(i>Ntip(tree)){		#if its a clade
				subclade<-extract.clade(tree,i)
				subtree<-bind.tree(subtree,subclade,where=dtip)
				subtree<-collapse.singles(subtree)
				}else{subtree$tip.label[dtip]<-tree$tip.label[i]}	#if its a tip
			}
		#replace original node with new, resolved, scaled node
		if(node!=(Ntip(tree)+1)){	#if it isn't the root node
			drtips<-prop.part(tree)[[node-Ntip(tree)]]
			tip_lab<-tree$tip.label[drtips[1]]	#cut out all but one tip, just to put it back together later
			droptree<-collapse.singles(drop.tip(tree,drtips[-1]))
			droptree<-bind.tree(droptree,subtree,where=which(droptree$tip.label==tip_lab))	#put in subtree at tip
			tree<-droptree
		}else{				#if it is the root node
			tree<-subtree
			}	
		nodes<-Ntip(tree)+1:Nnode(tree)
		polys<-nodes[sapply(nodes,function(x) sum(tree$edge[,1]==x)>2)]
		}
	tree$edge.length<-NULL
	return(tree)
	}
