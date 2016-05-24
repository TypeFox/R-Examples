#' Modify, Drop or Bind Terminal Branches of Various Types (Mainly for Paleontological Phylogenies)
#' 
#' These functions modify terminal branches or drop certain terminal branches
#' based on various criteria.

#' DropZLB drops tip-taxa that are attached to the tree via zero-length
#' terminal branches ("ZLBs"). This is sometimes useful for paleo-trees, as
#' various time-scaling methods often produce these ZLBs, taxa whose early
#' appearance causes them to be functionally interpreted as ancestors in some
#' time-scaling methods. Removing ZLBs is advised for analyses of
#' diversification/diversity, as these will appear as simultaneous
#' speciation/extinction events. Note this function only drops tips attached to
#' a terminal zero-length branch; if you want to collapse internal zero-length
#' branches, see the ape function \code{\link{di2multi}}.
#' 
#' DropExtinct drops all terminal branches which end before the modern (i.e.
#' extinct taxa). DropExtant drops all terminal branches which end at the
#' modern (i.e. extant/still-living taxa). In both cases, the modern is defined
#' based on tree$root.time if available, or the modern is inferred to be the
#' point in time when the tip furthest from the root (the latest tip)
#' terminates.
#' 
#' If the input tree has a $root.time element, as expected for most paleo-tree
#' objects handled by this library, that root.time is adjusted if the relative
#' time of the root divergence changes when terminal branches are dropped.
#' Adjusted root.times are only given if the input tree has root.times.
#' 
#' addTermBranchLength adds an amount equal to the argument 'addtime' to the
#' terminal branch lengths of the tree. If there is a $root.time element, this
#' is increased by an amount equal to addtime. A negative amount can be input
#' to reduce the length of terminal branches. However, if negative branch
#' lengths are produced, the function fails and a warning is produced. This function
#' does *not* call fixRootTime, so the root.time elements in the result tree may
#' be nonsensical, particularly if negative amounts are input.
#' 
#' When a tree is modified, such as having tips dropped or branches extended,
#' fixRootTime can be used to find the new $root.time. It is mainly used as a
#' utility function called by a majority of the other functions discussed
#' in this help file.
#'
#' dropPaleoTip is a wrapper for ape's \code{drop.tip} which also modifies the
#' $root.time element if necessary, using \code{fixRootTime}. Similarly,
#' bindPaleoTip is a wrapper for phytool's \code{bind.tip} which allows tip age
#' as input and modifies the $root.time element if necessary (i.e. if a tip
#' is added to edge leading up to the root).
#'
#' Note that for bindPaleoTip, tips added below the root are subtracted from
#' any existing $root.edge element, as per behavior of bind.tip and bind.tree.
#' However, bindPaleoTip will append a $root.edge of the appropriate length
#' if one does not exist (or is not long enough) to avoid an error. After
#' binding is finished, any $root.edge equal to 0 is removed before the
#' resulting tree is output.
#' 

#' @aliases dropZLB dropExtinct dropExtant addTermBranchLength fixRootTime

#' @param tree A phylogeny as a phylo object. dropPaleoTip requires this
#' input object to also have a $root.time element. If not provided for
#' bindPaleoTip, then the root.time will be presumed to be such that the
#' furthest tip from the root is at time=0.

#' @param tol Tolerance for determining modern age; used for distinguishing
#' extinct from extant taxa. Tips which end within 'tol' of the furthest
#' distance from the root will be treated as 'extant' taxa for the purpose of
#' keeping or dropping.

#' @param ignore.root.time Ignore root.time in calculating which tips are
#' extinct? root.time will still be adjusted

#' @param addtime Extra amount of time to add to all terminal branch lengths.

#' @param treeOrig A phylo object of a time-scaled phylogeny with a $root.time
#' element

#' @param treeNew A phylo object containing a modified form of treeOrig (with
#' no extra tip taxa added, but possibly with some tip taxa removed).

#' @param consistentDepth A logical, either TRUE or FALSE. If TRUE (the default)
#' the tree's root-to-furthest-tip depth is tested to make sure this depth is
#' not greater than the new root.time appended to the output tree.

#' @param nodeAgeTransfer A logical. If TRUE, the root.time of the new tree is determined by
#' comparing the clades of taxa between the two input trees. The new root age assigned is the age of
#' (\emph{1}) the treeOrig clade that contains *all* taxa present in treeNew and, if the set of (1)
#' contains multiple clades, (\emph{2}) the clade in the (1) set that contains the fewest taxa not in
#' treeNew. If FALSE, the root.time assigned to treeNew is the root.time of treeOrig, adjusted
#' based on the change in total tree depth between treeOrig and treeNew, as measured between the root and
#' the first matching taxon in both trees. The later is how fixRootTime functioned by default
#' prior to paleotree v2.3.

#' @param ... additional arguments passed to dropPaleoTip are passed to drop.tip

#' @param tipLabel A character string of \code{length = 1} containing the name of the new tip
#' to be added to \code{tree}.

#' @param tipAge The age of the tip taxon added to the tree, in time before present (i.e. where
#' present is 0), given in the same units as the edges of the tree are already scaled. Cannot be
#' given if \code{edgeLength} is given.

#' @param edgeLength The new \code{edge.length} of the terminal branch this tip is connected to.
#' Cannot be given if \code{tipAge} is given. 

#' @param nodeAttach Node or tip ID number (as given in tree$edge) at which to attach the new tip. 
#' See documentation of \code{bind.tip} for more details.

#' @param positionBelow The distance along the edge below the node to be attached to 
#' (given in code{nodeAttach} to add the new tip. Cannot be negative or greater than the length of the
#' edge below \code{nodeAttach}.

#' @param noNegativeEdgeLength Return an error if a negative terminal edge length is calculated
#' for the new tip.

#' @return Gives back a modified phylogeny as a phylo object, generally with a
#' modified $root.time element.

#' @author David W. Bapst. The functions dropTipPaleo and bindTipPaleo are modified imports of
#' \code{drop.tip} and \code{bind.tip} from packages \code{ape} and \code{phytools}.

#' @seealso 
#' \code{\link{compareTermBranches}}, \code{\link{phyloDiv}}, 
#' \code{\link{drop.tip}}, \code{\link{bind.tip}}

#' @examples
#' 
#' set.seed(444)
#' #Simulate some fossil ranges with simFossilRecord
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' #simulate a fossil record with imperfect sampling with sampleRanges
#' rangesCont<-sampleRanges(taxa,r=0.5)
#' #Now let's make a tree using taxa2phylo
#' tree <- taxa2phylo(taxa,obs_time=rangesCont[,2])
#' #compare the two trees
#' layout(1:2)
#' plot(ladderize(tree))
#' plot(ladderize(dropZLB(tree)))
#' layout(1)
#' 
#' #example using dropExtinct and dropExtant
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=c(10,20))
#' taxa<-fossilRecord2fossilTaxa(record)
#' tree<-taxa2phylo(taxa)
#' phyloDiv(tree)
#' tree1 <- dropExtinct(tree)
#' phyloDiv(tree1)
#' tree2 <- dropExtant(tree)
#' phyloDiv(tree2)
#'
#' #graphics.off()
#' 
#' #example using addTermBranchLength
#' set.seed(444)
#' treeA <- rtree(10)
#' treeB <- addTermBranchLength(treeA,1)
#' compareTermBranches(treeA,treeB)
#' 
#' #########################
#' #test dropPaleoTip
#' 	#(and fixRootTime by extension...)
#'
#' #simple example
#' tree<-read.tree(text="(A:3,(B:2,(C:5,D:3):2):3);")
#' tree$root.time<-10
#' plot(tree,no.margin=FALSE)
#' axisPhylo()
#'
#' # now a series of tests, dropping various tips
#' (test<-dropPaleoTip(tree,"A")$root.time) # =7
#' (test[2]<-dropPaleoTip(tree,"B")$root.time) # =10
#' (test[3]<-dropPaleoTip(tree,"C")$root.time) # =10
#' (test[4]<-dropPaleoTip(tree,"D")$root.time) # =10
#' (test[5]<-dropPaleoTip(tree,c("A","B"))$root.time) # =5
#' (test[6]<-dropPaleoTip(tree,c("B","C"))$root.time) # =10
#' (test[7]<-dropPaleoTip(tree,c("A","C"))$root.time) # =7
#' (test[8]<-dropPaleoTip(tree,c("A","D"))$root.time) # =7
#' 
#' # is it all good? if not, fail so paleotree fails...
#' if(!identical(test,c(7,10,10,10,5,10,7,7))){stop("fixRootTime fails!")}
#'
#' 
#' ##############
#' #testing bindPaleoTip
#' 
#' # simple example 
#' tree<-read.tree(text="(A:3,(B:2,(C:5,D:3):2):3);")
#' tree$root.time<-20
#' plot(tree,no.margin=FALSE)
#' axisPhylo()
#' 
#' \dontrun{
#' 
#' require(phytools)
#'
#' #bindPaleoTip effectively wraps bind.tip from phytools
#' # using a conversion like below
#' 
#' tipAge<-5
#' node<-6
#' 
#' #new length = the root time - tipAge - nodeheight(tree,node)
#'
#' newLength<-tree$root.time-tipAge-nodeheight(tree,node)
#' tree1<-bind.tip(tree,"tip.label",where=node,edge.length=newLength)
#' 
#' layout(1:2)
#' plot(tree);axisPhylo()
#' plot(tree1);axisPhylo()
#'
#' }
#'
#' # now with bindPaleoTip
#' 
#' 
#' tree1<-bindPaleoTip(tree,"new",nodeAttach=6,tipAge=5)
#' 
#' layout(1:2)
#' plot(tree);axisPhylo()
#' plot(tree1);axisPhylo()
#' layout(1)
#' 
#' #then the tip age of "new" should 5
#' test<-dateNodes(tree1)[which(tree1$tip.label=="new")]==5
#' if(!test){stop("bindPaleoTip fails!")}
#' 
#' # with positionBelow
#' 
#' tree1<-bindPaleoTip(tree,"new",nodeAttach=6,tipAge=5,positionBelow=1)
#' 
#' layout(1:2)
#' plot(tree);axisPhylo()
#' plot(tree1);axisPhylo()
#' layout(1)
#' 
#' # at the root
#' 
#' tree1<-bindPaleoTip(tree,"new",nodeAttach=5,tipAge=5)
#' 
#' layout(1:2)
#' plot(tree);axisPhylo()
#' plot(tree1);axisPhylo()
#' layout(1)
#' 
#' #then the tip age of "new" should 5
#' test<-dateNodes(tree1)[which(tree1$tip.label=="new")]==5
#' if(!test){stop("bindPaleoTip fails!")}
#' 
#' # at the root with positionBelow
#' 
#' tree1<-bindPaleoTip(tree,"new",nodeAttach=5,tipAge=5,
#' 	positionBelow=3)
#' 
#' layout(1:2)
#' plot(tree);axisPhylo()
#' plot(tree1);axisPhylo()
#' layout(1)
#' 
#' #then the tip age of "new" should 5
#' test<-dateNodes(tree1)[which(tree1$tip.label=="new")]==5
#' #and the root age should be 23
#' test1<-tree1$root.time==23
#' if(!test | !test1){stop("bindPaleoTip fails!")}
#' 


#' @name modifyTerminalBranches
#' @rdname modifyTerminalBranches
#' @export
dropZLB<-function(tree){
	#drops terminal branches that are zero length
		#adjusts tree$root.time if necessary
	#require(ape)
	#checks
	if(!inherits(tree,"phylo")){
		stop("tree must be of class 'phylo'")
		}
	drop_e<-(tree$edge[,2]<(Ntip(tree)+1)) & (tree$edge.length==0)
	drop_t<-(tree$edge[,2])[drop_e]
	if((Ntip(tree)-length(drop_t))>1){
		tree1<-drop.tip(tree,drop_t)
		if(!is.null(tree$root.time)){tree1<-fixRootTime(tree,tree1)}
		res<-tree1
	}else{res<-NA}
	return(res)
	}
	
#' @rdname modifyTerminalBranches
#' @export
dropExtinct<-function(tree,tol=0.01,ignore.root.time=FALSE){
	#drop all terminal taxa that are less than 0.001 from the modern
	#require(ape)
	#checks
	if(!inherits(tree,"phylo")){
		stop("tree must be of class 'phylo'")
		}
	if(is.null(tree$root.time)){
		message("No tree$root.time: Assuming latest tip is at present (time=0)")
		}
	dnode<-node.depth.edgelength(tree)[1:Ntip(tree)]
	dnode<-round(dnode,6)
	if(!is.null(tree$root.time) & !ignore.root.time){
		if(round(tree$root.time,6)>max(dnode)){
		stop("all tips are extinct based on tree$root.time!")}}
	droppers<-which((dnode+tol)<max(dnode))
	if((Ntip(tree)-length(droppers))<2){stop("Less than 2 tips are extant on the tree!")}
	stree<-drop.tip(tree,droppers)
	if(!is.null(tree$root.time)){
		#now need to add $root.time given the droppers
		#should be root.time MINUS distance from furthest tip in tree PLUS distance from latest tip to root of stree
		#stree$root.time<-tree$root.time-max(dnode)+max(dist.nodes(stree)[1:Ntip(stree),Ntip(stree)+1])
		stree<-fixRootTime(tree,stree)
		}
	return(stree)
	}
	
#' @rdname modifyTerminalBranches
#' @export
dropExtant<-function(tree,tol=0.01){
	#drop all terminal taxa that are more than 0.001 from the modern
	#require(ape)
	#checks
	if(!inherits(tree,"phylo")){
		stop("tree must be of class 'phylo'")
		}
	if(is.null(tree$root.time)){
		message("Warning: no tree$root.time! Assuming latest tip is at present (time=0)")
		}
	dnode<-node.depth.edgelength(tree)[1:Ntip(tree)]
	dnode<-round(dnode,6)
	if(!is.null(tree$root.time)){if(round(tree$root.time,6)>max(dnode)){stop("all tips are extinct based on tree$root.time!")}}
	droppers<-which((dnode+tol)>max(dnode))
	if((Ntip(tree)-length(droppers))<2){stop("Less than 2 tips extinct on the tree!")}
	stree<-drop.tip(tree,droppers)
	if(!is.null(tree$root.time)){
		#now need to add $root.time given the droppers
		#should be root.time MINUS distance from earliest tip in tree PLUS distance from earliest tip to root of stree
		#stree$root.time<-tree$root.time-min(dnode)+min(dist.nodes(stree)[1:Ntip(stree),Ntip(stree)+1])
		stree<-fixRootTime(tree,stree)
		}
	return(stree)
	}
	
#' @rdname modifyTerminalBranches
#' @export
addTermBranchLength<-function(tree,addtime=0.001){
	#require(ape)
	#checks
	if(!inherits(tree,"phylo")){
		stop("tree must be of class 'phylo'")
		}

	tree$edge.length[tree$edge[,2]<(Ntip(tree)+1)]<-tree$edge.length[tree$edge[,2]<(Ntip(tree)+1)]+addtime
	if(any(tree$edge.length<0)){stop("tree has negative branch lengths!")}
	if(!is.null(tree$root.time)){tree$root.time<-tree$root.time+addtime}
	return(tree)
	}
	
#' @rdname modifyTerminalBranches
#' @export
fixRootTime<-function(treeOrig,treeNew,consistentDepth=TRUE,nodeAgeTransfer=TRUE){
	# UNUSED FUNCTION?
	#treeDepth<-function(tree){
	#	#require(ape)
	#	max(dist.nodes(tree)[,Ntip(tree)+1])
	#	}
	#
	#require(ape)
	if(!inherits(treeOrig, "phylo")){
		stop("treeOrig is not of class phylo")}
	if(!inherits(treeNew, "phylo")){
		stop("treeNew is not of class phylo")}
	if(is.null(treeOrig$root.time)){
		stop("ERROR: treeOrig passed to fixRootTime with no $root.time??")}
	#also need a warning message if taxa present in treeNew that aren't in treeOrig
	taxaNewNM<-treeNew$tip.label[sapply(treeNew$tip.label,function(x) !any(x==treeOrig$tip.label))]
	if(length(taxaNewNM)>0){
		stop(paste("taxa:",taxaNewNM,"are present in treeNew but not treeOrig"))}
	#two different ways to fix the root time
	if(nodeAgeTransfer){
		##NEW WAY 11-28-14
			#If TRUE, the root.time of the new tree is determined by
			#comparing the clades of taxa between the two input trees. The new root age assigned is the age of
			#(\emph{1}) the treeOrig clade that contains *all* taxa present in treeNew and, if the set of (1)
			#contains multiple clades, (\emph{2}) the clade in the (1) set that contains the fewest taxa not in
			#treeNew.
		dates<-dateNodes(treeOrig,labelDates=FALSE,tolerance=0.001)
		treeDesc<-lapply(Descendants(treeOrig),function(x) sort(treeOrig$tip.label[x]))
		#treeRootNew<-sort(treeNew$tip.label[Descendants(treeNew)[[Ntip(treeNew)+1]]]) #no
		#the descendants of treeNew's root are ALL the taxa in treeNEW
		#So which treeOrig clades contain ALL taxa in treeNew?
		allNewTaxa<-sapply(treeDesc,function(x) all(sapply(treeNew$tip.label,function(y) any(y==x)))) #logical vector
		#now, if more than one contains all-new-taxa, which of these treeOrig clades minimizes not-shared taxa?
		if(sum(allNewTaxa)>1){
			nUnshared<-sapply(treeDesc,function(x) sum(sapply(x,function(y) all(y!=treeNew$tip.label)))) #numeric
			matchRootNew<-which(allNewTaxa & nUnshared==min(nUnshared[allNewTaxa]))
		}else{
			matchRootNew<-which(allNewTaxa)
			}
		if(length(matchRootNew)>1){stop("More than one node contains these taxa")} #maybe sort by age
		if(length(matchRootNew)<1){stop("No nodes match the new tree's root, a root age can not be obtained")}
		treeNew$root.time<-unname(dates[matchRootNew])
	}else{
		##OLD WAY
			#If FALSE, the root.time assigned to treeNew is the root.time of treeOrig, adjusted
			# based on the change in total tree depth between treeOrig and treeNew, as measured between the root and
			# the first matching taxon in both trees. The later is how fixRootTime functioned by default
			# prior to paleotree v2.3.
		orig_dist<-node.depth.edgelength(treeOrig)[which(treeNew$tip.label[1]==treeOrig$tip.label)]
		new_dist<-node.depth.edgelength(treeNew)[1]
		treeNew$root.time<-treeOrig$root.time-(orig_dist-new_dist)
		}
	if(consistentDepth){
		if(round(max(node.depth.edgelength(treeNew)) - treeNew$root.time)>0){
			stop("fixRootTime isn't fixing correctly, root.time less than max tip-to-root length!")}
		}
	return(treeNew)
	}
	
#' @rdname modifyTerminalBranches
#' @export
dropPaleoTip<-function(tree, ...){
	tree1<-drop.tip(phy=tree, ...)
	tree2<-fixRootTime(tree,tree1)
	return(tree2)
	}
	
#' @rdname modifyTerminalBranches
#' @export
bindPaleoTip<-function(tree, tipLabel, nodeAttach=NULL, tipAge=NULL,
		edgeLength=NULL, positionBelow=0, noNegativeEdgeLength=TRUE){
	# CHECKS
	tipLabel<-as.character(tipLabel)
	if(!is.character(tipLabel)){
		stop("cannot coerce tipLabel to a string value")
		}
	if(length(tipLabel)!=1){
		stop("A string of length=1 is needed for tipLabel (i.e. for a single tip) is required")
		}
	# positionBelow
	if(positionBelow<0){
		stop("bindTipPaleo does not accept negative positionBelow values")
		}
	if(nodeAttach==(Ntip(tree)+1)){
		if(positionBelow>0){
			if(is.null(tree$root.edge)){
				tree$root.edge<-positionBelow
			}else{
				if(tree$root.edge<positionBelow){
					tree$root.edge<-positionBelow
					}
				}	
			}
	}else{
		if(positionBelow>tree$edge.length[tree$edge[,2]==nodeAttach]){
			stop("positionBelow cannot be greater than the $edge.length of the edge below nodeAttach")
			}
		}
	# check root.time
	if(is.null(tree$root.time)){
		message("Warning: no tree$root.time! Setting root.time such that latest tip is at present (time=0)")
		tree$root.time<-max(node.depth.edgelength(tree))
		}
	#
	if(is.null(tree$edge.length)){stop("bindTipPaleo is for trees with edge lengths")}
	#
	if(is.null(edgeLength)){
		if(!is.null(tipAge)){
			#nodeHeight<-nodeheight(tree,nodeAttach)-position
			nodeHeight<-node.depth.edgelength(tree)[nodeAttach]
			modNodeHeight<-nodeHeight-positionBelow
			newLength<-tree$root.time-tipAge-modNodeHeight
			if(newLength<0){
				if(noNegativeEdgeLength){
					stop(paste0("Negative edge length created due to tipAge being",
						" older than age of nodeAttach + positionBelow"))
				}else{
					message(paste0("Warning: negative edge length created due to tipAge being",
						" older than age of nodeAttach + positionBelow"))
				}
			}
		}else{
			stop("either tipAge or edgeLength must be given")
			}
	}else{
		if(!is.null(tipAge)){
			stop("both tipAge or edgeLength cannot be given")
			}
		newLength<-edgeLength
		if(newLength<0 & noNegativeEdgeLength){
			stop("Negative edgeLength given ?!")
			}
		}
	tree1<-bind.tip(tree, tip.label=tipLabel, where=nodeAttach,
		position=positionBelow, edge.length=newLength)
	# fix root.time if nodeAttach=root ID of tree and positionBelow>0
	if(nodeAttach==(Ntip(tree)+1) & positionBelow>0){
		#adjust root.time by the positionBelow
		tree1$root.time<-tree1$root.time+positionBelow
		}
	if(!is.null(tree1$root.edge)){
		if(tree1$root.edge==0){tree1$root.edge<-NULL}
		}
	#return tree1
	return(tree1)
	}

