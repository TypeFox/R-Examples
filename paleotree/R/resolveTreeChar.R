#' Resolve Polytomies Using Parsimony-Based Reconstruction of a Discrete Character
#'
#' This function resolves a set of given topology with less than fully-binary phylogenetic resolution so that
#' lineages are shifted and internal nodes added that minimize the number of independent character transitions needed to explain
#' an observed distribution of discrete character states for the taxa on such a tree, under various maximum-parsimony algorithms of 
#' ancestral character reconstruction, powered ultimately by function \code{ancestral.pars} in library \code{phangorn}.
#' This function is mainly designed for use with poorly resolved trees which are being assessed with the function
#' \code{\link{minCharChange}}.

#' @details
#' As shown in the example code below, this function offers a wide variety of options for manipulating the
#' maximum-parsimony algorithm used (i.e. MPR versus ACCTRAN), the ordering (or not) of character states,
#' and potential biasing of uncertainty character state reconstructions (when ordered characters are
#' assessed). This allows for a wide variety of possible resolutions for a given tree with polytomies
#' and a discrete character. In general, the author expects that use of this function will be optimal
#' when applied to ordered characters using one of the \code{stateBias} options, perhaps
#' \code{stateBias = "primitive"} (based on theoretical expectations for slow evolving characters). However,
#' anecdotal use of this function with various simulation datasets suggests that the results are quite
#' variable, and so the best option needs to be assessed based on the prior assumptions regarding the
#' data and the performance of the dataset with the various arguments of this function.
#'

#' @inheritParams minCharChange

#' @param orderedChar Is the character of interest given for \code{trait} ordered or not?
#' If FALSE (the default), then for each polytomy, all child nodes that appear to have the
#' same state as the ancestor node will remain in the polytomy, and any additional states held by child nodes will
#' each be grouped into their own unique polytomy that forms from a descendant node of the original polytomy.
#' If TRUE, then the character will be reconstructed with a cost (step)
#' matrix of a linear, ordered character, and polytomies will be resolved so that lineages with different
#' states will be placed into a nested ladder that reflects the ordered character. As with the unordered option,
#' child nodes with a state equivalent to the ancestral node will remain in the polytomy, while more primitive
#' or more derived states will be sorted into their own separate ladders composed of paraphyletic groups, ordered
#' so to move 'away' state-by-state from the ancestral node's inferred character state.
#' This option is not applicable if \code{type = "ACCTRAN"}, as cost matrices cannot
#' be used with ACCTRAN in \code{ancestral.pars}, and an error will be returned if \code{orderedChar=TRUE} but
#' a cost matrix is given manually.

#' @param stateBias This argument controls how \code{resolveTreeChar} handles ancestral node reconstructions that have
#' multiple states competing for the maximum weight of any state (i.e. if states 0 and 1 both have 0.4 of the weight). The
#' default, where \code{stateBias = NULL} causes uncertainty at nodes among states to be treated as a single 'group' identical
#' to any states within it. Essentially, this means that for the example polytomy where the ancestor hax maximum weight for both 0 and 1, 
#' any child nodes with 0, 1 or both of these states will be considered to have an identical state for the purpose of grouping nodes
#' for the purpose of further resolving polytomies. If and only if \code{orderedChar = TRUE}, then additional options of
#' \code{stateBias = 'primitive'} and \code{stateBias = 'derived'} become available, which instead force uncertain node
#' assignments to either be the most primitive (i.e. the minimum) or the most derived (i.e. the maximum) among the
#' maximum-weight states. In particular, \code{stateBias = 'primitive'} should favor gains and bias any analysis of
#' character transitions against finding reversals.
#'
#' @param iterative A logical argument which, if TRUE (the default), causes the function to repeat the polytomy-resolving
#' functionality across the entire tree until the number of nodes stabilizes. If FALSE, polytomies are only passed a single
#' time.

#' @return
#' Returns the resulting tree, which may be fully resolved, partly more resolved or not more resolved at all
#' (i.e. have less polytomies) depending on what was possible, as constrained by ambiguities in character
#' reconstructions. Applying \code{\link{multi2di}} is suggested as a post-step to obtain a fully-resolved
#' cladogram, if one is desired.

#' @seealso 
#' \code{\link{ancPropStateMat}} which is used internally by this function. This function was
#' intentionally designed for use with \code{\link{minCharChange}}.

#' @author David W. Bapst

#' @references
#' Hanazawa, M., H. Narushima, and N. Minaka. 1995. Generating most parsimonious reconstructions on
#' a tree: A generalization of the Farris-Swofford-Maddison method. Discrete Applied Mathematics
#' 56(2-3):245-265.
#' 
#' Narushima, H., and M. Hanazawa. 1997. A more efficient algorithm for MPR problems in phylogeny.
#' Discrete Applied Mathematics 80(2-3):231-238.
#'
#' Schliep, K. P. 2011. phangorn: phylogenetic analysis in R. \emph{Bioinformatics} 27(4):592-593.
#'
#' Swofford, D. L., and W. P. Maddison. 1987. Reconstructing ancestral character states under
#' Wagner parsimony. Mathematical Biosciences 87(2):199-229.

#' @examples
#' 
#' \donttest{
#' 
#' # let's write a quick&dirty ancestral trait plotting function
#' 
#'  quickAncPlot<-function(tree,trait,cex,orderedChar=FALSE,type="MPR",cost=NULL){
#' 	ancData<-ancPropStateMat(tree=tree, trait=trait, orderedChar=orderedChar)
#' 	ancCol<-(1:ncol(ancData))+1
#'  	plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="upwards")
#'  	tiplabels(pch=16,pie=ancData[(1:Ntip(tree)),],cex=cex,piecol=ancCol,
#' 		col=0)
#'  	nodelabels(pie=ancData[-(1:Ntip(tree)),],cex=cex,piecol=ancCol)	
#'  	}
#' 
#' ##########
#' 
#' # examples with simulated data
#' 
#' set.seed(2)
#' tree<-rtree(50)
#' #simulate under a likelihood model
#' trait<-rTraitDisc(tree,k=3,rate=0.7)
#' tree<-degradeTree(tree,prop_collapse=0.6)
#' tree<-ladderize(tree,right=FALSE)
#' 
#' #a bunch of type=MPR (default) examples
#' treeUnord<-resolveTreeChar(tree,trait,orderedChar=FALSE)
#' treeOrd<-resolveTreeChar(tree,trait,orderedChar=TRUE,stateBias=NULL)
#' treeOrdPrim<-resolveTreeChar(tree,trait,orderedChar=TRUE,stateBias="primitive")
#' treeOrdDer<-resolveTreeChar(tree,trait,orderedChar=TRUE,stateBias="derived")
#' #and finally an unordered one with ACCTRAN
#' treeACCTRAN<-resolveTreeChar(tree,trait,orderedChar=FALSE,type="ACCTRAN")
#' 
#' #compare number of nodes
#' Nnode(tree)			#original
#' Nnode(treeUnord)		#unordered, biasStates=NULL, MPR
#' Nnode(treeOrd)		#ordered, biasStates=NULL
#' Nnode(treeOrdPrim)	#ordered, biasStates='primitive'
#' Nnode(treeOrdDer)	#ordered, biasStates='derived'
#' Nnode(treeACCTRAN)	#unordered, biasStates=NULL, ACCTRAN
#' 
#' #let's compare original tree with unordered-resolved tree
#' layout(1:2)
#' quickAncPlot(tree,trait,orderedChar=FALSE,cex=0.3)
#' text(x=43,y=10,"Original",cex=1.5)
#' quickAncPlot(treeUnord,trait,orderedChar=FALSE,cex=0.3)
#' text(x=43,y=10,"orderedChar=FALSE",cex=1.5)
#' #some resolution gained
#' 
#' #now let's compare the original and ordered, both biasStates=NULL
#' layout(1:2)
#' quickAncPlot(tree,trait,orderedChar=FALSE,cex=0.3)
#' text(x=43,y=10,"Original",cex=1.5)
#' quickAncPlot(treeOrd,trait,orderedChar=TRUE,cex=0.3)
#' text(x=43,y=10,"orderedChar=TRUE",cex=1.5)
#' 
#' #now let's compare the three ordered trees
#' layout(1:3)
#' quickAncPlot(treeOrd,trait,orderedChar=TRUE,cex=0.3)
#' text(x=41,y=8,"ordered, biasStates=NULL",cex=1.5)
#' quickAncPlot(treeOrdPrim,trait,orderedChar=TRUE,cex=0.3)
#' text(x=41.5,y=8,"ordered, biasStates='primitive'",cex=1.5)
#' quickAncPlot(treeOrdDer,trait,orderedChar=TRUE,cex=0.3)
#' text(x=42,y=8,"ordered, biasStates='derived'",cex=1.5)
#' 
#' #let's compare unordered with ordered, biasStates='primitive'
#' layout(1:2)
#' quickAncPlot(treeUnord,trait,orderedChar=FALSE,cex=0.3)
#' text(x=41,y=8,"orderedChar=FALSE",cex=1.5)
#' quickAncPlot(treeOrdPrim,trait,orderedChar=TRUE,cex=0.3)
#' text(x=40,y=11,"orderedChar=TRUE",cex=1.5)
#' text(x=40,y=4,"biasStates='primitive'",cex=1.5)
#' 
#' #let's compare unordered with MPR to unordered with ACCTRAN
#' layout(1:2)
#' quickAncPlot(treeUnord,trait,orderedChar=FALSE,type="MPR",cex=0.3)
#' text(x=41,y=8,"unordered: MPR",cex=1.5)
#' quickAncPlot(treeACCTRAN,trait,orderedChar=FALSE,type="ACCTRAN",cex=0.3)
#' text(x=41,y=8,"unordered: ACCTRAN",cex=1.5)
#' 
#' #these comparisons will differ greatly between datasets
#' 	# need to try them on your own
#' 
#' layout(1)
#' 
#' }
#'

#' @name resolveTreeChar
#' @rdname resolveTreeChar
#' @export
resolveTreeChar<-function(tree, trait, orderedChar=FALSE, stateBias=NULL, iterative=TRUE, type="MPR", cost=NULL,
			ambiguity= c(NA, "?"), dropAmbiguity=FALSE, polySymbol="&", contrast=NULL){
		#	orderedChar=TRUE; type="MPR"; cost=NULL; stateBias="primitive"
		#	orderedChar=FALSE; type="MPR"; cost=NULL; stateBias=NULL
	#orderedChar=TRUE : put clades together relative to character being ordered, 0 is most primitive state
		# nested paraphyletic grades
	#if FALSE, clades within a polytomy are formed for all groups which the ancestor has 0 prob of being
	#require(phangorn)
	if(Nnode(tree)==(Ntip(tree)-1)){stop("Input tree is fully resolved??")}
	#tests of orderedChar and stateBias
	if(!is.null(stateBias)){
		if(orderedChar){
			if(!any(stateBias==c("primitive","derived"))){
				stop("For 'orderedChar' analyses, stateBias must be set to NULL, 'primitive' or 'derived'")}
		}else{
			stop("stateBias cannot be used in analyses for unordered characters")
			}}
	#check iterative
	if(!is.logical(iterative)){
		stop("iterative must be a logical class element")
		}
	if(length(iterative)!=1){
		stop("iterative must be a single logical element")
		}
	#now run the resolution mechanism function, possibly in a loop
	if(iterative){
		tree2<-tree
		continueRes<-TRUE
		while(continueRes){
			tree1<-resolveTreeCharMechanism(tree2, trait, orderedChar=orderedChar, stateBias=stateBias, type=type, cost=cost)
			if(is.binary.tree(tree1) & is.rooted(tree1)){continueRes<-FALSE}
			if(Nnode(tree1)==Nnode(tree2)){continueRes<-FALSE}
			tree2<-tree1
			}
		treeFinal<-tree2
	}else{
		#only do it once
		treeFinal<-resolveTreeCharMechanism(tree, trait, orderedChar=orderedChar, stateBias=stateBias, type=type, cost=cost)
		}
	return(treeFinal)
	}

# internal function called by resolveTreeChar
resolveTreeCharMechanism<-function(tree, trait, orderedChar, stateBias, type, cost){
	nodeParts<-c(1:Ntip(tree),prop.part(tree))
	#pick nodes based on being a polytomy
	whichPoly<-sapply(1:max(tree$edge),function(x) sum(x==tree$edge[,1])>2)
	#polytomy parts
	polyParts<-nodeParts[whichPoly & node.depth(tree)>1]
	#sort polyParts by depth
	polyParts<-polyParts[order(node.depth(tree)[whichPoly & node.depth(tree)>1])]	
		#depthPoly<-cbind(1:max(tree$edge),node.depth(tree))[whichPoly & node.depth(tree)>1,]
		#chosen<-(depthPoly[order(depthPoly[,2]),]
	#replace polyParts number IDs with taxon labels
	polyParts<-lapply(polyParts,function(x) sort(tree$tip.label[x]))
	tree1<-tree
	for(i in 1:length(polyParts)){	
		#check to make sure taxon lists are still the same
		if(!identical(sort(tree$tip.label),sort(tree1$tip.label))){
			stop("taxa mysteriously lost/added?")}
		#chosen is the polytomy to be resolved
		labelClades<-lapply(prop.part(tree1),function(x) sort(tree1$tip.label[x]))
		chosen<-Ntip(tree1)+which(sapply(labelClades,function(x) identical(x,polyParts[[i]])))
		if(length(chosen)==0){stop("Cannot find a polytomy on the input tree??")}
		if(length(chosen)>1){stop("Selected more than one polytomy based on descendants??")}
		#
		ancMat<-ancPropStateMat(trait, tree1, orderedChar=orderedChar, type=type, cost=cost)
		#convert state/node matrix table from state weights to T/F, where T = max weight per node for a state
		maxWt<-t(apply(ancMat,1,function(x) x==max(x)))
		if(any(apply(maxWt,1,sum)==0)){stop("Error in ancestral state reconstruction, some have no state possible?")}
		#now get the chosen's children
		descChosen<-tree1$edge[tree1$edge[,1]==chosen,2]
		#now decide what state each desc and ancestor has
		if(is.null(stateBias)){				
			nodeChar<-apply(maxWt[sapply(rownames(maxWt),function(y) any(y==descChosen)),],1,function(x) which(x))
			ancChar<-which(maxWt[rownames(maxWt)==chosen,])				
		}else{
			if(orderedChar & (stateBias=="derived" | stateBias=="primitive")){
				#given max wt uncertainty for states, choose highest or lowest state if an ordered character
					#the lower would weight against reversals
				if(stateBias=="derived"){				
					nodeChar<-apply(maxWt[sapply(rownames(maxWt),function(y) any(y==descChosen)),],1,function(x) max(which(x)))
					ancChar<-max(which(maxWt[rownames(maxWt)==chosen,]))
					}
				if(stateBias=="primitive"){				
					nodeChar<-apply(maxWt[sapply(rownames(maxWt),function(y) any(y==descChosen)),],1,function(x) min(which(x)))
					ancChar<-min(which(maxWt[rownames(maxWt)==chosen,]))
					}
			}else{
				stop("Inconsistent orderedChar and stateBias arguments")
				}
			}
		#combine above into a single list
		foundStates<-c(nodeChar,list(ancChar))
		#find unique character groupings
			#if uncertainty, merge states into single state
		groupings<-unique(foundStates)
		#need to remove groupings that are short subsets of longer groupings
			#order by length first
		groupings<-groupings[order(-sapply(groupings,length))]
		k<-1
		#use while loop to go through groupings
		while(k<=length(groupings)){
			drop<-sapply(groupings,function(x) all(sapply(x,function(y) any(y==groupings[[k]]))))
			drop[k]<-FALSE
			#edit foundStates
			dropFS<-sapply(foundStates,function(x) any(sapply(groupings[drop],identical,x)))
			for(m in which(dropFS)){foundStates[[m]]<-groupings[[k]]}
			groupings<-groupings[!drop]
			k<-k+1
			}
		if(length(groupings)>1){
			#if only one group, just move on
			#
			#order groupings by average state
			groupings<-groupings[order(sapply(groupings,mean))]
			#assign groups, now in order
			groupAssign<-sapply(foundStates,function(x) which(sapply(groupings,identical,x)))
			#replace groupings with ranked grouping by taking unique of groupAssign
			groupings<-sort(unique(groupAssign))
			#
			#now we have a nodeChar vector with names = desc node IDs, values = ranked groupings
			nodeChar<-groupAssign[-length(groupAssign)]
			#and the ranked grouping for the ancestor
			ancChar<-groupAssign[length(groupAssign)]
			#
			#what if anc grouping isn't in nodeChar, add an artificial one
			if(!any(nodeChar==ancChar)){
				nodeChar<-c(nodeChar,'NaN'=ancChar)
				}
			#
			#now is time to build a tree
			#build a new edgeMat
			edgeMat<-matrix(,1,2)
			#for adding each level
			bottom<-chosen	#the ancestral node
			newsie<-max(tree1$edge)+1	#the new descendant node
			#
			# split algorithm into ordered and not ordered now
			if(orderedChar){
				#build a ladder tree with nested paraphyletic grades
				#first build anc to upwards ladder
				chainUp<-groupings[groupings>=ancChar]			
				#then down from anc ladder; reverse so going away from anc state
				chainDown<-rev(groupings[groupings<ancChar])
				#okay, chainUp
				for(j in 1:length(chainUp)){
					#which char is this 'level' of nodes going to use
					leveler<-chainUp[j]
					#let's make a polytomy for each unique level of nodeChar
					edgar<-cbind(bottom,as.numeric(names(nodeChar[nodeChar==leveler])))
					#add a new edge for chainDown (if such exists)
					if(length(chainDown)>0 & j==1){
						edgar<-rbind(edgar,c(bottom,newsie))
						savedBottom<-newsie	#this will be the new bottom
						if(j!=length(chainUp)){newsie<-newsie+1}
						}
					if(j!=length(chainUp)){
						#Need to add the edge for the polytomy of the next level: bottom,newsie
						edgar<-rbind(edgar,c(bottom,newsie))
						#need to update bottom and newsie
						bottom<-newsie   #the old newsie becomes the new ancestor (bottom)
						newsie<-newsie+1 #need to define a new descendant from that ancestor for attaching next level
						}
					edgeMat<-rbind(edgeMat,edgar)
					}
				if(length(chainDown)>0){ 
					bottom<-savedBottom
					newsie<-newsie+1
					for(j in 1:length(chainDown)){
						#which char is this 'level' of nodes going to use
						leveler<-chainDown[j]
						#let's make a polytomy for each unique level of nodeChar
						edgar<-cbind(bottom,as.numeric(names(nodeChar[nodeChar==leveler])))
						if(j!=length(chainDown)){
							#Need to add the edge for the polytomy of the next level: bottom,newsie
							edgar<-rbind(edgar,c(bottom,newsie))
							#need to update bottom and newsie
							bottom<-newsie   #the old newsie becomes the new ancestor (bottom)
							newsie<-newsie+1 #need to define a new descendant from that ancestor for attaching next level
							}
						edgeMat<-rbind(edgeMat,edgar)
						}
					}
				edgeMat<-edgeMat[-1,]
				if(nrow(edgeMat)!=(length(nodeChar)+length(unique(nodeChar))-1)){
					stop(paste("edgeMat is not the right size for polyParts=",i))}
				#
			}else{	#if unordered
				#
				#make every grouping but the ancestor a monophyletic cluster
					#old lolz: # stop("haven't done this yet")
				edgar<-cbind(bottom,as.numeric(names(nodeChar[nodeChar==ancChar])))
				#add to edgeMat
				edgeMat<-rbind(edgeMat,edgar)
				grouping<-groupings[groupings!=ancChar]
				for(j in 1:length(grouping)){
					#Need to add the edge for the polytomy of the next group: bottom,newsie
						#add directly to edgemat
					edgeMat<-rbind(edgeMat,c(bottom,newsie))	#bottom never changes
					#need to update bottom and newsie
					bottomG<-newsie   #the old newsie becomes the new ancestor (bottom)
					#which char is this 'level' of nodes going to use
					leveler<-grouping[j]
					#let's make a polytomy for each unique level of nodeChar
					edgar<-cbind(bottomG,as.numeric(names(nodeChar[nodeChar==leveler])))
					newsie<-newsie+1 #need to define a new descendant from that ancestor for attaching next level
					edgeMat<-rbind(edgeMat,edgar)
					}
				edgeMat<-edgeMat[-1,]
				}
			#drop any artificial ancestral taxa added to nodeChar
			edgeMat<-edgeMat[!is.nan(edgeMat[,2]),]
			#now need to clean old edge matrix, remove all with chosen as edge[,1], combine with edgeMat
			tree2<-tree1
			tree2$edge<-rbind(tree2$edge[tree2$edge[,1]!=chosen,],edgeMat)
			tree2$Nnode<-tree2$Nnode+length(unique(nodeChar))-1
			#tree2<-collapse.singles(reorder(tree2))
			#if(!testEdgeMat(tree2)){stop("Produced edge matrix has inconsistencies")}
			#tree3<-read.tree(text=write.tree(tree2))
			tree3<-cleanNewPhylo(tree2)
			if(Ntip(tree1)!=Ntip(tree3)){
				if(Ntip(tree1)!=Ntip(tree2)){
					stop("Taxa loss/added somehow in fixing polytomies?")
				}else{
					stop("Taxa loss/added somehow in fixing tree structure?")
					}
				}
			tree1<-tree3							
			}
		}
	#randRes remaining unresolved areas: yes, no?
		#no, leave this alone, can be a second line of code, for god sakes
	treeFinal<-ladderize(tree1,right=FALSE)
	return(treeFinal)
	}