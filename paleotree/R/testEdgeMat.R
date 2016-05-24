#' Test the Edge Matrix of a 'phylo' Phylogeny Object for Inconsistencies
#'
#' \code{testEdgeMat} is a small simple function which tests the $edge matrix of 'phylo' objects for
#' inconsistencies that can cause downstream analytical problems.
#' The associated function, \code{cleanNewPhylo} puts an input
#' phylo object, presumably freshly created or reconstituted by some function, through a series
#' of post-processing, This includes having singles collapsed,
#' nodes reordered and being written out as a Newick string and read back in,
#' to ensure functionality with ape functions
#' and ape-derived functions. 

#' @aliases cleanNewPhylo cleanTree

#' @details
#' Useful when doing complex manipulations and reconstitutions of 'phylo' objects (or their
#' de novo construction), and thus is used by a number of paleotree functions.

#' @param tree A phylogeny object of type phylo

# @param reorderTree A logical indicating whether a step of \code{reorder.phylo()} will be applied.
# Reordering may cause more problems than it is worth.

#' @return
#' For \code{testEdgeMat}, if all the checks in the function pass correctly, the logical TRUE is returned.
#'
#' For \code{cleanNewPhylo}, an object of class 'phylo' is returned.

#' @author 
#' David W. Bapst, with a large number of tests incorporated from Emmanuel Paradis's checkValidPhylo function,
#' provided at his github repo here, which was released GPL v>2:
#' 
#' https://github.com/emmanuelparadis/checkValidPhylo

#' @examples
#'
#' set.seed(444)
#' tree<-rtree(10)
#' # should return TRUE
#' testEdgeMat(tree)
#'
#' # should also work on star trees
#' testEdgeMat(stree(10))
#'
#' # should also work on trees with two taxa
#' testEdgeMat(rtree(2))
#'
#' # should also work on trees with one taxon
#' testEdgeMat(stree(1))
#'
#' #running cleanNewPhylo on this tree should have little effect
#' 		#beyond ladderizing it...
#' tree1<-cleanNewPhylo(tree)
#'
#' #compare outputs
#' layout(1:2)
#' plot(tree)
#' plot(tree1)
#' layout(1)


#' @name testEdgeMat
#' @rdname testEdgeMat
#' @export testEdgeMat
testEdgeMat<-function(tree){
	#CHECKS
	if(!is(tree,"phylo")){stop("tree is not of type 'phylo'")}
	if(any(is.na(match(c("edge","tip.label","Nnode"),names(tree))))){
		stop("Missing key required elements of a 'phylo' object")}
	#MASSIVE SET OF TESTS FROM PARADIS'S checkValidPhylo, added 06-15-15
	#check that there are no NAs in $edge
	if(any(is.na(tree$edge))){
		stop("NA values in $edge table")}
	#check that there is a $tip.label and its a vector of type character with length>0
	if(!is.vector(tree$tip.label)){
		stop("$tip.label must be a vector")}
	if(!is.character(tree$tip.label)){
		stop("$tip.label must be of type character")}
	if(length(tree$tip.label)<1){
		stop("$tip.label must be of length greater than 0")}
	#check than Nnode exists, is a vector of length 1, of type number, stored as an integer
	if(!is.vector(tree$Nnode)){
		stop("$Nnode must be a vector")}
	if(!is.numeric(tree$Nnode)){
		stop("$Nnode must be of type numeric")}
	if(length(tree$Nnode)!=1){
		stop("$Nnode must be of length 1")}
	if(tree$Nnode<1){
		stop("$Nnode must be at least 1")}
	#test edge matrix has two columns, is numeric
	if(!is.matrix(tree$edge)){
		stop("$edge must be of type matrix")}
	if(ncol(tree$edge)!=2){
		stop("$edge must have two columns")}
	if(!is.numeric(tree$edge)){
		stop("$edge must be of type numeric")}
	if(any(is.na(tree$edge))){
		stop("NAs found in $edge matrix")}
	#test that edge and Nnode is stored as integers
	if(storage.mode(tree$Nnode)!="integer"){
		stop("$Nnode is not stored as an integer")}
	if(storage.mode(tree$edge)!="integer"){
		stop("$Nnode is not stored as an integer")}
	#check values in $edge for bad values
	if(any(is.na(tree$edge))){
		stop("NAs found in $edge matrix")}	
	if(any(tree$edge<1)){
		stop("All elements of $edge must be integers of 1 or greater")}
	#check that root node and Nnode defined correctly
	rootID<-getRootID(tree)
	if(rootID!=(Ntip(tree)+1)){
		stop(paste0("The root node is numbered ",rootID," in $edge, should be Ntip(tree)+1 (",Ntip(tree)+1,")"))}
	#expected number of tips and nodes is Nnode+length($tip.label)
	expNodeNumber<-tree$Nnode+length(tree$tip.label)
	if(any(tree$edge>expNodeNumber)){
		stop(paste0("Some elements of edge are numbered greater than ",
			expNodeNumber,"\n calculated from Nnode+length(tree$tip.label)","\n Check these: ",
			paste(unique(tree$edge[tree$edge>expNodeNumber]),collapse=" ")))}
	#test if Nnode agrees with number of nodes in $edge
	if(Nnode(tree)!=(max(tree$edge[,1])-Ntip(tree))){
		stop("Nnode is lower than number implied by edge[,1]?")}
	#now switch to tabulate based checking of node IDs in $edge (stolen from Paradis)
	tabEdge<-tabulate(tree$edge)
	#are all expected node IDs found in tabEdge?
	if(length(tabEdge)<expNodeNumber){
		stop("Fewer node IDs found in $edge than expected from Nnode+length(tree$tip.label)")}
	#do all tips only appear once?
	if(any(tabEdge[1:length(tree$tip.label)]>1)){
		stop("Some tip IDs appear more than once in $edge")}
	if(any(tabEdge[1:length(tree$tip.label)]<1)){
		stop("Some tip IDs appear less than once in $edge")}
	#
	#All internal nodes should appear at least once (even if singleton nodes)
	if(tree$Nnode>1){
		if(any(tabEdge[length(tree$tip.label) + 2:tree$Nnode]<2)){
			stop("Some internal node IDs appear less than twice in $edge: ")}
		if(any(tabEdge[length(tree$tip.label) + 2:tree$Nnode]<2)){
			stop(paste0("Some internal node IDs appear less than twice in $edge: ",
				which(tabEdge[length(tree$tip.label) + 1:tree$Nnode]<2),
				collapse=" "))}
		}
	#check that tips do not appear in tree$edge[,1]
	if(any(tree$edge[,1]<(length(tree$tip.label) + 1))){
		stop(paste0("Apparent tip IDs appear in column 1 of $edge: \n",
		tree$edge[tree$edge[,1]<(length(tree$tip.label) + 1),1],collapse=" "))}
	#test edge matrix
	if(!testParentChild(parentChild=tree$edge)){stop("Edge matrix has inconsistencies")}
	#more tests of edge matrix
	if(Ntip(tree)>2){
		if(Nnode(tree)!=(max(tree$edge)-Ntip(tree))){
			stop("Number of nodes is incorrect based on edge numbering?")}
		if(Nnode(tree)>1){
			#is every internal node listed as a descendant and ancestor, in edge[,2] and edge[,1]?
			if(!all(sapply((1:Nnode(tree))+Ntip(tree),function(x) any(x==tree$edge[,1])))){
				stop("Not all internal nodes (including root) listed in edge[,1]?")}
			if(sum(!sapply((1:Nnode(tree))+Ntip(tree),function(x) any(x==tree$edge[,2])))>1){
				stop("Not all internal nodes (except root) listed in edge[,2]?")}
			}
		}
	#if(identical(sort(unique(tree$edge[,2])),c(1L,2L))){stop("Number of nodes is incorrect based on edge[,2]?")}
	return(TRUE)
	}
	
#hidden function
getRootID<-function(tree){
	uniqueNode<-unique(tree$edge[,1])
	whichRoot<-sapply(uniqueNode,function(x)
		(sum(x==tree$edge[,2])==0))
	if(sum(whichRoot)>1){
		stop("More than one apparent root in $edge matrix")}
	rootID<-uniqueNode[whichRoot]
	return(rootID)
	}

#' @rdname testEdgeMat
#' @export 
cleanNewPhylo<-function(tree){ 		#,reorderTree=TRUE
		#
		renumberRootID<-function(tree){
			rootID<-getRootID(tree)
			expRootID<-length(tree$tip.label)+1
			if(rootID!=expRootID){
				tree$edge[tree$edge==rootID]<-0
				tree$edge[tree$edge==expRootID]<-rootID
				tree$edge[tree$edge==0]<-expRootID
				storage.mode(tree$edge)<-"integer"
				}
			return(tree)
			}
		#
		#CHECKS
		if(!inherits(tree,"phylo")){
			stop("tree must be of class 'phylo'")
			}
		if(any(is.na(match(c("edge","tip.label","Nnode"),names(tree))))){
			stop("Missing key required elements of a 'phylo' object")}
		oldNtip<-length(tree$tip.label)
		#make it a good tree
		#coerce edge and Nnode to storage mode for integers
		if(storage.mode(tree$edge)!="integer"){
			storage.mode(tree$edge)<-"integer"
			}
		if(storage.mode(tree$Nnode)!="integer"){
			storage.mode(tree$Nnode)<-"integer"
			}
		#check it
		if(!testEdgeMat(tree)){stop("Edge matrix has inconsistencies")}
		#collapse singles
			#count number of single nodes
		Nsingle<-sum(sapply(unique(tree$edge[,1]),function(x) sum(x==tree$edge[,1])==1))
		if(Nsingle>0){
			treePrev<-tree
			while(Nsingle>0){
				tree1<-collapse.singles(treePrev)
				#renumber root if numbered incorrectly
				tree1<-renumberRootID(tree1)
				if(!testEdgeMat(tree1)){stop("Edge matrix has inconsistencies")}
				if((Nnode(treePrev)-Nnode(tree1))>Nsingle){
					stop("collapse.singles dropped too many nodes")}
				Nsingle<-sum(sapply(unique(tree1$edge[,1]),function(x) sum(x==tree1$edge[,1])==1))
				treePrev<-tree1
				#print("for counting how many times singles need to be dropped")
				}
		}else{
			tree1<-tree
			}
		#
		#reorder	#if(reorderTree){
		tree1<-reorder.phylo(tree1,"cladewise") 	#REORDER IT
		if(!testEdgeMat(tree1)){stop("Edge matrix has inconsistencies")}
		tree1<-read.tree(text=write.tree(tree1))
		if(!testEdgeMat(tree1)){stop("Edge matrix has inconsistencies")}
		tree1<-ladderize(tree1)
		#test tip numbers
		if(oldNtip!=Ntip(tree1)){stop("Final tip taxon number different from original number of tip taxon names")}
		return(tree1)
		}