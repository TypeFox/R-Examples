#' Extrapolating Lower-Level Taxon Phylogenies from Higher-Level Taxon Trees
#' 
#' This function takes a tree composed of higher-level taxa and a vector of
#' lower-level taxa belonging to the set of higher-level taxa included in the
#' input tree and produces a tree composed of the lower-level taxa, by treating
#' the higher-level taxa as unresolved monophyletic polytomies. A user can also
#' mark higher taxa as paraphyletic such that these are secondarily collapsed
#' and do not form monophyletic clades in the output tree.
#' 
#' @details The output tree will probably be a rough unresolved view of the
#' relationships among the taxa, due to the treatment of higher-level taxa as
#' polytomies. This is similar to the methods used in Webb and Donoghue (2005)
#' and Friedman (2009). Any analyses should be done by resolving this tree with
#' \code{\link{multi2di}} in the ape package or via the various time-scaling
#' functions found in this package (paleotree).
#' 
#' The taxaData vector should have one element per lower-level taxon that is to
#' be added to the tree. The name of each element in the vector should be the
#' names of the lower-level taxa, which will be used as new tip labels of the
#' output lower-taxon tree. There should be no empty elements! expandTaxonTree
#' won't know what to do with taxa that don't go anywhere.
#' 
#' By default, all higher-level taxa are treated as monophyletic clades if not
#' otherwise specified. The collapse vector can (and probably should) be used
#' if there is doubt about the monophyly of any higher-level taxa included in
#' the input taxon-tree, so that such a group would be treated as a paraphyletic
#' group in the output tree.
#' 
#' Also by default, the output tree will lack branch lengths and thus will not be
#' time-scaled. If keepBrLen is true, then the tree's edge lengths are kept and
#' new taxa are added as zero length branches attaching to a node that
#' represents the previous higher-taxon. This tree is probably not useful for
#' most applications, and may even strongly bias some analyses. USE WITH
#' CAUTION! The 'collapse' vector will cause such edges to be replaced by
#' zero-length branches rather than fully collapsing them, which could have odd
#' effects. If 'collapse' is not null and keepBrLen is true, a warning is
#' issued that the output probably won't make much sense at all.
#' 
#' @param taxonTree A phylo object where tips represent higher taxa
#' @param taxaData Character vector of higher taxa, with elements names equal
#' to the lower taxa. See below.
#' @param collapse Character vector of non-monophyletic higher taxa to be
#' collapsed
#' @param keepBrLen Logical, decides if branch lengths should be kept or
#' discarded. FALSE by default. See details below.
#' @param plot If true, plots a comparison between input and output trees
#' @return Outputs the modified tree as an object of class phylo, with the
#' higher-level taxa expanded into polytomies and the lower-level taxa as the
#' tip labels.
#' @author David W. Bapst
#' @seealso \code{\link{multi2di}}, \code{\link{bind.tree}}
#' @references Friedman, M. 2009 Ecomorphological selectivity among marine
#' teleost fishes during the end-Cretaceous extinction. \emph{Proceedings of
#' the National Academy of Sciences} \bold{106}(13):5218--5223.
#' 
#' Webb, C. O., and M. J. Donoghue. 2005 Phylomatic: tree assembly for applied
#' phylogenetics. \emph{Molecular Ecology Notes} \bold{5}(1):181--183.
#' @examples
#' 
#' set.seed(444)
#' #lets make our hypothetical simulated tree of higher taxa
#' taxtr <- rtree(10)
#' taxd <- sample(taxtr$tip.label,30,replace=TRUE)	#taxa to place within higher taxa
#' names(taxd) <- paste(taxd,"_x",1:30,sep="")
#' coll <- sample(taxtr$tip.label,3)		#what to collapse?
#' expandTaxonTree(taxonTree=taxtr,taxaData=taxd,collapse=coll,plot=TRUE)
#' 
#' @export expandTaxonTree
expandTaxonTree<-function(taxonTree,taxaData,collapse=NULL,keepBrLen=FALSE,plot=FALSE){
	#this function takes a higher-level taxon tree and
		#expands it to a lower level species-level tree
		#using a species list
	# term 'taxa' here represents the groups to be replaced on the taxonTree
	#taxonTree = tree with taxon IDs as tips
	#taxaData = character vector of higher taxon ids for each new tip, tip labels as vector names
	#collapse = if present, vector of taxa names to be collapsed
	#should be possible to take a tree of mixed species/genera
		#and just replace the genera
	#taxonTree<-rtree(10);taxonTree$tip.label<-as.character(1:10);collapse<-sample(taxonTree$tip.label,5)
	#taxaData<-as.character(sample(1:10,100,replace=TRUE));names(taxaData)<-paste("t",1:100,sep="")
	#require(ape)
	if(!inherits(taxonTree, "phylo")){
		stop("taxonTree is not of class phylo")
		}
	if(any(is.na(taxaData))){stop("some values of taxonData missing!!")}
	if(!is.null(collapse) & keepBrLen){
		message("Warning: collapsed branch lengths turned to zero when keepBrLen is TRUE!")}
	tree<-taxonTree
	if(!keepBrLen){tree$edge.length<-rep(1,Nedge(tree))}		#get rid of all branch lengths
	#first, expand all higher taxa to lower taxon polytomies
	for(i in unique(taxaData)){				#loop through all 	
		tip<-which(tree$tip.label==i)
		if(length(collapse)>0){if(any(collapse==i)){
			tree$edge.length[which(tree$edge[,2]==tip)]<-0
			}}
		cotaxa<-names(taxaData)[taxaData==i]	#which species do I want? These...
		repTree<-stree(length(cotaxa))		#replacement polytomy
		if(keepBrLen){
			repTree$edge.length<-rep(0,length(cotaxa))
		}else{
			repTree$edge.length<-rep(1,length(cotaxa))
			}
		repTree$tip.label<-cotaxa			#replace names,edge.lengths
		tree<-bind.tree(tree,repTree,tip)	#replace the right tip	
		}
		#now collapse non-monophyletic groupings
	if(!keepBrLen){
		tree1<-di2multi(tree)
		tree1$edge.length<-NULL
		#tree1<-collapse.singles(tree1)
		#if(!testEdgeMat(tree1)){stop("Edge matrix has inconsistencies")}
		#tree1<-read.tree(text=write.tree(tree1))
		tree1<-cleanNewPhylo(tree1)
	}else{
		tree1<-tree
		}
	if(plot==TRUE){layout(1:2);plot(taxonTree);plot(tree1);layout(1)}
	return(tree1)
	}
