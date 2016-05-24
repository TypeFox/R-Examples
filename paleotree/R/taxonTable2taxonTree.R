#' Create a Taxonomy-Based Phylogeny ('Taxon Tree') from a Hierarchical Table of Taxonomy Memberships
#'
#' This function takes a matrix of taxon names,
#' indicating a set of hierarchical taxonomic relationships
#' conveyed as nested placements for a set of tip-taxa (listed in
#' the last column of the matrix) and returns
#' a 'taxonomy-tree' phylogeny object of class 'phylo'.
 
#' @details
#' This function can deal with empty entries in cells of \code{taxonTable} by assuming these
#' are lower-level taxa which are 'floating' freely somewhere in
#' taxa several levels higher.

#' @inheritParams makePBDBtaxonTree

#' @param taxonTable A matrix of type character and multiple rows and columns, containing the
#' tip taxa in the last column, one per row, with progressively larger taxa listed in prior
#' columns (reading left-to-right). Invariant columns (i.e. taxa that all
#' tip taxa are in) are allowed, but all but the most 'shallow' of such invariant taxa are
#' dropped prior to transformation to a taxon-tree phylogeny object.

#' @return
#' A phylogeny of class 'phylo', where each tip is a taxon listed in the last column of the
#' input taxonTable. Edges are scaled so that
#' the distance from one taxon rank to another 1, then merged to remove singleton nodes. As not all
#' taxa have parents at the immediate taxon level above, this leads to some odd cases. For example,
#' two genera emanating from a node representing a class but with a very short (length=1) branch
#' and a long branch (length=3) means one genus is simply placed in the class, with no family or order listed
#' while the one on the long branch is within an order and family that is otherwise monogeneric.
#'
#' The names of higher taxa than the tips should be appended as the element $node.label for the internal nodes.

#' @seealso \code{\link{makePBDBtaxonTree}}, \code{\link{parentChild2taxonTree}}

#' @author David W. Bapst

#' @examples
#'
#' #let's create a small, really cheesy example
#' pokeTable<-rbind(cbind("Pokezooa","Shelloidea","Squirtadae",
#' 		c("Squirtle","Blastoise","Wartortle")),
#' 	c("Pokezooa","Shelloidea","","Lapras"),
#' 	c("Pokezooa","","","Parasect"),
#' 	cbind("Pokezooa","Hirsutamona","Rodentapokemorpha",
#' 		c("Linoone","Sandshrew","Pikachu")),
#' 	c("Pokezooa","Hirsutamona",NA,"Ursaring"))
#' 
#' pokeTree<-taxonTable2taxonTree(pokeTable)
#' plot(pokeTree);nodelabels(pokeTree$node.label)
#'

#' @name taxonTable2taxonTree
#' @rdname taxonTable2taxonTree
#' @export
taxonTable2taxonTree<-function(taxonTable,cleanTree=TRUE){
	# taxonTable<-taxonData
	#CHECKS
	if(length(dim(taxonTable))!=2 | !is.character(taxonTable)){
		stop("taxonTable must be a matrix of class character with multiple columns and rows")
		}
	if(all(dim(taxonTable)<2)){stop("taxonTable must have multiple rows and columns")}
	#
	#turn NAs to blanks
	taxonTable[is.na(taxonTable)]<-""
	#
	#remove constant columns, as they are meaningless
	constantCol<-apply(taxonTable,2,function(x) all(x==x[1]))
	# except for the last (which will be the root)
	constantCol[rev(which(constantCol))[1]]<-FALSE
	taxonTable<-taxonTable[,!constantCol,drop=FALSE]
	#check that enough columns remain
	if(ncol(taxonTable)<2){stop("Must have at least one varying taxonomy columns in taxonTable")}
	#
	#turn blanks to NAs
	taxonTable[taxonTable==""]<-NA
	#are there any missing names
	if(any(is.na(taxonTable[,ncol(taxonTable)]))){
		stop("Missing tip taxon names in taxonTable?")}
	#tip taxon labels
	labels<-taxonTable[,ncol(taxonTable)]
	edge<-matrix(NA,,2)
	#need to define columns BACKWARDS..
	levels<-rev(1:ncol(taxonTable))[-1]
	for(level in levels){
		newNodes<-unique(taxonTable[,level])[!is.na(unique(taxonTable[,level]))]
		for(node in newNodes){
			#if(node=="Nodosauridae"){stop("hey")}
			labels<-c(labels,node)
			nodeID<-length(labels)
			#find the descendant nodes
			whichDesc<-which(taxonTable[,level]==node)
			descTaxRow<-lapply(whichDesc,function(x) (taxonTable[x,-(1:level)]))
			descTax<-sapply(descTaxRow,function(x) x[!is.na(x)][1])
			descTax<-unique(descTax)
			descID<-sapply(descTax,function(x) which(x==labels))
			names(descID)<-NULL
			newEdge<-cbind(nodeID,descID)
			edge<-rbind(edge,newEdge)
			}
		}
	edge<-edge[-1,]
	edge.length<-rep(1,nrow(edge))
	Nnode<-length(unique(edge[,1]))
	tip.label<-taxonTable[,ncol(taxonTable)]
	node.label<-labels[-(1:length(tip.label))]
	Ntip<-length(tip.label)
	#need to flip node numbers
	nodes<-sort(unique(edge[,1]))
	nodes<-cbind(c(1:Ntip,nodes),c(1:Ntip,rev(nodes)))
	edge[,1]<-sapply(edge[,1],function(x) nodes[x==nodes[,1],2])
	edge[,2]<-sapply(edge[,2],function(x) nodes[x==nodes[,1],2])
	#check root
	nodes<-sort(unique(edge[,1]))
	root<-nodes[sapply(nodes,function(x) all(x!=edge[,2]))]
	if(root!=(Ntip+1)){stop("Root isn't renumbering correctly")}
	#reorder edge
	edge<-edge[order(edge[,1],edge[,2]),]
	#make the tree
	tree<-list(edge=edge,tip.label=tip.label,edge.length=edge.length,
		Nnode=Nnode,node.label=rev(node.label))	
	class(tree)<-"phylo"
	if(cleanTree){ #make it a good tree
		tree<-cleanNewPhylo(tree)
		}
	if(Ntip(tree)!=nrow(taxonTable)){stop("Taxa number changed while cleaning tree")}
	#plot(tree)
	return(tree)
	}