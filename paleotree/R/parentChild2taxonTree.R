#' Create a Taxonomy-Based Phylogeny ('Taxon Tree') from a Table of Parent-Child Taxon Relationships
#'
#' This function takes a two-column matrix of taxon names,
#' indicating a set of binary parent-taxon:child-taxon 
#' paired relationships with a common root, and returns
#' a 'taxonomy-tree' phylogeny object of class 'phylo'.

#' @details
#' All taxa listed must be traceble via their parent-child relationships to a single,
#' common ancestor which will act as the root node for output phylogeny. Additionally,
#' the root used will be the parent taxon to all tip taxa closest in terms of parent-child
#' relationships to the tip taxa: i.e., the most recent common ancestor. Ancestral taxa which
#' are singular internal nodes that trace to this root are removed, and a message
#' is printed.

#' @param parentChild A two-column matrix of type \code{character} where
#' each element is a taxon name. Each row represents a parent-child relationship
#' with first the parent (column 1) taxon name and then the child (column 2).

#' @param tipSet This argument controls which taxa are selected as tip taxa for the
#' output tree. The default \code{tipSet="nonParents"} selects all child taxa which
#' are not listed as parents in \code{parentChild}. Alternatively, \code{tipSet="all"}
#' will add a tip to every internal node with the parent-taxon name encapsulated in
#' parentheses.

# @param reorderTree A logical indicating whether a step of \code{reorder.phylo()} will be applied,
# if \code{cleanTree=TRUE}; has no effect if \code{cleanTree=FALSE}.
# Reordering may cause more problems than it is worth in older versions of \code{ape}.

#' @inheritParams makePBDBtaxonTree

#' @return
#' A phylogeny of class 'phylo', with tip taxa as controlled by argument \code{tipSet}.
#' The output tree is returned with no edge lengths.
#'
#' The names of higher taxa than the tips should be appended as the element $node.label for the internal nodes.

#' @seealso \code{\link{makePBDBtaxonTree}}, \code{\link{taxonTable2taxonTree}}

#' @author David W. Bapst

#' @examples
#' 
#' #let's create a small, really cheesy example
#' pokexample<-rbind(cbind("Squirtadae",c("Squirtle","Blastoise","Wartortle")),
#' 	c("Shelloidea","Lapras"),c("Shelloidea","Squirtadae"),
#' 	c("Pokezooa","Shelloidea"),c("Pokezooa","Parasect"),
#' 	c("Rodentapokemorpha","Linoone"),c("Rodentapokemorpha","Sandshrew"),
#' 	c("Rodentapokemorpha","Pikachu"),c("Hirsutamona","Ursaring"),
#' 	c("Hirsutamona","Rodentapokemorpha"),c("Pokezooa","Hirsutamona"))
#' 
#' #Default: tipSet='nonParents'
#' pokeTree<-parentChild2taxonTree(pokexample, tipSet="nonParents")
#' plot(pokeTree);nodelabels(pokeTree$node.label)
#'
#' #Get ALL taxa as tips with tipSet='all'
#' pokeTree<-parentChild2taxonTree(pokexample, tipSet="all")
#' plot(pokeTree);nodelabels(pokeTree$node.label)
#'
#' 
#' \dontrun{
#' 
#' # let's try a dataset where not all the taxon relationships lead to a common root
#' 
#' pokexample_bad<-rbind(cbind("Squirtadae",c("Squirtle","Blastoise","Wartortle")),
#' 	c("Shelloidea","Lapras"),c("Shelloidea","Squirtadae"),
#' 	c("Pokezooa","Shelloidea"),c("Pokezooa","Parasect"),
#' 	c("Rodentapokemorpha","Linoone"),c("Rodentapokemorpha","Sandshrew"),
#' 	c("Rodentapokemorpha","Pikachu"),c("Hirsutamona","Ursaring"),
#' 	c("Hirsutamona","Rodentapokemorpha"),c("Pokezooa","Hirsutamona"),
#' 	c("Umbrarcheota","Gengar"))
#' 
#' #this should return an error, as Gengar doesn't share common root
#' pokeTree<-parentChild2taxonTree(pokexample_bad)
#' 
#' 
#' # another example, where a taxon is listed as both parent and child
#' pokexample_bad2<-rbind(cbind("Squirtadae",c("Squirtle","Blastoise","Wartortle")),
#' 	c("Shelloidea",c("Lapras","Squirtadae","Shelloidea")),
#' 	c("Pokezooa","Shelloidea"),c("Pokezooa","Parasect"),
#' 	c("Rodentapokemorpha","Linoone"),c("Rodentapokemorpha","Sandshrew"),
#' 	c("Rodentapokemorpha","Pikachu"),c("Hirsutamona","Ursaring"),
#' 	c("Hirsutamona","Rodentapokemorpha"),c("Pokezooa","Hirsutamona"),
#' 	c("Umbrarcheota","Gengar"))
#'
#' #this should return an error, as Shelloidea is its own parent
#' pokeTree<-parentChild2taxonTree(pokexample_bad2)
#'
#' }
#' 
#' 
#' 
#' # note that we should even be able to do this with ancestor-descendent pairs from
#' 	 # simulated datasets from simFossilRecord, like so:
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' # need to reorder the columns so parents (ancestors) first, then children 
#' parentChild2taxonTree(taxa[,2:1])
#' # now note that it issues a warning that the input wasn't type character
#'    # and it will be coerced to be such
#'


#' @name parentChild2taxonTree
#' @rdname parentChild2taxonTree
#' @export
parentChild2taxonTree<-function(parentChild,tipSet="nonParents",cleanTree=TRUE){
	#,reorderTree=TRUE
	#
	#takes a two column matrix of character class taxon names
		#each row is a relationship: parent, then child
	#CHECKS
	if(length(tipSet)!=1 | !is.character(tipSet)){stop("tipSet must be a single character element")}
	if(!is.character(parentChild)){
		message("parentChild isn't of class character, attempting to coerce")
		parentChild<-apply(parentChild,2,as.character)}
	if(length(dim(parentChild))!=2){
		stop("parentChild must be a matrix of class character with two columns and multiple rows")
	}else{
		if(!(dim(parentChild)[2]==2 & dim(parentChild)[1]>1)){
			stop("parentChild must be a matrix of class character with two columns and multiple rows")
			}
		}
	#
	if(!testParentChild(parentChild=parentChild)){stop("parentChild relationships are inconsistent")}
	#
	#remove singular root edges
	#trace tips to ultimate ancestor (should be same for all, as this has already been checked)
	continue<-TRUE
	while(continue){
		unqIDs<-unique(c(parentChild[,1],parentChild[,2]))
		ultimateAnc<-sapply(unqIDs,getUltimateAnc,parentChild=parentChild)
		if(length(unique(ultimateAnc))==1){
			ultAnc1<-ultimateAnc[1]
		}else{
			stop("parentChild constructed improperly")}
		descEdge<-which(parentChild[,1]==ultAnc1)
		if(length(descEdge)==1){
			message(paste("Removing singular node leading to root:",ultAnc1))
			#remove from parentChild
			parentChild<-parentChild[-descEdge,,drop=FALSE]
			if(nrow(parentChild)<1){stop("No branching nodes found?!")}
		}else{
			continue=FALSE
			}
		}
	#
	#first, get nodeNames, with root name first
	nodeNames<-unique(parentChild[,1])
	whichRoot<-which(sapply(nodeNames,function(x) 
		!any(sapply(parentChild[,2],identical,unname(x)))))
	#check that there isn't more than one root
	if(length(whichRoot)>1){
		stop(paste("Not all taxable are traceable to a single common root \n",
			length(whichRoot),"possible roots found:",
			paste0(nodeNames[whichRoot],collapse=", ")))}
	#now resort nodeNames
	nodeNames<-c(nodeNames[whichRoot],nodeNames[-whichRoot])
	if(tipSet!="nonParents"){
		if(tipSet=="all"){
			parentChild<-rbind(parentChild,cbind(nodeNames,paste0("(",nodeNames,")")))
		}else{
			stop("tipSet must be one of either 'nonParents' or 'all'")}
			}
	#identify tip taxa, this will be all taxa who are not-parents
	notParents<-sapply(parentChild[,2],function(x) 
		!any(sapply(parentChild[,1],identical,unname(x))))
	tipNames<-parentChild[notParents,2]
	#now convert parentChild matrix to edge matrix
	edgeMat<-matrix(,nrow(parentChild),ncol(parentChild))
	taxonNames<-c(tipNames,nodeNames)
	#test that none have been lost
	if(length(taxonNames)!=length(unique(c(parentChild[,1],parentChild[,2])))){
		stop("Number of tip and node names doesn't sum to total number of unique names in parentChild")}
	#convert internal nodes to Ntip+nodeNames ID
	edgeMat[,1]<-sapply(parentChild[,1],function(x)
		which(sapply(taxonNames,identical,x)))
	edgeMat[,2]<-sapply(parentChild[,2],function(x) 
		which(sapply(taxonNames,identical,x)))
	#reorder edge
	edge<-edgeMat[order(edgeMat[,1],edgeMat[,2]),]
	#check edge
	if(!testParentChild(parentChild=edge)){
		stop("created edge relationships are inconsistent")}	
	#make the tree
	tree<-list(edge=edge,tip.label=tipNames,edge.length=NULL, #edge.length=rep(1,nrow(edge))
		Nnode=length(nodeNames),node.label=nodeNames)
	class(tree)<-"phylo"
	if(cleanTree){ #make it a good tree
		#reordering seems to cause errors?? 06-11-15
		tree<-cleanNewPhylo(tree)
		}
	if(Ntip(tree)!=length(tipNames)){stop("Taxa number changed while cleaning tree")}
	#plot(tree);nodelabels(tree$node.label)
	return(tree)
	}

getUltimateAnc<-function(taxa,parentChild){
	count<-0
	while(any(sapply(parentChild[,2],identical,unname(taxa)))){
		count<-count+1
		taxa<-parentChild[match(taxa,parentChild[,2]),1]
		if(length(taxa)>1){
			stop("Some parents are listed as children twice in parentChild")}
		if(count>(length(parentChild)*2)){
			stop("Breaking while() loop: cannot find ultimate ancestor")
			}
		}
	return(taxa)
	}

charEdge2numeric<-function(parentChild){
	if(!is.character(parentChild)){
		stop("parentChild has to be character for charEdge2numeric")}
	#get unique IDs
	unqIDs<-c(NA,sort(unique(c(parentChild[,1],parentChild[,2]))))
	parentChild2<-matrix(,nrow(parentChild),2)
	#convert internal nodes to Ntip+nodeNames ID
	parentChild2[,1]<-sapply(parentChild[,1],function(x)
		which(sapply(unqIDs,identical,x))-1)
	parentChild2[,2]<-sapply(parentChild[,2],function(x) 
		which(sapply(unqIDs,identical,x))-1)	
	if(!is.numeric(parentChild2)){
		stop("parentChild not coercing correctly to numeric for charEdge2numeric")}
	return(parentChild2)
	}

testParentChild<-function(parentChild){
	#check that its a matrix with two columns
	if(!is.matrix(parentChild)){
		stop("edge/parentChild matrix must be of type matrix")}
	if(ncol(parentChild)!=2){
		stop("edge/parentChild matrix must have two columns")}
	#convert to numeric
	if(!is.numeric(parentChild)){
		if(is.character(parentChild)){
			parentChild<-charEdge2numeric(parentChild)
		}else{
			stop("input must be of type numeric or type character")
			}			
		}
	#
	#replace NA values
	parentChild[is.na(parentChild)]<-min(parentChild,na.rm=TRUE)-1
	#test monophyly of parentChild
		#test that all but one node has an ancestor
	parentMatch<-match(unique(parentChild[,1]),parentChild[,2])
	if(sum(is.na(parentMatch))>1){
		stop(paste("More than one apparent root; \n",
			"more than one parent without their own parent listed"))}
	#
	#check that any ancestor is listed as its own descendant
	parentANDchild<-apply(parentChild,1,function(x){
		x<-unname(x)
		identical(x[1],x[2])
		})
	if(any(parentANDchild)){
		stop("Some pairs with same ID listed for both parent and child (?!)")
		}
	#
	#trace all tips to a single ancestor
	unqIDs<-unique(c(parentChild[,1],parentChild[,2]))
	ultimateAnc<-sapply(unqIDs,getUltimateAnc,parentChild=parentChild)
	if(length(unique(ultimateAnc))!=1){
		stop("IDs trace back to more than one unique common ancestor")}
	#
	#test for nodes listed as descendant twice
	if(any(table(parentChild[,2])>1)){
		stop(paste("Some IDs are listed as a descendant twice in the edge/parentChild matrix",
		paste0(parentChild[duplicated(parentChild[,2]),2],collapse=", ")))}
	#test that all but one node has an ancestor
	parentMatch<-match(unique(parentChild[,1]),parentChild[,2])
	if(sum(is.na(parentMatch))>1){
		stop(paste("More than one apparent root; \n",
			"more than ancestor without an ancestor listed"))}	
	return(TRUE)
	}
