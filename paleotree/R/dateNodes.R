#' Absolute Dates for Nodes of a Time-Scaled Phylogeny
#' 
#' This function returns the ages of nodes (both internal and terminal tips) for a given phylogeny
#' of class 'phylo'. Specialized for use with time-scaled trees from paleotree, see Details.

#' @details
#' This function is specialized for phylo objects time-scaled or simulated with functions from
#' paleotree, and thus have a $root.time element. This function will still work without such,
#' but users should see the details for the \code{rootAge} argument.

#' @param tree A phylogeny object of class 'phylo'. Must have edge.lengths!

#' @param rootAge The root age of the tree, assumed by default to be equal to the element
#' tree$root.time, which is a standard element for trees time-scaled by the paleotree
#' package. If not given by the user and if the $root.time element does not exist, then
#' the maximum depth of the tree will be taken as the root age, which implicitly assumes
#' the latest most terminal tip is an extant taxon at the modern day (time = 0). If rootAge
#' is so defined that some nodes may occur later than time = 0, this function may return
#' negative dates.

#' @param labelDates If FALSE (the default), the dates returned are labeled with the
#' tip/node numbers as in \code{tree$edge}. If TRUE, they are labeled with the tip labels
#' of every descendant tip, which for terminal tips means a single taxon label, and for
#' internal tips a label that might be very long, composed of multiple tip labels pasted
#' together. Thus, by default, this argument is FALSE.

#' @param tolerance The tolerance within which a node date has to be removed from time = 0 
#' (i.e. the modern) to issue a warning that there are 'negative' node dates.

#' @return 
#' Returns a vector of length \code{Ntip(tree) + Nnode(tree)} which contains the dates for
#' all terminal tip nodes and internal nodes for the tree, in that order, as numbered in the \code{tree$edge}
#' matrix. These dates are always on a descending scale (i.e. time before present); see rootAge for how
#' the present time is determined. If rootAge is so defined that some nodes may occur later than
#' time = 0 units before present, this function may (confusingly) return negative dates and a 
#' warning message will be issued.

#' @seealso \code{\link{compareTimescaling}}

#' @author 
#' David W. Bapst, based on a function originally written by Graeme Lloyd.

#' @examples
#' #let's simulate some example data
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' #get the true time-sclaed tree
#' tree1 <- taxa2phylo(taxa)
#'
#' #now let's try dateNodes
#' dateNodes(tree1)
#' 
#' #let's ignore $root.time
#' dateNodes(tree1,rootAge=NULL)
#' 
#' #with the lengthy tip-label based labels
#'    #some of these will be hideously long
#' dateNodes(tree1,labelDates=TRUE)

#' @export
dateNodes<-function(tree,rootAge=tree$root.time,labelDates=FALSE,tolerance=0.001){
	#based on date.nodes by Graeme Lloyd, but using node.depth.edgelength
	#checks
	if(!inherits(tree,"phylo")){
		stop("tree must be of class 'phylo'")
		}
	#test that it has edge lengths
	if(is.null(tree$edge.length)){stop("tree does not appear to have edge lengths?")}
	nodeRelTimes<-node.depth.edgelength(tree)
	if(is.null(rootAge)){
		rootAge <- max(nodeRelTimes)
		message("Root age not given; treating tree as if latest tip was at modern day (time=0)")}
	res<-rootAge-nodeRelTimes
	if(any(res<(-tolerance))){
		message("Warning: Some dates are negative? rootAge may be incorrectly defined or you are using a time-scaling method that warps the tree, like aba or zbla.")}
	if(labelDates){
		names(res)<-sapply(Descendants(tree),function(x) paste0(sort(tree$tip.label[x]),collapse=" "))
	}else{
		names(res)<-1:(Ntip(tree)+Nnode(tree))
		}
	return(res)
	}
