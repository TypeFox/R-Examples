## methods for computing various statistics of a tree's splits.

topNsplitvars <- function(tree,topn=0){
# Internal function that
# returns splits vars for topn depth splits. 0= root-only
	if(!inherits(tree, "itree")) stop("Not legitimate itree object")
	
	nn <- as.numeric(rownames(tree$frame))  #rowname of the frame is actually the nodenumber
	notleaf <- tree$frame[,1]!="<leaf>"
	depth <- floor(log(nn,base=2))[notleaf]

	tf <- tree$frame[notleaf,1] #get rid of the leafs[order(depth)]
	tf <- tf[order(depth)]  #order by depth
	depth <- depth[order(depth)]

	tf[(depth <= topn)]
}

splitstats <- function(tree,featlist=NULL){
# Externally available function.
# count of how many times each feature in featlist appears in the tree
# also total "inverse-depth"
# if featlist is NULL than it takes the features from the itree object.
     if(!inherits(tree, "itree"))
                stop("'tree' is not a legitimate itree object!")
                
	if(is.null(featlist)){
		featlist <- attr(tree$terms, "term.labels")
	}
	
	#compute non-leaf splitvars, depth, normalized inverse-depth
	ff  <- tree$frame
	nn <- as.numeric(rownames(ff))  #rowname of the frame is actually the nodenumber
	depth <- round(log(nn,base=2))
	depth <- depth[ff[,1]!="<leaf>"]
	md <- max(depth)  #lowest non-root
	inv.depth <- 1/(depth+1)   #higher = more important. 0 = not in the tree.
	inv.depth <- inv.depth/sum(inv.depth)

	#nodesize
	nodesize <- ff$n[ff[,1]!="<leaf>"]
	nodesize <- nodesize/sum(nodesize) 

	#splitvars
	splitvar <- as.character(ff[ff[,1]!="<leaf>",1]) #get rid of the leaves
	splitvars.unique <- unique(splitvar)
	names(inv.depth) <- names(nodesize) <-  splitvar
	
	#OK, now record stuff about which variables are used, where they occur, etc.
	statmat <- matrix(0,nrow=length(featlist),ncol=6)

	#first a useful fcn telling us where to record a given number.
	where <- function(varlist){
		apply(as.matrix(varlist,ncol=1),1,FUN=function(split){which(split==featlist)})
	}

	#counts
	ct <- table(splitvar)
	statmat[where(names(ct)),2] <- ct
	
	#total inverse depth for each variable
	str <- apply(as.matrix(splitvars.unique,ncol=1),1,FUN=function(splitvar){sum(inv.depth[names(inv.depth)==splitvar])})
	statmat[where(splitvars.unique),3] <- str
	
	#normalized sum of nodes sizes that each variable split. 
	total.ns <- apply(as.matrix(splitvars.unique,ncol=1),1,FUN=function(splitvar){sum(nodesize[names(nodesize)==splitvar])})			
	statmat[where(splitvars.unique),4] <- total.ns

	#root and depth=1
	topsplits <- topNsplitvars(tree,1)
	ww <- where(topsplits)
	statmat[ww[1],5] <- 1
	statmat[ww[2:length(ww)],6] <- 1
	
	statmat <- as.data.frame(statmat)
	statmat[,1] <- as.factor(featlist)
	colnames(statmat) <- c("var","split.ct","total.inv.depth","total.node.size","isroot","isdepth1")
	return(statmat)
}

