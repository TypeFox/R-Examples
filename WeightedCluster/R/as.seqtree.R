
###########################
## disstree main function
###########################
as.seqtree <- function(object, seqdata, diss, weighted=TRUE,...){
	UseMethod("as.seqtree")
}
as.seqtree.twins <- function(object, seqdata, diss, weighted=TRUE, ncluster, ...) {
	return(as.seqtree.hclust(object, seqdata=seqdata, diss=diss, weighted=weighted, ncluster=ncluster,...))
}
as.seqtree.hclust <- function(object, seqdata, diss, weighted=TRUE, ncluster, ...) {
	pred <- data.frame(Split2=factor(cutree(object, 2)))
	for(p in 3:ncluster){
		pred[, paste("Split", p, sep="")] <- factor(cutree(object, p))
	}
	object <- pred
	return(as.seqtree.default(object, seqdata=seqdata, diss=diss, weighted=weighted, ...))
}

as.seqtree.default <- function(object, seqdata, diss, weighted=TRUE, ...) {
	predictor <- object
	ncluster <- ncol(object)+1
	if(ncluster<2){
		stop(" [!] ncluster should be bigger than 2")
	}
	if (inherits(diss, "dist")) {
		diss <- as.matrix(diss)
 	}
	## Model matrix from forumla
	nobs= nrow(diss)
	
	## Allow integer weights for replicates
	if(weighted & !is.null(attr(seqdata,"weights") )){
		weights <- attr(seqdata,"weights") 
	} 
	else {
		weights <- as.double(rep(1,nobs))
	}
	pop <- sum(weights)
#	TraMineR:::.localstuffDissTree$DTNnodeCounter <- as.integer(1)
	
	vardis <- dissvar(diss, weights=weights)
	
	as.seqtreeDTNBuildNode <- function(ind, vardis, depth, current) {
		node <- TraMineRInternalNodeInit(ind=ind, vardis=vardis, depth=depth, dmat=diss, weights=weights)
		node$info$splitschedule <- depth
		SCtot <- vardis*node$info$n
		SCres <- SCtot
		if(current>ncluster ||length(ind)==1){
			return(node)
		}
		## print(SCtot)
		#varnames <- colnames(pred)
		
		for (p in current:ncluster) {
			clust <- predictor[ind, p-1]
			SplitLabels <- unique(clust)
			if(length(SplitLabels)==2){
				#print(p)
				#print(ncluster)
				bestSpl <- list()
				bestSpl$variable <- clust == SplitLabels[1]
				# print(ind)
				# print(as.integer(ind[bestSpl$variable]))
				# print(diss)
				lSCres <- .Call("tmrWeightedInertiaDist", diss, as.integer(nrow(diss)), 
					as.integer(FALSE), as.integer(ind[bestSpl$variable]), as.double(weights), 
					as.integer(FALSE), package="TraMineR")
				rSCres <- .Call("tmrWeightedInertiaDist", diss, as.integer(nrow(diss)), 
					as.integer(FALSE), as.integer(ind[!bestSpl$variable]), as.double(weights), 
					as.integer(FALSE), package="TraMineR")
				info <- list(
						lpop=sum(weights[ind[bestSpl$variable]]),
						rpop=sum(weights[ind[!bestSpl$variable]]),
						SCres=lSCres+rSCres
				)
				info$lvar=lSCres/info$lpop
				info$rvar=rSCres/info$rpop
				bestSpl$spl <- TraMineRInternalSplitInit(p-1, index = 1:2,  
					 prob = c(info$lpop, info$rpop)/node$info$n, info = info, labels=SplitLabels)
				SCres <- bestSpl$spl$info$SCres
				node$split <- bestSpl$spl
				node$split$info$R2 <- 1-(SCres/SCtot)
				node$surrogates <- bestSpl$sur
				
				#print(bestSpl)
				left <- as.seqtreeDTNBuildNode(ind=ind[bestSpl$variable], vardis=bestSpl$spl$info$lvar, depth=p, current=p+1)
				right <- as.seqtreeDTNBuildNode(ind=ind[!bestSpl$variable], vardis=bestSpl$spl$info$rvar, depth=p, current=p+1)
				node$kids <- list(left, right)
				## We have found the split, so leave the loop
				return(node)
			}
		}
		## Maximum depth reached
		return(node)
	}
	
	root <- as.seqtreeDTNBuildNode(ind=1:nobs, vardis=vardis, depth=1, current=2)
	#print(root)
	tree <- list()
	tree$fitted <- data.frame(disstreeleaf(root))
	names(tree$fitted) <- "(fitted)"
	tree$info <- list(method="disstree", n=pop, parameters= list(minSize=1, maxdepth=ncluster, R=0, pval=1), object=seqdata, weight.permutation="diss")
	if(!weighted) {
		tree$info$adjustment <- dissassoc(diss, tree$fitted[,1], R=0, weights=NULL)
	}
	else {
		tree$info$adjustment <- dissassoc(diss, tree$fitted[,1], R=0, weights=weights, weight.permutation="diss")
	}
	tree$data <- predictor
	tree$terms <- NULL
	tree$weights <- weights
	##tree <- party(root, data=predictor, fitted =fitted, terms = terms(formula.call),  info = info)
	tree$root <- root
	
	class(tree) <- c("seqtree", "disstree", class(tree))
	return(tree)
}

