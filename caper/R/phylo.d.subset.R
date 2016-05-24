phylo.d.subset <- function(data, phy, names.col, binvar, permut=1000, 
	                       rnd.bias=NULL, min.tips=1, max.tips=length(data$phy$tip.label), 
	                       min.nodes=1, max.nodes=data$phy$Nnode, verbose=FALSE) {

    # - test to see if there is a comparative data object and if not then
    #   retrofit the remaining arguments into a comparative data object.
	if(! missing(data)){
		if(! inherits(data, 'comparative.data')){
			if(missing(names.col)) stop('names column is missing')
			names.col <- deparse(substitute(names.col))
			data <- caicStyleArgs(data=data, phy=phy, names.col=names.col)
		}
	}
	
	# look for binary variable
	binvar <- deparse(substitute(binvar))
    bininds <- match(binvar, names(data$data))
    if (is.na(bininds)) (stop("'", binvar, "' is not a variable in data."))

	# get the variable out and do a general test for binarity
	ds <- data$data[ ,bininds]
	if(length(unique(ds)) != 2) stop("'", binvar, "' doesn't contain two states.")
	if(any(is.na(ds))) stop("'", binvar, "' contains missing values.")
	
	# check for a number
    if (!is.numeric(permut)) (stop("'", permut, "' is not numeric.")) 
	
	# look for probaility weights argument and get its value if found
	if(missing(rnd.bias)) rnd.bias<-NULL else {
    	rnd.bias <- deparse(substitute(rnd.bias))
    	rnd.ind <- match(rnd.bias, names(data$data))
    	if (is.na(rnd.ind)) (stop("'", rnd.bias, "' is not a variable in data."))
		rnd.bias<-data$data[ ,rnd.bias]}
	
	# make the subtrees
	phy.subtrees <- subtrees(data$phy)
	if(verbose) phy.max.subtrees <- length(phy.subtrees)
	
	# filter them
	phy.subtrees <- phy.subtrees[sapply(phy.subtrees, function(x) max.tips >= length(x$tip.label) & length(x$tip.label) >= min.tips & x$Nnode <= max.nodes & min.nodes <= x$Nnode)]
	
	# prepare output
	output.raw <- vector(mode="list", length(phy.subtrees))
	output.D <- output.depth <- output.P0 <- output.P1 <- output.tips <- output.nodes <- output.bin.freq <- numeric(length(phy.subtrees))
	
	##TO-DO:
	#  - make a 'clean' version where the raw output isn't saved (when used on a big phylogeny this could be big)
	#  - allow filtering according to presence:absence
	#  - record node labels/tip labels within each subset in a convenient format
	
	# run phylo.d on each subtree
	if(verbose) cat("\nCalculating D values for ", length(phy.subtrees), " out of a possible ", phy.max.subtrees, "clades")
	for(i in seq(along=phy.subtrees)){
		if(verbose) cat(".")
		# make a temporary comparative.data.frame and check it's worth passing through phylo.d - if not, fill it with NAs
		# - assuming most users would want to know about subsets that fail this test---I would!
		t.data <- data[rownames(data$data) %in% phy.subtrees[[i]]$tip.label,]
		if(length(unique(t.data$data[,bininds])) != 2){
			output.raw[[i]] <- output.depth[i] <- output.P0[i] <- output.P1[i] <- NA
			next}
		# get the raw output
		output.raw[[i]] <- eval(substitute(phylo.d(t.data, binvar=XXX), list(XXX=as.name(names(data$data)[bininds]))))
		# get the summaries
		output.D[i] <- output.raw[[i]]$DEstimate
		output.P0[i] <- output.raw[[i]]$Pval0
		output.P1[i] <- output.raw[[i]]$Pval1
		clade.mat <- clade.matrix(phy.subtrees[[i]])
		output.depth[i] <- sum(clade.mat$edge.length[as.logical(clade.mat$clade.matrix[,1])])
		output.tips[i] <- length(phy.subtrees[[i]]$tip.label)
		output.nodes[i] <- phy.subtrees[[i]]$Nnode
		output.bin.freq[i] <- sum(t.data$data[,bininds])
		}
	
	output <- list(raw=output.raw, DEstimate=output.D, Pval1=output.P1, Pval0=output.P0, phy.depth=output.depth, tips=output.tips, nodes=output.nodes, bin.freq=output.bin.freq)
	class(output)<-'phylo.d.subset'
	return(output)
}

print.phylo.d.subset <- function(x, ...){
    summary(x)
}

##############################################
##TO-DO:
##MAKE THIS MORE PLEASANT TO INTERPRET
##ADD IN PLOTS OF D vs. PHYLOGENETIC DEPTH
##############################################

summary.phylo.d.subset <- function(object, ...){
    cat('\nCalculation of D statistic for the phylogenetic structure of a binary variable')
    cat('\n...across multiple clades within a phylogeny\n')
    cat('\n  Data : ', object$raw[[1]]$data$data.name)
    cat('\n  Binary variable : ', object$raw[[1]]$binvar)
    cat('\n  Phylogeny : ', object$raw[[1]]$data$phy.name)
    cat('\n  Number of permutations : ', object$raw[[1]]$nPermut)
    
    cat("\n\nEstimated D values: \n")
	t<-format(round(object$DEstimate, digits=2), trim=TRUE)
	cat(format(object$DEstimate, width=nchar(t)))
	if(is.null(object$raw[[1]]$rnd.bias)) cat("\nProbabilities of E(D) resulting from no (random) phylogenetic structure : \n") else cat("\nProbability of E(D) resulting from no (biased random) phylogenetic structure : \n")
	cat(format(object$Pval0, width=nchar(t)))
    cat("\nProbabilities of E(D) resulting from Brownian phylogenetic structure    : \n")
	cat(format(object$Pval1, width=nchar(t)))
	cat("\nAges of clades                                                        : \n")
	cat(format(object$phy.depth, width=nchar(t)))
    cat("\n\n")
}