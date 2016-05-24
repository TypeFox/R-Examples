##############################################################
#
# GUniFrac2: Generalized UniFrac distances for comparing microbial
#						communities. (Modified version)
# Chong Wu (wuxx0845@umn.edu)
#
# The original version is contributed by
# Jun Chen (chenjun@mail.med.upenn.edu)
# Jun 12, 2015
#
###############################################################

rdirichlet<-function(n,alpha)
## generate n random deviates from the Dirichlet function with shape
## parameters alpha
{
    l<-length(alpha);
    x<-matrix(rgamma(l*n,alpha),ncol=l,byrow=TRUE);
    sm<-x%*%rep(1,l);
    x/as.vector(sm);
}

GUniFrac_cum <- function (otu.tab, tree) {
	# Calculate Generalized UniFrac distances. Unweighted and 
	# Variance-adjusted UniFrac distances will also be returned.
	#	
	# Args:
	#		otu.tab: OTU count table, row - n sample, column - q OTU
	#		tree: rooted phylogenetic tree of R class "phylo"
	#		alpha: parameter controlling weight on abundant lineages
	#
	# Returns:
	# 	unifracs: three dimensional array containing the generalized 
	#							UniFrac distances, unweighted UniFrac distance and 
	#							variance adjusted UniFrac distances. 
	#
	if (!is.rooted(tree)) stop("Rooted phylogenetic tree required!")
	
	# Convert into proportions
	otu.tab <- as.matrix(otu.tab)
	row.sum <- rowSums(otu.tab)
	otu.tab <- otu.tab / row.sum
	n <- nrow(otu.tab)
	
	# Construct the returning array
	if (is.null(rownames(otu.tab))) {
		rownames(otu.tab) <- paste("comm", 1:n, sep="_")
	}

	# Check OTU name consistency
	if (sum(!(colnames(otu.tab) %in% tree$tip.label)) != 0) {
		stop("The OTU table contains unknown OTUs! OTU names
					in the OTU table and the tree should match!" )
	}
	
	# Get the subtree if tree contains more OTUs
	absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
	if (length(absent) != 0) {
		tree <- drop.tip(tree, absent)
		warning("The tree has more OTU than the OTU table!")
	}
	
	# Reorder the otu.tab matrix if the OTU orders are different
	tip.label <- tree$tip.label
	otu.tab <- otu.tab[, tip.label]
	
    
	ntip <- length(tip.label)
	nbr <- nrow(tree$edge)	
	edge <- tree$edge
	edge2 <- edge[, 2]
    br.len <- tree$edge.length  # branch length
	

    #  Accumulate OTU proportions up the tree
	cum <- matrix(0, nbr, n)							# Branch abundance matrix
	for (i in 1:ntip) {
		tip.loc <- which(edge2 == i)
		cum[tip.loc, ] <- cum[tip.loc, ] + otu.tab[, i]	
		node <- edge[tip.loc, 1]						# Assume the direction of edge 
		node.loc <- which(edge2 == node)
		while (length(node.loc)) {
			cum[node.loc, ] <- cum[node.loc, ] + otu.tab[, i]		
			node <- edge[node.loc, 1]
			node.loc <- which(edge2 == node)
		}
	}

	out = list(cum = cum, br.len= br.len)
	return(out)
}


GUniFrac <-function (otu.tab, tree, alpha= c(0,0.5,1)) {
    alpha = matrix(alpha,1,length(alpha))
    alpha2 =matrix(c(0,0.5,1),1,3)
    
    GuniF.cum = GUniFrac_cum(otu.tab,tree)
    cum = GuniF.cum$cum
    br.len = GuniF.cum$br.len
    br.len = as.matrix(br.len)
    
    if (identical(alpha,alpha2)) # if true, use a faster version
    {
        GUniF = GUniFracCpp(cum,br.len)
        d0 = GUniF$d0
        dw <- GUniF$d1  # Weighted UniFrac
        d5 <- GUniF$d5
    
        output = list(d0 = d0,dw=dw, d5=d5)
    } else {
        GUniF = GUniFracCpp2(cum,br.len, alpha)
        output = list(GUniF = GUniF, alpha = alpha)
    }
    output
}





