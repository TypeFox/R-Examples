####################################################
#
#		samplingProbs <- function(...)
#
#		tree = object of class phylo
#		cladeTable = a dataframe with 1 column of species names and a second column of group assignment
#			Must either be a table for all species in tree, or a table of greater species richness, including those species in the tree.
#		cladeRichness = either NULL or a vector of species counts, named by group names.
#		globalSampling = percent sampling of the backbone of the phylogeny
#		output = path + output file name (.txt)
#		writeToDisk = boolean, should the table be written to disk, defaults to TRUE


samplingProbs <- function(tree, cladeTable, cladeRichness = NULL, globalSampling, output, writeToDisk = TRUE) {
	
	if (length(intersect(tree$tip.label,cladeTable[,1])) != length(tree$tip.label)) {
		stop("Not all species from tree are in cladeTable.");
	}
	
	if (nrow(cladeTable) == length(tree$tip.label)) {
		if (is.null(cladeRichness)) {
			stop("If cladeTable only contains species from tree, then cladeRichness must be provided.");
		}
	}
	
	if (nrow(cladeTable) > length(tree$tip.label)) {
		if (!is.null(cladeRichness)) {
			warning("cladeTable contains more species than in tree, so cladeRichness vector will be ignored.");
		}
	}
	
	if (!is.null(cladeRichness)) {
		if (length(cladeRichness) != length(unique(cladeTable[,2]))) {
			stop("The cladeRichness vector must contain the same number of clades as are described in the cladeTable.");
		}
	}
	
	if (!is.vector(cladeRichness) & !is.null(cladeRichness)) {
		stop("Error: cladeRichness must either be NULL or an integer vector named with clade names.");
	}
	
	if (ncol(cladeTable) > 2) {
		stop("cladeTable must contain 2 columns: one of species, and one of clade assignment.");
	}
	
	if (is.matrix(cladeTable)) {
		cladeTable <- as.data.frame(cladeTable, stringsAsFactors=FALSE);
	}
	
	if (nrow(cladeTable) > length(tree$tip.label)) {
		probs <- as.data.frame(matrix(nrow=length(tree$tip.label),ncol=3));
		colnames(probs) <- c('sp','clade','prob');
		for (i in 1:length(tree$tip.label)) {
			probs[i,1] <- tree$tip.label[i];
			clade <- cladeTable[cladeTable[,1] == tree$tip.label[i],2];
			inTree <- intersect(cladeTable[cladeTable[,2] == clade,1],tree$tip.label);
			probs[i,2] <- clade;
			probs[i,3] <- length(inTree) / length(cladeTable[cladeTable[,2] == clade,1]);
		}
		probs <- rbind(c(globalSampling,'','',''),probs);
	}
	
	if (nrow(cladeTable) == length(tree$tip.label) & !is.null(cladeRichness)) {
		probs <- as.data.frame(matrix(nrow=length(tree$tip.label),ncol=3));
		colnames(probs) <- c('sp','clade','prob');
		for (i in 1:length(tree$tip.label)) {
			probs[i,1] <- tree$tip.label[i];
			clade <- cladeTable[cladeTable[,1] == tree$tip.label[i],2];
			probs[i,2] <- clade;
			probs[i,3] <- nrow(cladeTable[cladeTable[,2] == clade,]) / cladeRichness[clade];
		}
		probs <- rbind(c(globalSampling,'','',''),probs);
	}
	if (writeToDisk) {
		write.table(probs, file=output, quote=F, col.names=F, row.names=F, sep='\t');
	} else {
		return(probs);
	}
}








