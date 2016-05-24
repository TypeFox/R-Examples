loadTargetsToGenes <- function(targToGenesFile, targRegexpReplace=list(c("\\.\\.", "-")), geneDelim="\\|", targColumn=1, geneColumn=3, CHROMOSOMES_START_WITH_CHR=FALSE) {
	if (!file.exists(targToGenesFile)) {
		stop(paste("Cannot find file '", targToGenesFile, "'", sep=""))
	}

	if (!CHROMOSOMES_START_WITH_CHR) {
		if (is.null(targRegexpReplace)) {
			targRegexpReplace = list()
		}
		targRegexpReplace[[length(targRegexpReplace) + 1]] = c("^chr", "")
	}

	g = read.table(targToGenesFile, header=FALSE, check.names=FALSE, row.names=targColumn, stringsAsFactors=FALSE)

	if (geneColumn > targColumn) { # Correct for the targColumn being used as the row names [if it comes after targColumn]:
		geneColumn = geneColumn - 1
	}
	targetsToGenes = g[, geneColumn]

	names(targetsToGenes) = rownames(g)
	if (!is.null(targRegexpReplace)) {
		for (i in 1:length(targRegexpReplace)) {
			names(targetsToGenes) = sub(targRegexpReplace[[i]][1], targRegexpReplace[[i]][2], names(targetsToGenes), perl=TRUE)
		}
	}

	targetsToGenes = strsplit(targetsToGenes, geneDelim, perl=TRUE)

	return(targetsToGenes)
}
