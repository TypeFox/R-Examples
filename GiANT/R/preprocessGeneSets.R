##########################################################################
#formatGenesets
#
# dat:		Matrix of data with one row per gene.
# 			Only the row-(gene-)names are used.
# genesets: A list of gene sets (e.g. GO categories or KEGG pathways)
#
# The function removes all genes from the gene sets which are not part
# of the experiment (not in rownames(dat)). All names are set to lower 
# case.
##########################################################################
preprocessGs <- function(
		dat,
		geneSets){
	gs_sizes <- sapply(geneSets, length)

	allgenes <- tolower(rownames(dat))
	geneSets <- lapply(geneSets, tolower)
	new_geneSets <- lapply(geneSets, intersect, allgenes)

	return(new_geneSets)
}

##########################################################################
#filterGenesets
#
# genesets: A list of gene sets (e.g. GO categories or KEGG pathways)
# includedGenes:	A vector of gene names or a single gene name. All
#					gene sets without these genes are filtered out.
#					NULL if this type of filtering should not be applied. 
#					Also needed as 'start point' for filtering according
#					to interactions included in gene sets.
# minIntersectSize:	If not NULL, all gene sets with an intersect to
#					includedGenes smaller than this number are removed.
# adjMatrix:		An adjacence matrix showing the interactions within
#					the analyzed genes or a subset of them or NULL.
# steps:			If adjMatrix is not NULL also steps must be set to
#					the number of interaction steps which should be
#					added to 'includedGenes'. 'steps' = 1 means that all
#					genes which are, according to 'adjMatrix', directly
#					interacting with the initial gene(s) in
#					'includedGenes' are added to 'includedGenes'.
#
# Gene sets which do not fulfill the constraints given to this function
# are removed from 'genesets'. Possible constraints are a set of genes
# which have to be in each gene set, a minimal intersect to a 'coreset'
# or a set of genes interacting with each other (directly or with 
# intermediate steps).
##########################################################################
filterGeneSets <- function(
		geneSets,
		includedGenes = NULL,
		minIntersectSize = length(includedGenes),
		adjMatrix = NULL,
		steps = NULL){

	gs_size <- length(geneSets)
	tmp_genesets <- geneSets
	if(!is.null(adjMatrix)){
		if(is.null(includedGenes)){
			stop("If an adjacence matrix is given also a starting point (gene) must be given in 'includedGenes'.\n")
		}
		if(is.null(steps)){
			stop("If an adjacence matrix is given also a number of interaction steps must be given.\n")
		}
		if(is.null(rownames(adjMatrix)) || is.null(colnames(adjMatrix))){
			stop("Rows and cols of adjMatrix must be named.")
		}else if(sum(rownames(adjMatrix) != colnames(adjMatrix)) > 0){
			stop("adjMatrix is named incorrect. rownames(adjMatrix) != colnames(adjMatrix)")
		}

		recInteractions <- function(adjMatrix,set,step){
			if(step == 0){
				return(set)
			}else if(is.null(set)){
				return(set)
			}else{
				nms <- rownames(adjMatrix)
				new_set <- nms[rowSums(adjMatrix[,set, drop=FALSE]) > 0]
				return(c(set, recInteractions(adjMatrix, set = new_set, step = step-1)))
			}
		}

		interact <- recInteractions(adjMatrix = adjMatrix, set = includedGenes, step = steps)

		includedGenes <- unique(interact)
	}
	if(!is.null(includedGenes)){
		tmp <- sapply(tmp_genesets, function(x){
				length(intersect(x, includedGenes)) >= minIntersectSize
			})
		tmp_genesets <- tmp_genesets[tmp]
	}
	return(tmp_genesets)
}
