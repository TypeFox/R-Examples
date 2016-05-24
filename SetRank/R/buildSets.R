# Project: SetRank
# 
# Author: cesim
###############################################################################
utils::globalVariables(c("keys", "select"))

uniqueCount <- function(x) {
	if (class(x) == "factor") length(levels(x)) else length(unique(x))
}

expandWithTermOffspring <- function(subTable, tableSplit, offspringList) {
	termID = unique(as.character(subTable$termID))
	termName = unique(as.character(subTable$termName))
	dbName = unique(as.character(subTable$dbName))
	offspring = offspringList[[termID]] %i% names(tableSplit)
	extension = do.call( rbind, lapply(offspring, 
					function(x) data.frame(geneID = tableSplit[[x]]$geneID, 
								termID = termID, termName = termName, 
								dbName = dbName, stringsAsFactors = FALSE)))
	expandedTable = rbind(subTable, extension)
	expandedTable[!is.na(expandedTable$geneID),]
}

#' Create a gene set collection
#' 
#' Builds an object containing the collection of all gene sets to be used 
#' by the \code{\link{setRankAnalysis}} function.
#' 
#' @section Execution time:
#' This function typically takes some time to execute as it pre-calculates all
#' significant intersections between pairs of gene sets in the collection. An 
#' intersection between two gene sets is considered significant if it contains 
#' more elements than expected by chance, given the sizes of both sets. 
#' Computation time can be sped up dramatically by running this function on
#' multiple CPU-cores. To do so, simply set the \code{mc.cores} option to the
#' desired number of cores to use, like so: \code{options("mc.cores=4")}
#' Performing this calculation beforehand allows to re-use the same 
#' setCollection object for different analysis. It is therefore recommended to 
#' separate the creation of the setCollection object and the actual analysis in
#' different scripts. Once the collection is created, it can be stored on disk 
#' using the \code{save} command. The analysis script can then load the 
#' collection using the \code{load} command.
#' 
#' @section Creation of custom annotation tables:
#' \describe{
#'   \item{geneID}{The gene identifier. Can be any type of identifier, but one
#'     must make sure that all annotation frames passed to 
#'     \code{buildSetCollection} use the same identifier. As the packages 
#'     created by the GeneSets package use Entrez Gene identifiers, it
#'     is best to use these in your own annotation frames as well. Also, make 
#'     sure the identifiers as passed as character and not as integer values.}
#'   \item{termID}{Pathway identifier. Make sure each pathway identifier is 
#'     unique across all pathway databases used. You can do this by prefixing 
#'     the IDs with a namespace identifier like \code{"REACTOME:"}.}
#'   \item{termName}{Name of the pathway. A string describing the pathway, e.g.
#'     "negative regulation of sterol metabolism"}
#'   \item{description}{Pathway description. A longer description of the 
#'     pathway. This field can be a full paragraph describing what this pathway
#'     does.}
#'   \item{dbName}{A short string given the name of the pathway database used for
#'     the annotation. E.g. "KEGG".}
#' }
#' 
#' @param \ldots One or more data frame objects containing the annotation of 
#' genes with pathway identifiers and descriptions. The idea is to provide one
#' data frame per pathway database. Several gene set databases are provided in
#' the organism-specific GeneSets packages. Alternatively, you can specify your
#' own annotation tables. See the \emph{Details} section for more information.
#' @param referenceSet Optional but very strongly, recommended. A vector of 
#' geneIDs specifying the background gene set against which to test for 
#' over-representation of genesets. The default is to use all genes present in
#' the supplied gene annotation tables. However, many experiments are 
#' intrinsically biased for certain pathways e.g. because they only contain 
#' samples from a specific tissue. Supplying a suitable reference set will 
#' remove this bias. See the vignette for more details.
#' @param maxSetSize The maximum number of genes in a gene set. Any gene sets
#' with more genes will not be considered during the analysis.
#' 
#' @return A gene set collection which is a list object containing the following
#' fields:
#' \itemize{
#'   \item{maxSetSize}{The maximum set size applied when constructing the
#'     collection.}
#'   \item{referenceSet}{A vector listing all gene IDS that are part of the
#'     reference.}
#'   \item{sets}{A list of vectors. The list names are the pathway IDs as supplied
#'     in the \code{termID} column of the annotation frame(s) supplied..}  Each
#'     vector contains all geneIDs of the gene set and has three attributes set: 
#'     \code{ID}, \code{name}, and \code{db} which correspond respectively to the
#'     \code{termID}, \code{termName}, and \code{dbName} fields of the annotation
#'     frame.
#'   \item{g}{The size of the reference set.}
#'   \item{bigSets}{A list of pathway IDs of gene sets with sizes bigger than the
#'     specified maximum set size.}
#'   \item{intersection.p.cutoff}{The p-value cutoff used to determine which
#'     intersections of pairs of gene sets (see \emph{Details}) are significant.}
#'   \item{intersections}{A data frame listing all significant intersections
#'     together with the p-value.}  
#' }
#' @examples
#' options(mc.cores=1)
#' referenceSet = sprintf("gene_%02d", 1:50)
#' geneSets = lapply(1:9, function(i) sample(referenceSet[((i-1)*5):((i+1)*5)], 5))
#' annotationTable = data.frame(termID=sprintf("set_%02d", rep(1:9, each=5)), 
#'         geneID=unlist(geneSets),
#'         termName = sprintf("dummy gene set %d", rep(1:9, each=5)),
#'         dbName = "dummyDB",
#'         description = "A dummy gene set DB for testing purposes")
#' collection = buildSetCollection(annotationTable, referenceSet=referenceSet)
#' @author Cedric Simillion
#' @import parallel
#' @export
buildSetCollection <- function(..., referenceSet = NULL, maxSetSize = 500) {
	annotationTable = do.call(rbind, list(...))
	if	(!is.null(referenceSet)) {
		annotationTable = 
				annotationTable[annotationTable$geneID %in% referenceSet,]
		annotationTable$termID = factor(annotationTable$termID)
	} else {
		referenceSet = unique(referenceSet)
		referenceSet = unique(as.character(annotationTable$geneID))
	}
	collection = list(maxSetSize = maxSetSize, referenceSet=referenceSet)
	collection$sets = if (nrow(annotationTable) == 0) list() else
				by(annotationTable, annotationTable$termID, createSet, 
						simplify=FALSE)
	collection$sets[sapply(collection$sets, is.null)] = NULL
	collection$g =  length(referenceSet)
	collection$bigSets = sapply(collection$sets, length) > maxSetSize
	message(uniqueCount(annotationTable$dbName), " gene set DBs, ", 
			length(collection$sets), " initial gene sets, ", 
			length(collection$sets) - length(which(collection$bigSets)), 
			" sets remaining and ", collection$g, " genes in collection")
	collection$intersection.p.cutoff = 0.01
	cluster = if (.Platform$OS.type == "windows") 
            makePSOCKcluster(getOption("mc.cores")) else makeForkCluster()
	collection$intersections = getSignificantIntersections(collection$sets, 
			annotationTable, collection$g, collection$intersection.p.cutoff,
			collection$bigSets, cluster)
	message("Pre-calculating critical Fisher-test values...")
	collection$iMatrix = fisherCriticalValues(collection$g, maxSetSize, 0.05,
			cluster)
	stopCluster(cluster)
	collection
}

#' @import data.table
getSignificantIntersections <- function(collectionSets, annotationTable, g, 
		pValueCutoff, bigSets, cluster) {
	setIDs = names(collectionSets[!bigSets])
	setIndicesPerGene = by(annotationTable, annotationTable$geneID, 
			function(x) which(setIDs %in% as.character(x$termID)))
	setCount = unlist(lapply(setIndicesPerGene, length))
	geneOrder = names(sort(setCount[(setCount > 1)], decreasing=TRUE))
	intersectionsPerGene = lapply(setIndicesPerGene[geneOrder],
			function(x) apply(t(combn(x, 2)), 1, pack, length(setIDs)))
	intersectionIndices = unlist(intersectionsPerGene, use.names=FALSE)
	intersectionIndices = unique(intersectionIndices)
	message(length(intersectionIndices), " intersections to test...", 
			appendLF=FALSE)
	intersectionPValues = parSapply(cluster, intersectionIndices, 
			getIntersectionPValue, collectionSets[!bigSets], g)
	significantIndices = which(intersectionPValues <= pValueCutoff)
	pValueFrame = rbindlist(parLapply(cluster, 
					intersectionIndices[significantIndices], unpack, setIDs))
	if (nrow(pValueFrame) > 0) {
		pValueFrame$pValue = intersectionPValues[significantIndices]
		message(nrow(pValueFrame), " intersections significant")
	} else {
		pValueFrame = data.frame(setA=c(),setB=c())
	}
	pValueFrame
}

createSet <- function(x) {
	geneSet = unique(as.character(x$geneID))
	attr(geneSet, "ID") <- 
			unique(as.character(x$termID))
	attr(geneSet, "name") <- 
			unique(as.character(x$termName))
	attr(geneSet, "db") <- 
			unique(as.character(x$dbName))
	attr(geneSet, "description") <- 
			unique(as.character(x$description))
	geneSet
}


pack <- function(indexPair, n) {
	if (indexPair[1] > n || indexPair[2] > n) stop("Index higher than n")
	if (indexPair[1] > indexPair[2]) {
		a = indexPair[2]
		b = indexPair[1]
	} else {
		a = indexPair[1]
		b = indexPair[2]
	}
	(a-1)*(n-1)+(b-1)
}

unpack <- function(packed, setIDs) {
	n = length(setIDs)
	n_ = n-1
	a = ceiling(packed/n_)
	r = (packed %% n_)
	b = if (r == 0) n else r+1
	if (a >= b) stop("Invalid packed value")
	data.frame(setA=setIDs[a], setB=setIDs[b], stringsAsFactors=FALSE)
}

getIntersectionPValue <-function(intersectionIndex, setCollection, g) {
	setIndexPair = unpack(intersectionIndex, names(setCollection))
	setA = setCollection[[setIndexPair$setA]]
	setB = setCollection[[setIndexPair$setB]]
	i = length(setA %i% setB)
	m = length(setA)
	n = length(setB)
	return(intersectionTest(g, m, n, i)$p.value)
}	

intersectionTest <- function(g, m, n, i) {
	fisher.test(rbind(c(i, n-i), c(m, g-(m+n-i))), alternative="greater")
}

fisherCriticalValues <- function(g,maxSize,criticalP, cluster) {
	do.call(rbind, parLapply(cluster, 1:g, function(s) {
						sapply(1:maxSize, minimalI, s, g, criticalP)
					}))
}

minimalI <- function(m,s,g, maximalP) {
	iVector=min(s,m):0
	pValues = rev(cumsum(dhyper(s-iVector,g-m,m,s)))
	minimalI = which(pValues < maximalP)[1] - 1
}

#' Create a function to convert gene or protein IDs
#' 
#' Creates a function based on an \code{AnnotationDb} package. This 
#' package accepts a vector of input IDs and returns a vector of output IDs. 
#' If an input ID cannot be mapped to an output ID, to output vector will be 
#' one element shorter. This behaviour can be changed by setting the additional
#' \code{na.rm} argument to \code{FALSE}. Likewise, if an input ID maps to 
#' multiple output IDs, the output vector will contain all of the latter. If 
#' you really need the output vector to have the same length as the input
#' vector, you can set the \code{drop.ambiguous} argument to \code{TRUE}
#' 
#' @param annotationPackageName The name of the \code{AnnotationDb} package 
#' that will be used to create the conversion function. The package will be 
#' loaded automically if necessary.
#' @param from The ID type to convert from. This should be one of the available
#' keytypes in the \code{AnnotationDb} package. Use the \code{keytypes}
#' function to find out which keytypes can be used.
#' @param to The ID type to convert to. This should be one of the available 
#' columns in the \code{AnnotationDb} package. Use the \code{cols}
#' function to find out which column names can be used.
#' 
#' @return A function which takes a vector of input IDs as single argument and
#' returns another vector with the converted IDs.
#' 
#' @author Cedric Simillion
#' @export
createIDConverter <- function(annotationPackageName, from, to) {
	require(annotationPackageName, character.only=TRUE)
	keySet = keys(eval(parse(text=annotationPackageName)), from)
	conversionTable = select(eval(parse(text=annotationPackageName)), 
			keys=keySet, keytype=from, columns=to)
	conversion = as.list(by(conversionTable, as.factor(conversionTable[[from]]), 
					function(t) unique(t[[to]]), simplify=FALSE))
	outputFunction = function(x, na.rm=TRUE, drop.ambiguous=FALSE) {
		knownX = unique(x)
		outputList = as.list(rep(NA, length(knownX)))
		names(outputList) = knownX
		knownX = intersect(knownX, names(conversion))
		if (drop.ambiguous) {
			ambiguous = sapply(conversion, function(x) length(x) > 1)
			knownX = setdiff(knownX, names(conversion[ambiguous]))
		}
		outputList[knownX] = conversion[knownX]
		output = unlist(outputList[x], use.names=FALSE)
		if (na.rm) {
			output = output[!is.na(output)]
		}
		output
	} 
	return(outputFunction)
}


