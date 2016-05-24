# Project: SetRank
# 
# Author: cesim
###############################################################################
#' Advanced gene set enrichment analysis.
#' 
#' Performs advanced gene set enrichment analysis on a set of genes.
#' 
#' @param geneIDs A vector containing the set of gene IDs to test for gene set 
#' enrichment. This is typically the list of significant genes returned by the 
#' analysis of an omics dataset.
#' @param setCollection A gene set collection object, generated with the
#' \code{\link{buildSetCollection}} function.
#' @param use.ranks Logical value indicating if the \code{geneIDs} vector is in
#' ranked order or not. When \code{TRUE}, a ranked analysis will be performed.
#' @param setPCutoff The p-value cutoff to be used to consider a gene set
#' significant. Recommended value: 0.01
#' @param fdrCutoff The cutoff to be applied on the corrected p-value after
#' false-positive sets have been removed.
#' @param delete A flag indicating if non-significant gene sets should be 
#' deleted, which is the behaviour that you want. It is best to ignore this
#' argument, it was only added for debugging purposes.
#' 
#' @return An igraph object. Use the igraph \code{\link{get.data.frame}} 
#' function to get a data frame with all the significant gene sets.
#' @examples
#' options(mc.cores=1)
#' reference = sprintf("gene_%03d", 1:50)
#' geneSets = lapply(1:9, function(i) sample(reference[((i-1)*5):((i+1)*5)], 5))
#' annotationTable = data.frame(termID=sprintf("set_%02d", rep(1:9, each=5)), 
#'         geneID=unlist(geneSets),
#'         termName = sprintf("dummy gene setet %d", rep(1:9, each=5)),
#'         dbName = "dummyDB",
#'         description = "A dummy gene set DB for testing purposes")
#' collection = buildSetCollection(annotationTable, referenceSet=reference)
#' genes = reference[sample(c(TRUE, FALSE), 50, TRUE)]
#' network = setRankAnalysis(genes, collection, TRUE)
#'
#' @author Cedric Simillion
#' @import igraph
#' @import parallel
#' @export 
setRankAnalysis <- function(geneIDs, setCollection, use.ranks = TRUE,
		setPCutoff = 0.01, fdrCutoff = 0.05, delete = TRUE) {
	testSet = if (use.ranks) as.RankedTestSet(geneIDs, setCollection) else
				as.UnrankedTestSet(geneIDs, setCollection)
	message(Sys.time(), " - calculating primary set p-values")
	pValues = getPrimarySetPValues(testSet, setCollection)
	edgeTable = buildEdgeTable(testSet, setCollection, pValues,	setPCutoff)
	setNet = buildSetNet(edgeTable, setCollection, testSet, pValues, 
			setPCutoff)
	initialNodes = vcount(setNet)
	message(Sys.time(), " - 1st round of node removal")
	toDelete = getNodesToDelete(edgeTable)
	if (length(toDelete) > 0 && delete) setNet = setNet - toDelete
	message(Sys.time(), " - calculating SetRank values")
	setNet = calculateSetRank(setNet)
	message(Sys.time(), " - secondary delete")
	if (delete) {
		setNet = sinkDelete(setNet, setCollection, setPCutoff, testSet)
	}
	message(Sys.time(), " discarded ", initialNodes-vcount(setNet), " out of ", 
			initialNodes, " gene sets.")
	setNet = set.graph.attribute(setNet, "analysis", 
			if (use.ranks) "ranked" else "unranked")
	setNet = pAdjustAndCorrect(setNet)
	notSignificant = which(V(setNet)$adjustedPValue > fdrCutoff)
	message(length(notSignificant), " gene sets removed after FDR correction.")
	setNet - V(setNet)[notSignificant]
}

as.RankedTestSet <- function(geneIDs, setCollection) {
	ranks = 1:length(geneIDs)
	names(ranks) <- geneIDs
	geneIDs = geneIDs %i% setCollection$referenceSet
	testSet = rank(ranks[geneIDs])
	class(testSet) <- "RankedTestSet"
	testSet
}

as.UnrankedTestSet <- function(geneIDs, setCollection) {
	testSet = as.character(geneIDs) %i% setCollection$referenceSet
	class(testSet) <- "UnrankedTestSet"
	testSet
}

getPrimarySetPValues <- function(testSet, setCollection) {
	pValues = rep(1, length(setCollection$sets))
	names(pValues) <- names(setCollection$sets)
	smallSetIndices = !setCollection$bigSets
	pValues[smallSetIndices] = 
			unlist(mclapply(setCollection$sets[smallSetIndices], getSetPValue,
							testSet, setCollection))
	pValues	
}

getSetPValue <- function(geneSet, testSet, setCollection)
	UseMethod("getSetPValue", testSet)

#' @S3method getSetPValue UnrankedTestSet
getSetPValue.UnrankedTestSet <- function(geneSet, testSet, setCollection) {
	m = length(geneSet)
	i = length(geneSet %i% testSet)
	s = length(testSet)
	fisherPValue(setCollection, m, i, s)
}

#' @S3method getSetPValue RankedTestSet
getSetPValue.RankedTestSet <- function(geneSet, testSet, setCollection) {
	ranks = sort(testSet[geneSet])
	ranks = ranks[!is.na(ranks)]
	if (length(ranks) == 0) return(1)
	m = length(ranks)
    setCollection$g = length(testSet)
	min(sapply(1:length(ranks), function(i) 
								fisherPValue(setCollection, m, i, ranks[i])))
}

fisherPValue <- function(setCollection, m, i, s) {
	g = setCollection$g
	if (i == 0) {
		return(1)
	}
	minimalI = setCollection$iMatrix[s,m]
	if (is.na(minimalI) || i < minimalI) {
		return(1)
	} else {
		return(fisher.test(rbind(c(i, m-i),c(s-i, g-(m+s-i))), 
						alternative="greater")$p.value)
	}
}

buildEdgeTable  <- function(testSet, setCollection, setPValues, setPCutoff) {
	significantSetIDs = setPValues <= setPCutoff
	message(Sys.time(), " - ", length(which(significantSetIDs)), " significant sets")
	intersectionTable = setCollection$intersections
	bothSignificant = significantSetIDs[intersectionTable$setA] & 
			significantSetIDs[intersectionTable$setB]
	intersectionTable = intersectionTable[bothSignificant,]
	if (nrow(intersectionTable) == 0) {
		return(data.frame(source=NA, sink=NA, type=NA))
	}
	n = length(which(significantSetIDs)) * 
			(length(which(significantSetIDs))-1)/2
	nNodes = if (nrow(intersectionTable) < getOption("mc.cores")) 1 else getOption("mc.cores")
	message(Sys.time(), " - ", nrow(intersectionTable),
			" intersections to test out of ", n, " possible intersections",
			" using ", nNodes, " cores.")
	cluster = if (.Platform$OS.type == "windows") makePSOCKcluster(nNodes) else 
            makeForkCluster(nNodes)
	pairStats = parRapply(cluster, intersectionTable, getSetPairStatistics, 
			testSet, setCollection)
	stopCluster(cluster)
	message(Sys.time(), " - intersections evaluated.")
	edgeTable = rbindlist(pairStats)
	message(Sys.time(), " - edge table constructed.")
	edgeTable$discardSource = FALSE
	edgeTable$discardSink = FALSE
	discardSourceIndices = edgeTable$source_pDiff > setPCutoff & 
			edgeTable$type != "subset"
	if (any(discardSourceIndices)) {
		edgeTable[discardSourceIndices,]$discardSource = TRUE
	}
	discardSinkIndices = edgeTable$sink_pDiff > setPCutoff
	if (any(discardSinkIndices)) {
		edgeTable[discardSinkIndices,]$discardSink = TRUE
	}
	if (any(edgeTable$discardSource & edgeTable$discardSink)) {
		edgeTable[edgeTable$discardSource & 
						edgeTable$discardSink,]$type = "intersection"
		edgeTable[edgeTable$type == "intersection",]$discardSource = FALSE
		edgeTable[edgeTable$type == "intersection",]$discardSink = FALSE
	}
	edgeTable = edgeTable[edgeTable$significantJaccard > 0,]
	edgeTable
}

getSetPairStatistics <- function(row, testSet, setCollection)
	UseMethod("getSetPairStatistics", testSet)

getSetPairStatistics_base <- function(row, testSet, setCollection) {
	setIDA = row[1]
	setIDB = row[2]
	setA = setCollection$sets[[setIDA]]
	setB = setCollection$sets[[setIDB]]
	intersection = setA %i% setB
	a = list(id = setIDA, size = length(setA))
	b = list(id = setIDB, size = length(setB))
	unionSet = setA %u% setB
	diffA = setA %d% setB
	diffB = setB %d% setA
	type = "overlap"
	pIntersection = getSetPValue(intersection, testSet, setCollection)
	a$pDiff = getSetPValue(diffA, testSet, setCollection)
	a$diffSize = length(diffA)
	b$pDiff = getSetPValue(diffB, testSet, setCollection)
	b$diffSize = length(diffB)
	if (length(diffA) == 0 || length(diffB) == 0) {
		type = "subset"
		if (length(diffA) == 0) {
			source = a
			sink = b
		} else {
			source = b
			sink = a
		}
		source$pDiff = 1.0
	} else if (a$pDiff > b$pDiff) {
		source = a
		sink = b
	} else {
		source = b
		sink = a
	}
	jaccard = length(intersection)/length(unionSet)
	data.frame(source=source$id, sink=sink$id, type=type, 
			source_pDiff = source$pDiff, source_ppDiff = -log10(source$pDiff), 
			source_diffSize = source$diffSize, 
			sink_pDiff = sink$pDiff, sink_ppDiff = -log10(sink$pDiff), 
			sink_diffSize = sink$diffSize, 
			deltaP = -log10(sink$pDiff) + log10(source$pDiff),
			intersectionSize = length(intersection), 
			pIntersection = pIntersection, 
			ppIntersection = -log10(pIntersection), jaccard=jaccard, 
			stringsAsFactors=FALSE)	
}

#' @S3method getSetPairStatistics RankedTestSet
getSetPairStatistics.RankedTestSet <- function(row, testSet, setCollection) {
	stats = getSetPairStatistics_base(row, testSet, setCollection)
	stats$significantJaccard = stats$jaccard
	stats
}

#' @S3method getSetPairStatistics UnrankedTestSet
getSetPairStatistics.UnrankedTestSet <- function(row, testSet, setCollection) {
	stats = getSetPairStatistics_base(row, testSet, setCollection)
	source = setCollection$sets[[stats$source]]
	sink = setCollection$sets[[stats$sink]]
	intersection = source %i% sink
	unionSet = source %u% sink
	stats$intersectionSignificant = length(intersection %i% testSet)
	stats$source_diffSignificant = length(source %i% testSet)
	stats$sink_diffSignificant = length(sink %i% testSet)
	stats$significantJaccard = length(intersection %i% testSet) /
			length(unionSet %i% testSet)
	stats
}

getNodesToDelete <- function(edgeTable) {
	if (is.na(edgeTable[1,]$sink)) {return(c())}
	superSetsToDelete = edgeTable[edgeTable$type == "subset" & 
					edgeTable$discardSink,]$sink
	killTable = edgeTable[edgeTable$type == "overlap",]
	while (TRUE) {
		killTable = killTable[killTable$discardSource,]
		allKillers = unique(killTable$sink)
		allKilled = unique(killTable$source)
		topKillers = allKillers %d% allKilled
		intermediateKillers = 
				unique(killTable[killTable$sink %in% topKillers,]$source) %i%
				allKillers
		invalidKills = killTable$sink %in% intermediateKillers
		if (any(invalidKills)) {
			killTable[invalidKills,]$discardSource = FALSE
		} else {
			break;
		}
	}
	return(unique(superSetsToDelete %u% killTable$source))
}

buildSetNet <- function(edgeTable, setCollection, testSet, pValues, setPCutoff)
{
	vertexTable = data.frame(name=sapply(setCollection$sets, attr, "ID"), 
			description = sapply(setCollection$sets, attr, "name"), 
			database=sapply(setCollection$sets, attr, "db"),  pValue = pValues,
			pp = -log10(pValues), size = sapply(setCollection$sets, length),
			stringsAsFactors=FALSE)
	vertexTable$nSignificant = sapply(setCollection$sets, 
			function(x) length(x %i% testSet))
	vertexTable = vertexTable[vertexTable$pValue <= setPCutoff,]
    edgeTable = edgeTable[edgeTable$pIntersection <= setPCutoff,]
	setNet = if (is.na(edgeTable[1,]$source)) {
				add.vertices(graph.empty(), nrow(vertexTable),
						attr=as.list(vertexTable))
			} else {
				graph.data.frame(edgeTable, directed=TRUE, 
						vertices=vertexTable)
			}
	setNet
}

pAdjustAndCorrect <- function(setNet) {
	if (vcount(setNet) == 0) {
		return(setNet)
	}
	edgeTable = get.data.frame(setNet, what="edges")
	nodeTable = get.data.frame(setNet, what="vertices")
	nodeTable$correctedPValue = nodeTable$pValue
	if (ecount(setNet) > 0) {
		correctedPValues = getCorrectedPValues(edgeTable)
		nodeTable[names(correctedPValues),]$correctedPValue = correctedPValues
	}
	nodeTable$pp = -log10(nodeTable$correctedPValue)
	setNet = graph.data.frame(edgeTable, directed=TRUE, vertices=nodeTable)
	setNet = getAdjustedPValues(setNet)
	return(setNet)
}

getMaxP <- function(edgeTable, type, idAttribute, pAttribute) {
	subTable = edgeTable[edgeTable$type == type,]
	maxPList = by(subTable, as.factor(subTable[[idAttribute]]),
			function(x) max(x[[pAttribute]]), simplify=FALSE)
	unlist(maxPList)
}

getCorrectedPValues <- function(edgeTable) {
	maxPValues = list(
			subset = getMaxP(edgeTable, "subset", "to", "sink_pDiff"),
			overlap = getMaxP(edgeTable, "overlap", "from", "source_pDiff"),
			intersectFrom = getMaxP(edgeTable, "intersection", "from", 
					"pIntersection"),
			intersectTo = getMaxP(edgeTable, "intersection", "to", 
					"pIntersection"))
	setIDs = unique(unlist(lapply(maxPValues, names), use.names=FALSE))
	maxPTable = do.call(cbind, lapply(maxPValues, function(x) x[setIDs]))
	rownames(maxPTable) <- setIDs
	return(apply(maxPTable, 1, max, na.rm=TRUE))
}

getAdjustedPValues <- function(setNet) {
	componentData = clusters(setNet)
	componentMembers = lapply(1:componentData$no, 
			function(x) which(componentData$membership == x))
	componentPValues = p.adjust(sapply(componentMembers, 
					function(x) min(V(setNet)[x]$correctedPValue)))
	for (i in 1:componentData$no) {
		V(setNet)[componentMembers[[i]]]$adjustedPValue = componentPValues[i]
	}
	setNet
}

sinkDelete <- function(setNet, setCollection, setPCutoff, testSet) {
	edgeTable = get.data.frame(setNet)
	edgeTable = edgeTable[edgeTable$type=="overlap",]
	pValues = unlist(by(edgeTable, edgeTable$from, sinkDeletePValue, 
					setCollection, testSet))
	toDelete = which(pValues > setPCutoff)
	setNet - V(setNet)[toDelete]
}

sinkDeletePValue <- function(subTable, setCollection, testSet) {
	source = setCollection$sets[[unique(subTable$from)]]
	discardMatrix = do.call(cbind, lapply(subTable$to, function(s) {
						sink = setCollection$sets[[s]]
						source %in% sink
					}))
	discard = apply(discardMatrix, 1, any)
	return(if (all(discard)) 1 else getSetPValue(source[!discard], testSet, 
								setCollection))
}

calculateSetRank <- function(setNet) {
	subsetEdges = which(E(setNet)$type == "subset")
	subNet = setNet - E(setNet)[subsetEdges]
	setRank = page.rank(subNet, options=list(maxiter=50000))
	setNet = set.vertex.attribute(setNet, "setRank", 
			index=names(setRank$vector), value=setRank$vector)
	setNet = set.vertex.attribute(setNet, "pSetRank", 
			index=names(setRank$vector), value=setRankPValue(setNet))
	setNet
}

setRankPValue <- function(net, minN=10000) {
	nNodes = vcount(net)
	if (nNodes == 0) {
		return(c())
	}
	nEdges = ecount(net)
	nReplicates = 1 #round(minN/nNodes)
	simSetRank = unlist(sapply(1:nReplicates, function(n) {
						g = erdos.renyi.game(nNodes, nEdges, type="gnm", 
								directed=TRUE)
						log(page.rank(g, options=list(maxiter=50000))$vector)
					}))
	p.adjust(pnorm(log(V(net)$setRank), mean(simSetRank), sd(simSetRank), 
					lower.tail=FALSE))
}
