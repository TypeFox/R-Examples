# Project: SetRankVisualization
# 
# Author: cesim
###############################################################################
#' Export multiple SetRank networks and accompanying tables.
#' 
#' Given a list of SetRank networks, writes out the following files for each
#' network:
#' \enumerate{
#' \item SetRank network in GML format called \code{<n>.gml} where \code{<n>}
#' is the network name.
#' \item A TAB-delimited file listing the signficant pathways in the network,
#' called \code{<n>_pathways.txt} with \code{<n>} the network name.
#' \item A TAB-delimited file listing which significant genes belong to which
#' pathway, called  \code{<n>_membership.txt} with \code{<n>} again the network
#' name.
#' }
#' The network names will be taken from the names of the input list.
#' Additionally, two global files will be created as well:
#' \enumerate{
#' \item A Cytoscape VizMap visualisation file, called\code{setrank.xml}.
#' \item A TAB-delimited file listing which pathways are found in which
#' networks, called "pathways.txt".
#' }
#' 
#' @param networkList A named list of SetRank networks.
#' @param selectedGenesList A named list with the same names as the 
#' \code{networkList} argument. Each list should be a vector with the set of 
#' significant genes as Entrez Gene IDs used to construct the SetRank network
#' with the same name.
#' @param collection The set collection used for the SetRank analysis.
#' @param IDConverter Optional. By default, Entrez Gene IDs will be displayed
#' in the output tables. This argument can be used to convert these into more 
#' human-friendly gene symbols. When supplied, should be a function that takes 
#' a vector of Entrez Gene IDs as single argument and returns the values of the
#' corresponding gene symbols or whatever identifier you wish to have displayed
#' in the output tables.
#' @param outputPath The name of the directory where the results should be 
#' written. If the last element of the path doesn't exist, a directory will be
#' created.
#' @return None. Files are written out as a side effect. 
#' @export
exportMultipleResults <- function(networkList, selectedGenesList, collection, 
		IDConverter = NULL, outputPath="./") {
	if (!file.exists(outputPath)) {
		dir.create(outputPath)
	}
	for (n in names(networkList)) {
		write.table(getNodeTable(networkList[[n]]), 
				sprintf("%s/%s_pathways.txt", outputPath, n), sep="\t", 
				row.names=FALSE, quote=FALSE)
		writeMembership(sprintf("%s/%s_membership.txt", outputPath, n),
				selectedGenesList[[n]], collection, networkList[[n]], 
				IDConverter)
	}
	cytoscapeExport(networkList, outputPath)
	pathways = createPathwayTable(networkList, collection)
	write.table(pathways, sprintf("%s/pathways.txt", outputPath), sep="\t",
			row.names=FALSE, quote=FALSE)
}

#' Export a SetRank network and accompanying tables.
#' 
#' Given a single SetRank analysis result, writes out the following files:
#' \enumerate{
#' \item SetRank network in GML format called \code{<n>.gml} where \code{<n>}
#' is the specified network name.
#' \item A Cytoscape VizMap visualisation file, called \code{setrank.xml}.
#' \item A TAB-delimited file listing the signficant pathways in the network,
#' called \code{<n>_pathways.txt} with \code{<n>} again the network name.
#' \item A TAB-delimited file listing which significant genes belong to which
#' pathway, called  \code{<n>_membership.txt} with \code{<n>} again the network
#' name.
#' }
#' 
#' @param network A SetRank network.
#' @param selectedGenes A vector with the set of significant genes as Entrez
#' Gene IDs. This should be the same set passed to the 
#' \code{\link[SetRank]{setRankAnalysis}} function.
#' @param collection The set collection used for the SetRank analysis.
#' @param networkName A name used to name the different output files.
#' @param IDConverter Optional. By default, Entrez Gene IDs will be displayed
#' in the output tables. This argument can be used to convert these into more 
#' human-friendly gene symbols. When supplied, should be a function that takes 
#' a vector of Entrez Gene IDs as single argument and returns the values of the
#' corresponding gene symbols or whatever identifier you wish to have displayed
#' in the output tables.
#' @param outputPath The name of the directory where the results should be 
#' written. If the last element of the path doesn't exist, a directory will be
#' created.
#' @return None. Files are written out as a side effect. 
#' genes = sprintf("gene_%03d", 1:100)
#' geneSets = lapply(1:9, function(i) sample(genes[((i-1)*10):((i+1)*10)], 10))
#' annotationTable = data.frame(termID=sprintf("set_%02d", rep(1:9, each=10)), 
#'         geneID=unlist(geneSets),
#'         termName = sprintf("dummy gene set %d", rep(1:9, each=10)),
#'         dbName = "dummyyDB",
#'         description = "A dummy gene set DB for testing purposes")
#' collection = buildSetCollection(annotationTable, referenceSet=genes)
#' network = setRankAnalysis(genes, collection, TRUE)
#' exportSingleResult(network, genes, collection, "example", function(x) x, "example_dir")
#' @author Cedric Simillion
#' @export
exportSingleResult <- function(network, selectedGenes, collection, networkName,
		IDConverter = NULL, outputPath="./") {
	if (!file.exists(outputPath)) {
		dir.create(outputPath)
	}
	write.table(getNodeTable(network), sprintf("%s/%s_pathways.txt", outputPath, 
					networkName), sep="\t", row.names=FALSE, quote=FALSE)
	writeMembership(sprintf("%s/%s_membership.txt", outputPath, networkName),
			selectedGenes, collection, network, IDConverter)
	write.graph(network, sprintf("%s/%s.net.xml", outputPath, networkName),
			format="graphML")
	createCytoscapeVizMap(network=network, outputFile=sprintf("%s/setrank.xml",
					outputPath))
}

writeMembership <- function(outputFile, selectedGenes, collection, network, 
		IDConverter=NULL) {
	membershipBool = membershipTable(selectedGenes, collection, network)
	if (is.null(membershipBool)) {
		return(NULL)
	}
	membership = matrix(".", nrow=nrow(membershipBool), 
			ncol=ncol(membershipBool), dimnames=dimnames(membershipBool))
	membership[membershipBool == TRUE] = "X"
	if (!is.null(IDConverter)) {
		rownames(membership) = IDConverter(rownames(membership))
	}
	membership = orderTable(membership)
	prettyTable = cbind(rownames(membership),membership)
	colnames(prettyTable)[1] = "gene"
	write.table(prettyTable, outputFile, sep="\t", row.names=FALSE, 
			quote=FALSE)
}

orderTable <- function(membershipTable){
    membershipTableMod <- ifelse(membershipTable=="X",1,0)
	if (nrow(membershipTableMod) > 1)  {
		clustering = hclust(dist(membershipTableMod))
		geneOrder = clustering$labels[clustering$order]
	} else {
		geneOrder = 1
	}
	if (ncol(membershipTableMod) > 1) {
		clustering = hclust(dist(t(membershipTableMod)))
		catOrder = clustering$labels[clustering$order]
	} else {
		catOrder = 1
	}
    return(membershipTable[geneOrder,catOrder])
}

getNodeTable <- function(network) {
	table = get.data.frame(network, what="vertices")
	table$pValue = NULL
	table$nSignificant = NULL
	table$pp = NULL
	table[order(table$pSetRank, table$adjustedPValue, table$correctedPValue),]
}

#' Creates a table of all significant pathways in different conditions.
#' 
#' @param networkList A list of SetRank networks created using the same
#' set collection.
#' @param setCollection The set collection used to perform the SetRank 
#' analysis.
#' @return A data frame with column names being the names of the networks in
#' the \code{networkList} argument and the rownames gene set IDs. The cells
#' contain the adjusted p-values of each gene set in each network. When a gene 
#' set is not present in a network, the p-value will be set to 1. Two
#' additional columns called "\code{description}" and "\code{score}" are added.
#' The former is simply the description of the gene set. The latter is a score
#' which attempts to reflect the importance of a gene set across the difference
#' networks. The higher the score, the more important the network. This score
#' is a combination of the number of networks where the gene set is observed
#' and the geometric mean of the p-values of that set in these networks.
#' @author Cedric Simillion
#' @export 
createPathwayTable <- function(networkList, setCollection) {
	nodeTables = lapply(networkList, getNodeTable)
	rankScores = lapply(nodeTables, function(t) {
				scores = rev(seq_along(t$name))/nrow(t)
				names(scores) <- t$name
				scores
			})
	allPathIDs = unique(do.call(rbind, nodeTables)$name)
	pathwayTable = as.data.frame(do.call(cbind, lapply(rankScores, 
							function(x) x[allPathIDs])))
	dimnames(pathwayTable) <- list(allPathIDs, names(networkList))
	pathwayTable[is.na(pathwayTable)] = 0
	pathwayTable$score = apply(pathwayTable, 1, 
			function(row) length(which(row > 0)) + mean(row))
	pathwayTable$name = rownames(pathwayTable)
	pathwayTable$description = sapply(rownames(pathwayTable), 
			function(x) attr(setCollection$sets[[x]], "name"))
	tableOrder = order(pathwayTable$score, decreasing=TRUE)
	pathwayTable[tableOrder,]
}

#' Create a gene set membership table.
#' 
#' Creates a table showing which significant genes belong to which significant
#' gene sets. This table allows to investigate the results of a SetRank 
#' analysis in more detail.
#' 
#' @param selectedGenes The set of significant genes used during the SetRank
#' analysis.
#' @param collection The setCollection used during the SetRank analysis.
#' @param network A network generated by the 
#' \code{\link[SetRank]{setRankAnalysis}} function.
#' @references A matrix of boolean values. The rownames are the geneIDs, the
#' column names are the gene set IDs. A \code{TRUE} value in a cell indicates
#' that the gene of the row is present in the gene set of the column; 
#' \code{FALSE} means otherwise.
#' @author Cedric Simillion
#' @export 
membershipTable <- function(selectedGenes, collection, network) {
	pathwayIDs = V(network)$name
	pathways = lapply(collection$sets[pathwayIDs], 
			function(p) p %i% selectedGenes)
	genesInPathways = unique(unlist(pathways, use.names=FALSE))
	table = do.call(rbind, lapply(genesInPathways, function(g) 
						sapply(pathways, function(p) g %in% p)))
	if (!is.null(table)) dimnames(table) <- list(genesInPathways, pathwayIDs)
	table
}

getPValueTable <- function(topTables, pColumn="adj.P.Val", idColumn="ENTREZID") {
	allGenes = topTables[[1]][[idColumn]]
	for (n in names(topTables)) {
		rownames(topTables[[n]]) <- topTables[[n]][[idColumn]]
	}
	pValueTable = do.call(cbind, lapply(topTables, 
					function(x) x[allGenes,]$adj.P.Val))
	rownames(pValueTable) <- allGenes
	return(pValueTable)
}

getLogFCTable <- function(topTables, pValueTable, idColumn="ENTREZID", 
		valColumn="logFC", pValueCutoff = 0.05) {
	for (tableName in names(topTables)) {
		rownames(topTables[[tableName]]) <- topTables[[tableName]][[idColumn]]
	}
	movingIndices = apply(pValueTable, 1, function(row) 
				any(row <= pValueCutoff))
	movingGenes = rownames(pValueTable[movingIndices,])
	logFCTable = do.call(cbind, lapply(topTables, 
					function(x) x[movingGenes,]$logFC))
	dimnames(logFCTable) <- list(movingGenes, names(topTables))
	return(logFCTable)
}

getConditionOrder <- function(logFCTable) {
	hclust(as.dist(1-cor(logFCTable)))$order
}

boolJaccard <- function(x,y) {
	nominator = length(which(x|y))
	if (nominator == 0) {
		return(0.0)
	} else {
		return(length(which(x&y))/nominator)
	}
}

jaccardDistanceMatrix <- function(membershipTable) {
	matrix = apply(membershipTable, 2, 
			function(x) apply(membershipTable, 2, boolJaccard, x))
	as.dist(matrix)
}
