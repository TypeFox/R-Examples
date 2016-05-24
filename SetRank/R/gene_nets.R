# TODO: Add comment
# 
# Author: cesim
###############################################################################


setGlobalAttributes <- function(network, topTableList, fields) {
	rankField = "rank"
	geneIDField = fields["geneID"]
	symbolField = fields["SYMBOL"]
	logFCField = fields["logFC"]
	pField = fields["p"]
	minP = 1e-100
	maxStrength = -log10(minP)
	minBetweenness = -10
	nodeTable = get.data.frame(network, what="vertices")
	for (name in names(topTableList)) {
		message(name)
		topTable = topTableList[[name]]
		commonGenes = topTable[[geneIDField]] %i% rownames(nodeTable)
		pAttribute = sprintf("%s:p", name)
		logFCAttribute = sprintf("%s:logFC", name)
		ppAttribute = sprintf("%s:pp", name)
		rankAttribute = sprintf("%s:rank", name)
		nodeTable[[pAttribute]] = 1
		nodeTable[commonGenes,][[pAttribute]] = 
				topTable[commonGenes,][[pField]]
		nodeTable[[rankAttribute]] = 0
		nodeTable[commonGenes,][[rankAttribute]] = 
				topTable[commonGenes,][[rankField]]
		if (any(nodeTable[[pAttribute]] == 0)) {
			nodeTable[nodeTable[[pAttribute]] < minP,][[pAttribute]] = minP
		}
		nodeTable[[ppAttribute]] = -log10(nodeTable[[pAttribute]])
		nodeTable[[logFCAttribute]] = NA
		nodeTable[commonGenes,][[logFCAttribute]] = 
				topTable[commonGenes,][[logFCField]]
		network = set.vertex.attribute(network, pAttribute, 
				value=nodeTable[[pAttribute]])
		network = set.vertex.attribute(network, logFCAttribute,
				value=nodeTable[[logFCAttribute]])
		network = set.vertex.attribute(network, ppAttribute,
				value=nodeTable[[ppAttribute]])
		network = set.vertex.attribute(network, rankAttribute,
				value=nodeTable[[rankAttribute]])
		edgeTable = get.data.frame(network, what="edges")
		costAttribute = sprintf("%s:cost", name)
		fromIDs = edgeTable$from
		toIDs = edgeTable$to
		strength = (nodeTable[fromIDs,][[ppAttribute]] + 
					nodeTable[toIDs,][[ppAttribute]])/2
		strength[strength > maxStrength] = maxStrength
		cost = (maxStrength-strength)+1
		network = set.edge.attribute(network, costAttribute, value=cost) 
		betweennessAttribute = sprintf("%s:globalBetweenness", name)
		nodeBetweenness = log10(betweenness(network, directed=FALSE, 
						weights=cost, normalized=TRUE))
		nodeBetweenness[nodeBetweenness < minBetweenness] = minBetweenness
		network = set.vertex.attribute(network, betweennessAttribute, 
				value=nodeBetweenness)
	}
	network
}

getMissingNodes <- function(net, missingGenes, topTables, fields) {
	rankField = "rank"
	geneIDField = fields["geneID"]
	symbolField = fields["symbol"]
	logFCField = fields["logFC"]
	pField = fields["p"]
	mergedTable = do.call(rbind, lapply(topTables, 
					function(x) x[,c(geneIDField,symbolField)]))
	uniqueTable = mergedTable[!duplicated(mergedTable[[geneIDField]]),]             
	nodeTable = data.frame(name=uniqueTable[[geneIDField]], 
			symbol=uniqueTable[[symbolField]], stringsAsFactors=FALSE)
	rownames(nodeTable) <- nodeTable$name
	nodeTable = nodeTable[missingGenes %i% rownames(nodeTable),]
	if (nrow(nodeTable) == 0) {
		return(net)
	}
	for (name in names(topTables)) {
		topTable = topTables[[name]]
		pAttribute = sprintf("%s:p", name)
		logFCAttribute = sprintf("%s:logFC", name)
		ppAttribute = sprintf("%s:pp", name)
		betweennessAttribute = sprintf("%s:globalBetweenness", name)
		rankAttribute = sprintf("%s:rank", name)
		nodeTable[[pAttribute]] = 1
		nodeTable[[logFCAttribute]] = NA
		commonGenes = rownames(nodeTable) %i% rownames(topTable)
		nodeTable[commonGenes,][[pAttribute]] = 
				topTable[commonGenes,][[pField]]
		nodeTable[commonGenes,][[logFCAttribute]] = 
				topTable[commonGenes,][[logFCField]]
		nodeTable[[betweennessAttribute]] = NA
		nodeTable[[rankAttribute]] = 0.0
	}
	net = add.vertices(net, nv=nrow(nodeTable), attr=as.list(nodeTable))
	net
}

writeGeneSetNetwork = function(geneSet, interactome, topTables, outPath, fields) {
	minBetweenness = -10
	nodeTable = get.data.frame(interactome, what="vertices")
	commonGenes = geneSet %i% rownames(nodeTable)
	missingGenes = geneSet %d% commonGenes
	setInteractome = induced.subgraph(interactome, 
			which(rownames(nodeTable) %in% commonGenes))
	setNodeTable = get.data.frame(setInteractome, what="vertices")
	for (name in names(topTables)) {
		globalAttribute = sprintf("%s:globalBetweenness", name)
		localAttribute = sprintf("%s:localBetweenness", name)
		ratioAttribute = sprintf("%s:betweennessRatio", name)
		costAttribute = sprintf("%s:cost", name)
		costs = get.edge.attribute(setInteractome, costAttribute)
		localBetweenness = log10(betweenness(setInteractome, directed=FALSE, 
						weights=costs, normalized=TRUE))
		localBetweenness[localBetweenness < minBetweenness] = minBetweenness
		betweennessRatio = localBetweenness - 
				get.vertex.attribute(setInteractome, globalAttribute)
		setInteractome = set.vertex.attribute(setInteractome, localAttribute, 
				value=localBetweenness)
		setInteractome = set.vertex.attribute(setInteractome, ratioAttribute,
				value=betweennessRatio)
	}
	setInteractome = getMissingNodes(setInteractome, missingGenes, topTables,
			fields)
	outFileName = sprintf("%s/%s.%s.net.xml", outPath, 
			make.names(attr(geneSet, "name")), attr(geneSet, "db"))
	write.graph(setInteractome, outFileName, format="graphml")
	fixGraphML(outFileName)
	setInteractome
}

createGeneNetVizFile <- function(names, outDir) {
	templateFile = system.file("extdata", "gene_net_style_template.xml", 
			package = "SetRank")
	templateLines = readLines(templateFile)
	header = templateLines[1:2]
	footer = templateLines[length(templateLines)]
	visStyleTemplate = templateLines[3:(length(templateLines)-1)]
	visStyles = sapply(names, function(n) gsub("setrank", n, visStyleTemplate))
	writeLines(c(header, visStyles, footer), 
			sprintf("%s/gene_net_styles.viz.xml", outDir))
}

#' Fix the igraph graphML export
#' 
#' This function only exists to fix some issues with the graphML output of the 
#' igraph module. 
#'
#' @param fileName the name of the file to fix.
#' @author Cedric Simillion
fixGraphML <- function(fileName) {
	netName = sub("\\.net\\.xml$", "", basename(fileName))
	lines = readLines(fileName)
	lines = sub("graph id=\"G\"",sprintf("graph id =\"%s\"", netName), 
			lines)
	lines = gsub("\\>inf\\<", ">1<", lines)
	writeLines(lines, fileName)
}

#' Create gene interaction networks for all significant gene sets.
#' 
#' Creates for every gene set present in one or more gene set networks
#' a gene interaction network. This network shows all known or predicted 
#' protein-protein interactions between all genes in the gene set. Each network
#' is written out to a file in GraphML format, with the extension .net.xml. 
#' Additionally, a file called gene_net_styles.viz.xml is created as well. This
#' file contains for each gene set network a Cytoscape visualisation style
#' which can be used to overlay the original expression data on top of the
#' gene set-specific interaction network for. See the vignette for more 
#' details.
#' 
#' @param topTables A named list object containing one or more data frames. The
#' names should be the same as those of \code{networks} argument. Each data 
#' frame should contain the expression analysis data used as input for the 
#' corresponding network in the \code{networks} arguments. They should at 
#' least contain four columns listing the NCBI Entrez Gene identifier, the gene
#' symbol, the observed log-fold change and the (adjusted) p-value of each 
#' gene. See the description of the \code{fields} argument for details.
#' @param networks A named list object containing one or more gene set networks
#' as returned by the \code{\link{setRankAnalysis}} function, based on the data
#' in the \code{topTables} argument.
#' @param collection The gene set collection object used to create the gene set
#' networks in the \code{networks} argument.
#' @param string An igraph object containing a species-specific protein-protein
#' interaction network from which to retrieve. You can use the SetRankTools
#' set of scripts to generate this object for your species of interest.
#' @param outDir The directory where to write the output files. If this 
#' directory doesn't exist, it will be created.
#' @param geneSetIDs The list of gene set identifiers for which to create gene
#' interaction networks. When omitted, the union of all gene sets found in the 
#' \code{networks} argument will be used.
#' @param fields A named vector of strings specifying the column names of the
#' data frames in the \code{topTables} argument. 
#' @author Cedric Simillion
#' @export
exportGeneNets <- function(topTables, networks, collection, string, outDir, 
		geneSetIDs = NULL, fields = c(geneID = "ENTREZID", symbol = "SYMBOL", 
				logFC = "logFC", p = "adj.P.Val")) {
	rankField = "rank"
	geneIDField = fields["geneID"]
	symbolField = fields["SYMBOL"]
	logFCField = fields["logFC"]
	pField = fields["p"]
	topTables = lapply(topTables, function(x) {
				x = x[!is.na(x[[geneIDField]]),]
				rownames(x) <- x[[geneIDField]]
				x$rank = 1 - 0:(nrow(x)-1)/nrow(x)
				x
			})
	stopifnot(names(topTables) == names(networks))
	if (!file.exists(outDir)) {
		dir.create(outDir, recursive=TRUE)
	}
	createGeneNetVizFile(names(networks), outDir)
	nodeTables = lapply(networks, get.data.frame, what="vertices")
	if (is.null(geneSetIDs)) {
		geneSetIDs = unique(unlist(lapply(nodeTables, rownames)))
	}
	message("Calculating global attributes...")
	string = setGlobalAttributes(string, topTables, fields)
	message("Writing pathway networks...")
	dir.create(outDir, showWarnings=FALSE)
	setNets = list()
	for (setID in geneSetIDs) {
		geneSet = collection$sets[[setID]]
		message(setID, ": ", attr(geneSet, "name"))
		setNets[[setID]] = writeGeneSetNetwork(geneSet, string, topTables, 
				outDir, fields)
	}
}
