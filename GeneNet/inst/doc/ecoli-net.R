# /*
# This is an R script containing R markdown comments.  It can be run as is in R.
# To generate a document containing the formatted R code, R output and markdown 
# click the "Compile Notebook" button in R Studio, or run the command
# rmarkdown::render() - see http://rmarkdown.rstudio.com/r_notebook_format.html
# */


#' ---
#' title: "Escherichia Coli Network"
#' output: pdf_document
#' author: ""
#' date: Example for GeneNet 1.2.13 (August 2015) or later
#' ---

#' This note reproduces the “Escherichia coli” network example from J. Schäfer and 
#' K. Strimmer. 2005. *A shrinkage approach to large-scale covariance 
#' estimation and implications for functional genomics.* 
#' Statist. Appl. Genet. Mol. Biol. **4**: 32 
#' (http://dx.doi.org/10.2202/1544-6115.1175)


#'
#' # Load GeneNet package

library("GeneNet")

#' E. Coli data set (9 time points for 102 genes):
data(ecoli)
dim(ecoli)


#'
#' # Estimation of partial correlations

#' Estimate matrix of partial correlation using a shrinkage estimator:
pc = ggm.estimate.pcor(ecoli)
dim(pc)

#' Assign p-values, q-values and empirical posterior probabilities to all
#' 5151 potential edges in the network:
ecoli.edges = network.test.edges(pc, direct=TRUE, fdr=TRUE)
dim(ecoli.edges)

#' The table lists all edges in the order strength of partial correlations:
ecoli.edges[1:5,]

#'
#' # Decide which edges to include in the network

#' To obtain a graph you need to select top ranking edges according to 
#' a suitable criterion.  Here are some suggestions:
#'
#' 1. Use local fdr cutoff 0.2, i.e. include all edges with posterior 
#' probability of at least 0.8.
ecoli.net = extract.network(ecoli.edges)
dim(ecoli.net)

#' 2. Use local fdr cutoff 0.1, i.e. i.e. include all edges with posterior 
#' probability of at least 0.9.

ecoli.net = extract.network(ecoli.edges, cutoff.ggm=0.9, cutoff.dir=0.9)
dim(ecoli.net)

#' 3. Include a fixed number of edges, say the 70 strongest edges
ecoli.net = extract.network(ecoli.edges, method.ggm="number", cutoff.ggm=70)
dim(ecoli.net)

#' 
#' Plot network

#' For plotting we use the graph and Rgraphviz packages from Bioconductor.
library("Rgraphviz") 

#' Create graph object from the list of edges:
node.labels = colnames(ecoli)
gr = network.make.graph(ecoli.net, node.labels, drop.singles=TRUE)
table(  edge.info(gr)$dir )
sort( node.degree(gr), decreasing=TRUE)


#' Set node and edge attributes for more beautiful graph plotting:
globalAttrs = list()
globalAttrs$edge = list(color = "black", lty = "solid", lwd = 1, arrowsize=1)
globalAttrs$node = list(fillcolor = "lightblue", shape = "ellipse", fixedsize = FALSE)
 
nodeAttrs = list()
nodeAttrs$fillcolor = c('sucA' = "yellow")

edi = edge.info(gr)
edgeAttrs = list()
edgeAttrs$dir = edi$dir # set edge directions 
edgeAttrs$lty = ifelse(edi$weight < 0, "dotted", "solid") # negative correlation -> dotted
edgeAttrs$color = ifelse(edi$dir == "none", "black", "red")
edgeAttrs$label = round(edi$weight, 2) # use partial correlation as edge labels

#+ fig.width=8, fig.height=7
plot(gr, attrs = globalAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs, "fdp")

