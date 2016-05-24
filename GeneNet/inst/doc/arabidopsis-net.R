# /*
# This is an R script containing R markdown comments.  It can be run as is in R.
# To generate a document containing the formatted R code, R output and markdown 
# click the "Compile Notebook" button in R Studio, or run the command
# rmarkdown::render() - see http://rmarkdown.rstudio.com/r_notebook_format.html
# */

#' ---
#' title: "Arabidopsis Thaliana Network"
#' output: pdf_document
#' author: ""
#' date: Example for GeneNet 1.2.13 (August 2015) or later
#' ---

#' This note reproduces the “Arabidopsis thaliana” network example from
#' R. Opgen-Rhein and K. Strimmer. 2007. *From correlation to causation 
#' networks: a simple approximate learning algorithm and its application 
#' to high-dimensional plant gene expression data.*
#' BMC Syst. Biol. **1**: 37.
#' (http://dx.doi.org/10.1186/1752-0509-1-37)

#'  
#'The original source of the data is 
#' Smith et al. 2004. *Diurnal changes in the transcriptom encoding 
#' enzymes of starch metabolism provide evidence for both transcriptional
#' and posttranscriptional regulation of starch metabolism in Arabidopsis 
#' leaves.*  Plant Physiol. **136**: 2687-2699.
#'  

#' This example was suggested by
#' Papapit Ingkasuwan, Division of Biotechnology,
#' School of Bioresources and Technology,
#' King Mongkut's University of Technology Thonburi, Bangkok, Thailand.


#'
#' # Inspect Data

library("GeneNet")
data("arth800")
summary(arth800.expr)

#' Plot time series:
plot(arth800.expr, 1:9)

#' Inspect pairwise scatter plots:
panel.cor = function(x, y, digits=2, prefix="", cex.cor)
{
    usr = par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r = abs(cor(x, y))
    txt = format(c(r, 0.123456789), digits=digits)[1]
    txt = paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex = 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex * r)
}
#+ fig.width=8, fig.height=10
pairs(arth800.expr[,1:9], lower.panel=panel.smooth, upper.panel=panel.cor)


#'
#' # Compute Partial Correlations and Select Relevant Edges

pcor.dyn = ggm.estimate.pcor(arth800.expr, method = "dynamic")
arth.edges = network.test.edges(pcor.dyn,direct=TRUE)
dim(arth.edges)

#' We use the strongest 150 edges:
arth.net = extract.network(arth.edges, method.ggm="number", cutoff.ggm=150)


#'
#' # Construct Graph

library("graph")

node.labels = as.character(1:ncol(arth800.expr))
gr = network.make.graph( arth.net, node.labels, drop.singles=TRUE) 

#' Some information about the graph

#' Number of nodes:
num.nodes(gr)

#' Correlations:
edge.info(gr)$weight

#' Number of directed ("forward") and undirected ("none") edges:
table(  edge.info(gr)$dir )
# /*
# forward    none 
#     55      95
# */

#'
#' # Well-Connected Nodes

#' Nodes connected with many edges:
sort(node.degree(gr), decreasing=TRUE)[1:10]
# /*
# 570  81 783  47 422 558 452 539 738 272
# 20  17  10   9   9   9   8   8   8   7
# */

arth800.descr[570]
# /* [1] "AP2 transcription factor - like protein" */

arth800.descr[81]
# /* [1] "ATRPAC43; DNA binding / DNA-directed RNA polymerase;  */

arth800.descr[558]
# /* [1]"structural constituent of ribosome; */

arth800.descr[539]
# /* [1] "DNA binding / transcription factor;  */

arth800.descr[783]
# /* [1] "RNA binding / RNA methyltransferase; */



#'
#' # Plot Network

library("Rgraphviz")

#' For a more beautiful plot of the network set node and edge parameters:

#' Set global node and edge attributes:
globalAttrs = list()
globalAttrs$edge = list(color = "black", lty = "solid", lwd = 1, arrowsize=1)
globalAttrs$node = list(fillcolor = gray(.95), shape = "ellipse", fixedsize = FALSE)

#' Set attributes of some particular nodes:
nodeAttrs = list()
nodeAttrs$fillcolor = c('570' = "red", "81" = "red") # highlight hub nodes

#' Set edge attributes:
edi = edge.info(gr) # edge directions and correlations
edgeAttrs = list()
edgeAttrs$dir =  edi$dir # set edge directions 
cutoff = quantile(abs(edi$weight), c(0.2, 0.8)) # thresholds for line width / coloring
edgeAttrs$lty = ifelse(edi$weight < 0, "dotted", "solid") # negative correlation
edgeAttrs$color = ifelse( abs(edi$weight <= cutoff[1]), "grey", "black") # lower 20% quantile
edgeAttrs$lwd = ifelse(abs(edi$weight >= cutoff[2]), 2, 1) # upper 20% quantile

#+ fig.width=8, fig.height=7
plot(gr, attrs = globalAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs, "fdp")


