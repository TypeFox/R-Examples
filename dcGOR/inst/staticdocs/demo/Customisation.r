# This is a demo for customised analysis using InterPro2GO mapping
# 
# InterPro2GO mapping is derived from manual annotations of InterPro from GO (<a href="http://www.ncbi.nlm.nih.gov/pubmed/22301074" target="22301074">http://www.ncbi.nlm.nih.gov/pubmed/22301074</a>.
# This demo assumes that the user has collected three bits (files) of information:
## 1) InterPro domains (see <a href="http://dcgor.r-forge.r-project.org/data/InterPro/InterPro.txt">InterPro.txt</a>);
## 2) GO Molecular Function (see <a href="http://dcgor.r-forge.r-project.org/data/onto/igraph_GOMF_nodes.txt">igraph_GOMF_nodes.txt</a> for GOMF term information and <a href="http://dcgor.r-forge.r-project.org/data/onto/igraph_GOBP_edges.txt">igraph_GOMF_edges.txt</a> for parent-child relations between GOMF terms);
## 3) Annotations of InterPro domains by GOMF terms (see <a href="http://dcgor.r-forge.r-project.org/data/InterPro/GO.txt">GO.txt</a> for all GO term information and <a href="http://dcgor.r-forge.r-project.org/data/InterPro/Domain2GOMF.txt">Domain2GOMF.txt</a> for annotations)
#
# From this demo, the user can learn how to customise analysis based on their own data (domains, ontology and annotations). The key to this customised analysis is to prepare input files exactly the same as those mentioned above. 
###############################################################################
library(dcGOR)

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Build (domain, ontology and annotations) objects from user-input files
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

## 1) InterPro domains
domain <- dcBuildInfoDataFrame(input.file="http://dcgor.r-forge.r-project.org/data/InterPro/InterPro.txt", output.file="domain.RData")
domain
## 2) ontology (GO Molecular Function, GOMF)
g <- dcBuildOnto(relations.file="http://dcgor.r-forge.r-project.org/data/onto/igraph_GOMF_edges.txt", nodes.file="http://dcgor.r-forge.r-project.org/data/onto/igraph_GOMF_nodes.txt", output.file="ontology.RData")
g
## 3) annotations (between InterPro domains and GOMF terms)
Anno <- dcBuildAnno(domain_info.file="http://dcgor.r-forge.r-project.org/data/InterPro/InterPro.txt", term_info.file="http://dcgor.r-forge.r-project.org/data/InterPro/GO.txt", association.file="http://dcgor.r-forge.r-project.org/data/InterPro/Domain2GOMF.txt", output.file="annotations.RData")
Anno
## In your working directory, you should see three RData-formatted files: "domain.RData", "ontology.RData", "annotations.RData"

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Enrichment analysis
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

## prepare data and background for enrichment analysis
### randomly select 100 domains as a list of domains of interest
data <- sample(rowNames(domain), 100)
length(data)
### randomly select 1000 domains as background
background <- sample(rowNames(domain), 1000)
length(background)

## perform enrichment analysis, producing an object of S4 class 'Eoutput'
eoutput <- dcEnrichment(data, background=background, domain.RData='domain.RData', ontology.RData='ontology.RData', annotations.RData='annotations.RData')
eoutput

## view the top 10 significance terms 
view(eoutput, top_num=10, sortBy="pvalue", details=FALSE)
### visualise the top 10 significant terms in the ontology hierarchy
### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
visEnrichment(eoutput)

#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Pair-wise semantic similarity between domains
#---------------------------------------------------------------------
#---------------------------------------------------------------------

## prepare for ontology appended with annotation information
dag <- dcDAGannotate(g, annotations=Anno, path.mode="shortest_paths", verbose=FALSE)
dag
## calculate pair-wise semantic similarity between 8 randomly chosen domains
alldomains <- unique(unlist(nInfo(dag)$annotations))
domains <- sample(alldomains,8)
dnetwork <- dcDAGdomainSim(g=dag, domains=domains, method.domain="BM.average", method.term="Resnik", parallel=FALSE, verbose=TRUE)
dnetwork
## convert it to an object of class 'igraph'
ig <- dcConverter(dnetwork, from='Dnetwork', to='igraph')
ig
## visualise the domain network
### extract edge weight (with 2-digit precision)
x <- signif(E(ig)$weight, digits=2)
### rescale into an interval [1,4] as edge width
edge.width <- 1 + (x-min(x))/(max(x)-min(x))*3
### do visualisation
dnet::visNet(g=ig, vertex.shape="sphere", edge.width=edge.width, edge.label=x, edge.label.cex=0.7)
