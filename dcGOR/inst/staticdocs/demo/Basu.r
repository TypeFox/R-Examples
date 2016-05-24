# This is a demo for analysing promiscuous Pfam domains by Basu et al
# 
# This list of promiscuous Pfam domains was reported in <a href="http://www.ncbi.nlm.nih.gov/pubmed/18230802" target="18230802">http://www.ncbi.nlm.nih.gov/pubmed/18230802</a>. There are 215 domains identified as strongly promiscuous (a tendency to occur in diverse domain architectures), in which 76 were taken from Pfam domains and thus used for this demo. It can be found in the file <a href="http://dcgor.r-forge.r-project.org/data/demo/domain_promiscuity_Basu_et_al_2008.txt">domain_promiscuity_Basu_et_al_2008.txt</a>, containing two columns: 1st column for Pfam domain ID, and 2nd column for promiscuity value (the higher the more promiscuous).
###############################################################################
library(dcGOR)

# Read promiscuous Pfam domains
domains <- read.delim("http://dcgor.r-forge.r-project.org/data/demo/domain_promiscuity_Basu_et_al_2008.txt",header=T)
domains[1:5,]

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Enrichment analysis for Pfam domains
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

## define the input data
data <- as.character(domains$Pfam)

### load Pfam domain informtion (as 'InfoDataFrame' object)
Pfam <- dcRDataLoader('Pfam')
Pfam

## 1) GOBP enrichment analysis, producing an object of S4 class 'Eoutput'
### By default, using all annotatable domains as the background
eoutput <- dcEnrichment(data, domain="Pfam", ontology="GOBP")
eoutput
### write into a local file <a href="Basu_GOBP_enrichments.txt">Basu_GOBP_enrichments.txt</a>
write(eoutput, file='Basu_GOBP_enrichments.txt')
### view the top 5 significant terms
view(eoutput, top_num=5, sortBy="pvalue", details=TRUE)
### visualise the top 4 significant terms (adjp < 0.05) in GOBP DAG
#### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
visEnrichment(eoutput, num_top_nodes=4, layout.orientation="top_bottom", zlim=c(0,4))
#### look at Pfam domains annotated by the most signficant term
tmp <- as.character(view(eoutput, top_num=1, sortBy="pvalue", details=T)$members)
tmp <- unlist(strsplit(tmp,","))
Data(Pfam)[match(tmp,rowNames(Pfam)),]

## 2) GOMF enrichment analysis, producing an object of S4 class 'Eoutput'
eoutput <- dcEnrichment(data, domain="Pfam", ontology="GOMF")
eoutput
### write into a local file <a href="Basu_GOMF_enrichments.txt">Basu_GOMF_enrichments.txt</a>
write(eoutput, file='Basu_GOMF_enrichments.txt')
### view the top 5 significant terms
view(eoutput, top_num=5, sortBy="pvalue", details=TRUE)
### visualise the top 5 significant terms in GOMF DAG
#### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
visEnrichment(eoutput, layout.orientation="top_bottom", zlim=c(0,4))
#### look at Pfam domains annotated by the most signficant term
tmp <- as.character(view(eoutput, top_num=1, sortBy="pvalue", details=T)$members)
tmp <- unlist(strsplit(tmp,","))
Data(Pfam)[match(tmp,rowNames(Pfam)),]

## 3) GOCC enrichment analysis, producing an object of S4 class 'Eoutput'
eoutput <- dcEnrichment(data, domain="Pfam", ontology="GOCC")
eoutput
### write into a local file <a href="Basu_GOCC_enrichments.txt">Basu_GOCC_enrichments.txt</a>
write(eoutput, file='Basu_GOCC_enrichments.txt')
### view the top 5 significant terms
view(eoutput, top_num=5, sortBy="pvalue", details=FALSE)
### visualise the top 5 significant terms in GOMF DAG
#### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
visEnrichment(eoutput)
#### look at Pfam domains annotated by the most signficant term
tmp <- as.character(view(eoutput, top_num=1, sortBy="pvalue", details=T)$members)
tmp <- unlist(strsplit(tmp,","))
Data(Pfam)[match(tmp,rowNames(Pfam)),]


#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Pair-wise semantic similarity between promiscuous Pfam domains
#---------------------------------------------------------------------
#---------------------------------------------------------------------

## 1) GOBP-based semantic similarity
### load onto.GOBP (as 'Onto' object)
g <- dcRDataLoader('onto.GOBP')
g
### load Pfam domains annotated by GOBP (as 'Anno' object)
Anno <- dcRDataLoader('Pfam2GOBP')
Anno
### prepare for ontology appended with annotation information
dag <- dcDAGannotate(g, annotations=Anno, path.mode="shortest_paths",verbose=FALSE)
dag
### calculate GOBP-based pair-wise semantic similarity between domains
dnetwork <- dcDAGdomainSim(g=dag, domains=as.character(domains$Pfam), method.domain="BM.average", method.term="Resnik", parallel=FALSE, verbose=TRUE)
dnetwork
### heatmap the adjacency matrix of the domain network
Adj_GOBP <- as.matrix(adjMatrix(dnetwork))
visHeatmapAdv(Adj_GOBP, Rowv=F, Colv=F, dendrogram="none", colormap="white-lightpink-darkred", zlim=c(0,1.5), cexRow=0.7, cexCol=0.7, KeyValueName="GOBP semantic similarity")

## 2) GOMF-based semantic similarity
### load onto.GOMF (as 'Onto' object)
g <- dcRDataLoader('onto.GOMF')
g
### load Pfam domains annotated by GOMF (as 'Anno' object)
Anno <- dcRDataLoader('Pfam2GOMF')
Anno
### prepare for ontology appended with annotation information
dag <- dcDAGannotate(g, annotations=Anno, path.mode="shortest_paths",verbose=FALSE)
dag
### calculate GOMF-based pair-wise semantic similarity between domains
dnetwork <- dcDAGdomainSim(g=dag, domains=as.character(domains$Pfam), method.domain="BM.average", method.term="Resnik", parallel=FALSE, verbose=TRUE)
dnetwork
### heatmap the adjacency matrix of the domain network
Adj_GOMF <- as.matrix(adjMatrix(dnetwork))
visHeatmapAdv(Adj_GOMF, Rowv=F, Colv=F, dendrogram="none", colormap="white-lightpink-darkred", zlim=c(0,1.5), cexRow=0.7, cexCol=0.7, KeyValueName="GOMF semantic similarity")

## 3) GOCC-based semantic similarity
### load onto.GOCC (as 'Onto' object)
g <- dcRDataLoader('onto.GOCC')
g
### load Pfam domains annotated by GOCC (as 'Anno' object)
Anno <- dcRDataLoader('Pfam2GOCC')
Anno
### prepare for ontology appended with annotation information
dag <- dcDAGannotate(g, annotations=Anno, path.mode="shortest_paths",verbose=FALSE)
dag
### calculate GOCC-based pair-wise semantic similarity between domains
dnetwork <- dcDAGdomainSim(g=dag, domains=as.character(domains$Pfam), method.domain="BM.average", method.term="Resnik", parallel=FALSE, verbose=TRUE)
dnetwork
### heatmap the adjacency matrix of the domain network
Adj_GOCC <- as.matrix(adjMatrix(dnetwork))
visHeatmapAdv(Adj_GOCC, Rowv=F, Colv=F, dendrogram="none", colormap="white-lightpink-darkred", zlim=c(0,1.5), cexRow=0.7, cexCol=0.7, KeyValueName="GOCC semantic similarity")

## 4) Obtain GO-based overall semantic similarity via merging all three subontology (GOBP, GOMF and GOCC) based semantic similarity
allnodes <- sort(unique(c(rownames(Adj_GOBP), rownames(Adj_GOMF), rownames(Adj_GOCC))))
D <- matrix(0, nrow=length(allnodes), ncol=length(allnodes))
colnames(D) <- rownames(D) <- allnodes
### add Adj_GOBP
ind <- match(rownames(Adj_GOBP), allnodes)
D[ind,ind] <- D[ind,ind]+Adj_GOBP
### add Adj_GOMF
ind <- match(rownames(Adj_GOMF), allnodes)
D[ind,ind] <- D[ind,ind]+Adj_GOMF
### add Adj_GOCC
ind <- match(rownames(Adj_GOCC), allnodes)
D[ind,ind] <- D[ind,ind]+Adj_GOCC
### heatmap the GO-based overall semantic similarity
visHeatmapAdv(D, Rowv=T, Colv=T, dendrogram="none", colormap="white-lightpink-darkred", zlim=c(0,2), cexRow=0.5, cexCol=0.5, KeyValueName="GO overall semantic similarity")
