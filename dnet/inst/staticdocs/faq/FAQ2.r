# <span style="font-weight:bold; color:#F87217; text-decoration:underline">These built-in data are the backend of various analytical utilities supported in the dnet package, spanning a wide range of the known gene-centric knowledge across well-studied organisms.</span> They are provided as RData-formatted files which are regularly updated. Also, we will populate them by adding new knowledge, for example, upon request by users. The built-in RData are summarised in brief and available in the <a href="http://supfam.org/dnet/rdata.html">RData</a> page.

# Usually, the users do not need to download them by self for use. Instead, the users are encouraged to understand what they want to use by simply looking up the keywords in the <a href="http://supfam.org/dnet/docs.html">Documentations</a> page. The package has functions to import them or deal with them directly.

# <span style="font-weight:bold; color:#F87217; text-decoration:underline">The function dRDataLoader allows the users to import what they want to use.</span>

# <span style="font-weight:bold; color:#F87217; text-decoration:underline">For the ease to use, organism-specific data start with 'org', followed by the specific organim ('Hs' for human), and the data content: only 'eg' means information about Entrez Genes, and further appendix (for example, 'GOBP') means information about their annotations by Gene Ontology Biological Process (GOBP).</span>

## load human Entrez Genes (EG), and list the first 3 genes
org.Hs.eg <- dRDataLoader(RData='org.Hs.eg')
org.Hs.eg$gene_info[1:3,]
## load annotations of human Entrez Genes by Gene Ontology Biological Process (GOBP), inspect the content and list the first 3 terms
org.Hs.egGOBP <- dRDataLoader(RData='org.Hs.egGOBP')
names(org.Hs.egGOBP)
org.Hs.egGOBP$set_info[1:3,]
## load annotations of human Entrez Genes by Disease Ontology (DO), inspect the content and list the first 5 terms
org.Hs.egDO <- dRDataLoader(RData='org.Hs.egDO')
org.Hs.egDO$set_info[1:3,]
## load phylostratific age (PS) information of human Entrez Genes, inspect the content and list all our ancestors
org.Hs.egPS <- dRDataLoader(RData='org.Hs.egPS')
org.Hs.egPS$set_info
## load domain superfamily (SF) information of human Entrez Genes, inspect the content and list the first 3 superfamilies
org.Hs.egSF <- dRDataLoader(RData='org.Hs.egSF')
org.Hs.egSF$set_info[1:3,]
## load KEGG pathways for human Entrez Genes, inspect the content and list the first 3 pathways
org.Hs.egMsigdbC2KEGG <- dRDataLoader(RData='org.Hs.egMsigdbC2KEGG')
org.Hs.egMsigdbC2KEGG$set_info[1:3,]
## load the network for human Entrez Genes as an 'igraph' object
org.Hs.string <- dRDataLoader(RData='org.Hs.string')
org.Hs.string
## This network is extracted from the STRING database. Only those associations with medium confidence (score>=400) are retained. And the users can restrict to those edges with high confidence (score>=700, for example)
network <- igraph::subgraph.edges(org.Hs.string, eids=E(org.Hs.string)[combined_score>=700])
network

# <span style="font-weight:bold; color:#F87217; text-decoration:underline">In addition to data import, the package has also functions (see below) to deal with them directly. In these functions, the users only need to specify which genome/organism and which ontology to use.</span>

# Here, we use human TCGA mutation dataset as an example
data(TCGA_mutations)
symbols <- as.character(fData(TCGA_mutations)$Symbol)

## Enrichment analysis using Disease Ontology (DO)
data <- symbols[1:100] # select the first 100 human genes
eTerm <- dEnricher(data, identity="symbol", genome="Hs", ontology="DO")

## gene set enrichment analysis (GSEA) using KEGG pathways
tol <- apply(exprs(TCGA_mutations), 1, sum) # calculate the total mutations for each gene
data <- data.frame(tol=tol)
eTerm <- dGSEA(data, identity="symbol", genome="Hs", ontology="MsigdbC2KEGG")
