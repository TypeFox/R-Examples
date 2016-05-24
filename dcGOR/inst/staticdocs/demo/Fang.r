# This is a demo for analysing ancestral superfamily domainome in Eukaryotes by Fang et al
# 
# This ancestral domain-ome dataset (available from <a href="http://www.ncbi.nlm.nih.gov/pubmed/23778980" target="23778980">http://www.ncbi.nlm.nih.gov/pubmed/23778980</a>) is stored as an 'Anno' object (S4 class). It contains information about domain repertoires (a complete set of domains: domain-ome) in Eukaryotes (including extant and ancestral genomes):
## annoData(Ancestral_domainome): a sparse matrix of 2019 domains X 875 terms/genomes (including 438 extant genomes and 437 ancestral genomes), with each entry telling how many different architectures a domain has in a genome. Note: zero entry also means that this domain is absent in the genome
## termData(Ancestral_domainome): variables describing terms/genomes (i.e. columns in annoData), including extant/ancestral genome information: "left_id" (unique and used as internal id), "right_id" (used in combination with "left_id" to define the post-ordered binary tree structure), "taxon_id" (NCBI taxonomy id, if matched), "genome" (2-letter genome identifiers used in SUPERFAMILY, if being extant), "name" (NCBI taxonomy name, if matched), "rank" (NCBI taxonomy rank, if matched), "branchlength" (branch length in relevance to the parent), and "common_name" (NCBI taxonomy common name, if matched and existed)
## domainData: variables describing domains (i.e. rows in annoData), including information about domains: "sunid" for SCOP id, "level" for SCOP level, "classification" for SCOP classification, "description" for SCOP description
###############################################################################
library(dcGOR)

#----------------------------------------------------------
#----------------------------------------------------------
# Derive domains unique/gained in human compared to Metazoa
#----------------------------------------------------------
#----------------------------------------------------------

# load data as an 'Anno' object
Ancestral_domainome <- dcRDataLoader("Ancestral_domainome")
Ancestral_domainome

# extract a list of domains that are present at Metazoa
flag_genome <- which(tData(Ancestral_domainome)$name=="Metazoa")
flag_domain <- annoData(Ancestral_domainome)[,flag_genome]!=0
domains_metazoa <- domainData(Ancestral_domainome)[flag_domain,]
domains_metazoa

# extract a list of domains that are present at human
flag_genome <- which(tData(Ancestral_domainome)$name=="Homo sapiens")
flag_domain <- annoData(Ancestral_domainome)[,flag_genome]!=0
domains_human <- domainData(Ancestral_domainome)[flag_domain,]
domains_human

# calculate domains unique in human
domains_human_unique <- setdiff(rowNames(domains_human), rowNames(domains_metazoa))
length(domains_human_unique)

# write into a local file <a href="Domains_unique_in_human.txt">Domains_unique_in_human.txt</a>
df <- dData(Ancestral_domainome)
ind <- match(domains_human_unique, rownames(df))
out <- df[ind,]
write.table(out, file="Domains_unique_in_human.txt", col.names=T, row.names=F, sep="\t")

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Enrichment analysis for domains unique/gained in human compared to Metazoa
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

## define the input data (domains gained in human) and the background (all domains in Metazoa)
data <- domains_human_unique
background <- rowNames(domains_metazoa)

## 1a) GOBP enrichment analysis, producing an object of S4 class 'Eoutput'
### By default, using all annotatable domains as the background
eoutput <- dcEnrichment(data, domain="SCOP.sf", ontology="GOBP")
eoutput
### write into a local file <a href="GOBP_enrichments.txt">GOBP_enrichments.txt</a>
write(eoutput, file='GOBP_enrichments.txt')
### view the top 5 significant terms
view(eoutput, top_num=5, sortBy="pvalue", details=FALSE)
### visualise the top 5 significant terms in GOMF DAG
#### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
visEnrichment(eoutput)

## 1b) GOBP enrichment analysis, producing an object of S4 class 'Eoutput'
### Alternatively, using all domains in Metazoa as the background (customised)
eoutput <- dcEnrichment(data, background, domain="SCOP.sf", ontology="GOBP")
eoutput
### write into a local file <a href="GOBP_enrichments.txt">GOBP_enrichments.txt</a>
write(eoutput, file='GOBP_enrichments_customised.txt')
### view the top 5 significant terms
view(eoutput, top_num=5, sortBy="pvalue", details=FALSE)
### visualise the top 5 significant terms in GOMF DAG
#### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
visEnrichment(eoutput)

visEnrichment(eoutput,height=10, width=10,graph.node.attrs=list(fontsize=20),node.attrs=list(fontsize=28), node.info = c("both", "none","term_id", "term_name", "full_term_name")[3])

## 2) GOMF enrichment analysis, producing an object of S4 class 'Eoutput'
eoutput <- dcEnrichment(data, background, domain="SCOP.sf", ontology="GOMF")
eoutput
### write into a local file <a href="GOMF_enrichments.txt">GOMF_enrichments.txt</a>
write(eoutput, file='GOMF_enrichments.txt')
### view the top 5 significant terms
view(eoutput, top_num=5, sortBy="pvalue", details=FALSE)
### visualise the top 5 significant terms in GOMF DAG
#### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
visEnrichment(eoutput)

## 3) HPPA enrichment analysis, producing an object of S4 class 'Eoutput'
eoutput <- dcEnrichment(data, domain="SCOP.sf", ontology="HPPA")
eoutput
### write into a local file <a href="HPPA_enrichments.txt">HPPA_enrichments.txt</a>
write(eoutput, file='HPPA_enrichments.txt')
### view the top 5 significant terms
view(eoutput, top_num=5, sortBy="pvalue", details=TRUE)
### visualise the top 5 significant terms in GOMF DAG
#### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
visEnrichment(eoutput, path.mode="all_paths", node.info="full_term_name")

## 4) DO enrichment analysis, producing an object of S4 class 'Eoutput'
eoutput <- dcEnrichment(data, domain="SCOP.sf", ontology="DO")
eoutput
### write into a local file <a href="DO_enrichments.txt">DO_enrichments.txt</a>
write(eoutput, file='DO_enrichments.txt')
### view the top 5 significant terms
view(eoutput, top_num=5, sortBy="pvalue", details=TRUE)
### visualise the top 5 significant terms in GOMF DAG
#### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
visEnrichment(eoutput, path.mode="all_paths")

#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Pair-wise semantic similarity between domains unique/gained in human
#---------------------------------------------------------------------
#---------------------------------------------------------------------

## 1) load onto.DO (as 'Onto' object)
g <- dcRDataLoader('onto.DO')
g
## 2) load SCOP superfamilies annotated by DO (as 'Anno' object)
Anno <- dcRDataLoader('SCOP.sf2DO')
Anno

## 3) prepare for ontology appended with annotation information
dag <- dcDAGannotate(g, annotations=Anno, path.mode="shortest_paths",verbose=FALSE)
dag

## 4) calculate pair-wise semantic similarity between domains unique/gained in human
domains <- domains_human_unique
dnetwork <- dcDAGdomainSim(g=dag, domains=domains, method.domain="BM.average", method.term="Resnik", parallel=FALSE, verbose=TRUE)
dnetwork

## 5) heatmap the adjacency matrix of the domain network
D <- as.matrix(adjMatrix(dnetwork))
supraHex::visHeatmapAdv(D, Rowv=F, Colv=F, dendrogram="none", colormap="white-lightpink-darkred", zlim=c(0,1.2), KeyValueName="DO semantic similarity")

## 6) visualise the domain network as a graph
### convert it to an object of class 'igraph' (for subsequent visualisation)
ig <- dcConverter(dnetwork, from='Dnetwork', to='igraph')
ig
### extract edge weight (with 2-digit precision)
x <- signif(E(ig)$weight, digits=2)
### rescale into an interval [1,4] as edge width
edge.width <- 1 + (x-min(x))/(max(x)-min(x))*3
### prepare the node labels (including domain id and description)
ind <- match(V(ig)$name,domainNames(Anno))
vertex.label <- paste(V(ig)$name, '\n', as.character(dData(Anno)[ind,]$description), sep='')
### do visualisation
dnet::visNet(g=ig, vertex.label=vertex.label, vertex.label.color="black", vertex.label.cex=0.7, vertex.shape="sphere", edge.width=edge.width, edge.label=x, edge.label.cex=0.7)

#------------------------------------------------------------
#------------------------------------------------------------
# Random Walk with Restart (RWR)-based contacts between terms
#------------------------------------------------------------
#------------------------------------------------------------

## 1) define sets of seeds: each seed with equal weight (i.e. all non-zero entries are '1')
### seeds are terms from GOMF
Anno <- dcRDataLoader('SCOP.sf2GOMF')
flag <- match(domainNames(Anno), nodeNames(dnetwork))
ind <- which(!is.na(flag))
data <- as.matrix(annoData(Anno)[ind,])
### GOMF terms having at least 3 annotatable domains that are also in human-unique domains above
ind <- apply(data,2,sum)>=3
data <- data[,ind]

## 2) assess their contact significance and strength based on RWR
coutput <- dcRWRpipeline(data=data, g=dnetwork, permutation="degree", num.permutation=2000, adjp.cutoff=0.1, parallel=FALSE)
coutput

# 3) write into the files in your local directory
## <a href="Coutput_adjp.txt">Coutput_adjp.txt</a> for contact significant adjusted p-values
write(coutput, file='Coutput_adjp.txt', saveBy="adjp")
## <a href="Coutput_zscore.txt">Coutput_zscore.txt</a> for contact strength z-scores
write(coutput, file='Coutput_zscore.txt', saveBy="zscore")
## <a href="Coutput_node_info.txt">Coutput_node_info.txt</a> for node info
df <- tData(Anno)
ind <- match(rownames(zscore(coutput)), rownames(df))
out <- df[ind,]
write.table(out, file="Coutput_node_info.txt", col.names=T, row.names=F, sep="\t")

## 4) extract contact network
cnet <- cnetwork(coutput)
cnet

## 5) visualise the contact network as a graph
### convert it to an object of class 'igraph' (for subsequent visualisation)
ig <- dcConverter(cnet, from='Cnetwork', to='igraph')
ig
### extract edge weight (with 2-digit precision)
x <- signif(E(ig)$weight, digits=2)
### rescale into an interval [1,4] as edge width
edge.width <- 1 + (x-min(x))/(max(x)-min(x))*3
### prepare the node labels (including term id and name)
ind <- match(V(ig)$name,termNames(Anno))
vertex.label <- paste(V(ig)$name, '\n', tData(Anno)[ind,]$Name, sep='')
### do visualisation
dnet::visNet(g=ig, vertex.label=vertex.label, vertex.label.color="blue", vertex.label.cex=0.7, vertex.shape="sphere", vertex.color="blue", edge.width=edge.width, edge.label=x, edge.label.cex=0.7)
