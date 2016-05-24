# This is a demo for multilayer replication timing dataset from Hiratani et al
# 
# This multilay omics dataset (available from <a href="http://www.ncbi.nlm.nih.gov/pubmed/19952138" target="19952138">http://www.ncbi.nlm.nih.gov/pubmed/19952138</a>) involves genome-wide replication-timing profiles of 22 cell lines from early mouse embryogenesis. These cell lines can be categorised into: 1) pluripotent cells, including ESCs (ESC_46C, ESC_D3 and ESC_TT2) and iPSCs (iPSC, iPSC_1D4 and iPSC_2D4); 2) partially-reprogrammed iPSCs (piPSC_1A2, piPSC_1B3 and piPSC_V3); 3) early epiblast (EPL and EMB3_D3); 4) late epiblast (EpiSC5 and EpiSC7); 5) Ectoderm (EBM6_D3, EBM9_D3, NPC_46C and NPC_TT2); 6) Mesoderm and Endoderm; and 7) late Mesoderm (Myoblast, MEF_female and MEF_male).
#
# The dataset is extracted for RefSeq gene TSS locations, including:
## RT: a replication timing matrix of 17,292 genes X 22 samples;
## CpG: a matrix of 17,292 genes X 1 containing gene additional information on promoter CpG classification (see <a href="http://www.ncbi.nlm.nih.gov/pubmed/17603471" target="17603471">http://www.ncbi.nlm.nih.gov/pubmed/17603471</a>), with '1' for HCP (high CpG density promoters), '-1' for LCP (low CpG density promoters), '0' for ICP (intermediate CpG density promoters), and 'NA' for unclassified;
## EX: an expression matrix of 17,292 genes X 8 samples, and samples include pluripotent cells (ESC_D3); early epiblast (EMB3_D3); late epiblast (EpiSC7); Ectoderm (EBM6_D3 and EBM9_D3); Mesoderm and Endoderm.
###############################################################################

# Load this multilayer dataset
load(url("http://dnet.r-forge.r-project.org/RData/1.0.7/Hiratani_TableS1.RData"))
ls() # you should see three variables: 'RT', 'CpG' and 'EX'

# Load the package 'dnet'
library(dnet)

# Load or/and install packages "Biobase" and "limma" that are specifically used in this demo
for(pkg in c("Biobase","limma")){
    if(!require(pkg, character.only=T)){
        source("http://bioconductor.org/biocLite.R")
        biocLite(pkg)
        lapply(pkg, library, character.only=T)
    }
}


# Here, we are interested to analyse replication timing data and their difference between different sample groups
# To this end, it is better to create the 'eset' object including sample grouping indication information
group <- c(rep("ESC",3), rep("iPSC",3), rep("eEpiblast",2), rep("lEpiblast",2), rep("Ectoderm",4), rep("Mesoderm",1), rep("Endoderm",1), rep("piPSC",3), rep("Myoblast",3))
pdata <- data.frame(group=group, row.names=colnames(RT))
esetGene <- new("ExpressionSet", exprs=as.matrix(RT), phenoData=as(pdata,"AnnotatedDataFrame"))
esetGene

# Look at the samples and their groups belonging to
pData(esetGene)

# Now, load the gene network in mouse
# As part of dnet package, this network has been prepared and stored as an igraph object
# The network is extracted from the STRING database (version 10). Only those associations with medium confidence (score>=400) are retained.
org.Mm.string <- dRDataLoader(RData='org.Mm.string')
org.Mm.string

# Look at the first 5 node information (gene symbols)
V(org.Mm.string)$symbol[1:5]

# Focus on the part of 'org.Mm.string' that only contains genes in 'esetGene'
ind <- match(V(org.Mm.string)$symbol, rownames(esetGene))
## this part of 'org.Mm.string' is called 'network'
nodes_mapped <- V(org.Mm.string)$name[!is.na(ind)]
network <- dNetInduce(g=org.Mm.string, nodes_query=nodes_mapped, knn=0, remove.loops=T, largest.comp=T)
V(network)$name <- V(network)$symbol
network

# Identification of gene-active subnetwork
# 1) obtain the information associated with nodes/genes, such as the p-value significance as node information 
# Here, we use the package 'limma' to identify differential Replication timing
## define the design matrix in an order manner
all <- as.vector(pData(esetGene)$group)
level <- levels(factor(all))
index_level <- sapply(level, function(x) which(all==x)[1])
level_sorted <- all[sort(index_level, decreasing=F)]
design <- sapply(level_sorted, function(x) as.numeric(all==x)) # Convert a factor column to multiple boolean columns
design
## define a contrast matrix: the pairwise comparisons of sample groups
contrasts <- dContrast(level_sorted, contrast.type="pairwise")
contrast.matrix <- makeContrasts(contrasts=contrasts$each, levels=design)
colnames(contrast.matrix) <- contrasts$name
colnames(contrast.matrix)
## a linear model is fitted for every gene by the function lmFit
fit <- lmFit(exprs(esetGene), design)
## computes moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
## for p-value
pvals <- as.matrix(fit2$p.value)
## for adjusted p-value
adjpvals <- sapply(1:ncol(pvals),function(x) {
    p.adjust(pvals[,x], method="BH")
})
colnames(adjpvals) <- colnames(pvals)
## num of differentially expressed genes
apply(adjpvals<1e-2, 2, sum)
## only for the comparisons of piPSC against iPSC
my_contrast <- "piPSC_iPSC"
## get the p-values and calculate the scores thereupon
pval <- pvals[,my_contrast]
## look at the distribution of p-values
hist(pval)

# 2) identification of gene-active subnetwork
## restrict the identified subnetwork to have the node size of 40 or so
#g <- dNetPipeline(g=network, pval=pval, nsize=40)
## corresponding to fdr=5.50e-07
g <- dNetPipeline(g=network, pval=pval, significance.threshold=5.50e-07)
g

# 3) visualisation of the gene-active subnetwork itself
## the layout of the network visualisation (fixed in different visuals) 
glayout <- layout.fruchterman.reingold(g)
## color nodes according to communities (identified via a spin-glass model and simulated annealing)
com <- spinglass.community(g, spins=25)
com$csize <- sapply(1:length(com),function(x) sum(com$membership==x))
vgroups <- com$membership
colormap <- "yellow-darkorange"
palette.name <- visColormap(colormap=colormap)
mcolors <- palette.name(length(com))
vcolors <- mcolors[vgroups]
com$significance <- dCommSignif(g, com)
## node sizes according to degrees
vdegrees <- igraph::degree(g)
## highlight different communities
mark.groups <- communities(com)
mark.col <- visColoralpha(mcolors, alpha=0.2)
mark.border <- visColoralpha(mcolors, alpha=0.2)
edge.color <- c("grey", "black")[crossing(com,g)+1]
## visualise the subnetwrok
visNet(g, glayout=glayout, vertex.label=V(g)$geneSymbol, vertex.color=vcolors, vertex.frame.color=vcolors, vertex.shape="sphere", mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

# 4) visualisation of the gene-active subnetwork overlaid by the node/gene score
max_colorbar <- ceiling(quantile(abs(V(g)$score),0.75))
visNet(g, glayout=glayout, pattern=V(g)$score, zlim=c(-1*max_colorbar,max_colorbar), vertex.shape="circle")

# 5) visualisation of the gene-active subnetwork overlaid by the differential replication timing
colormap <- "darkgreen-lightgreen-lightpink-darkred"
logFC <- fit2$coefficients[V(g)$name,my_contrast]
visNet(g, glayout=glayout, pattern=logFC, colormap=colormap, vertex.shape="circle")

# 6) Network-based sample classifications and visualisations on 2D sample landscape
# it uses the gene-active subnetwork overlaid by all replication timing data
data <- exprs(esetGene)[V(g)$name,]
sReorder <- dNetReorder(g, data, feature="edge", node.normalise="degree", amplifier=3, metric="none")
visNetReorder(g=g, data=data, sReorder=sReorder, height=ceiling(sqrt(ncol(data)))*2, newpage=T, glayout=glayout, colormap=colormap, vertex.label=NA,vertex.shape="sphere", vertex.size=16,mtext.cex=0.4,border.color="888888")

# 7) heatmap of replication timing data in the subnetwork
visHeatmapAdv(data, colormap=colormap, KeyValueName="log2(Early/Late)")

# 8) output the subnetwork and their replication timing data
## Write the subnetwork into a SIF-formatted file (Simple Interaction File)
sif <- data.frame(source=get.edgelist(g)[,1], type="interaction", target=get.edgelist(g)[,2])
write.table(sif, file=paste(my_contrast,".sif", sep=""), quote=F, row.names=F,col.names=F,sep="\t")
## Output the corresponding replication timing data
hmap <- data.frame(Symbol=rownames(data), data)
write.table(hmap, file=paste(my_contrast,".txt", sep=""), quote=F, row.names=F,col.names=T,sep="\t")

# 9) enrichment analysis for genes in the subnetwork
## get a list of genes in the subnetwork
data <- V(g)$name
data

## 9a) GOBP enrichment analysis
eTerm <- dEnricher(data, identity="symbol", genome="Mm", ontology="GOBP")
## write into the file called 'enrichment_GOBP.txt'
output <- dEnricherView(eTerm, top_num=NULL, sortBy="adjp", details=TRUE)
write.table(output, file="enrichment_GOBP.txt", quote=F, row.names=F,col.names=T,sep="\t")
## visualise the top significant terms in the GOBP heirarchy
## first, load the GOBP ontology
ig.GOBP <- dRDataLoader(RData='ig.GOBP')
g <- ig.GOBP
## select the top most significant 10 terms
top <- dEnricherView(eTerm, top_num=10, details=TRUE)
top
nodes_query <- rownames(top)
nodes.highlight <- rep("red", length(nodes_query))
names(nodes.highlight) <- nodes_query
## induce the shortest paths (one for each significant term) to the ontology root
subg <- dDAGinduce(g, nodes_query, path.mode="shortest_paths")
## color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
visDAG(g=subg, data=-1*log10(eTerm$adjp[V(subg)$name]), node.info="both", node.attrs=list(color=nodes.highlight))

## 9b) GOMF enrichment analysis
eTerm <- dEnricher(data, identity="symbol", genome="Mm", ontology="GOMF")
## write into the file called 'enrichment_GOMF.txt'
output <- dEnricherView(eTerm, top_num=NULL, sortBy="adjp", details=TRUE)
write.table(output, file="enrichment_GOMF.txt", quote=F, row.names=F,col.names=T,sep="\t")
## visualise the top significant terms in the GOMF heirarchy
## first, load the GOMF ontology
ig.GOMF <- dRDataLoader(RData='ig.GOMF')
g <- ig.GOMF
## select the top most significant 10 terms
top <- dEnricherView(eTerm, top_num=10, details=TRUE)
top
nodes_query <- rownames(top)
nodes.highlight <- rep("red", length(nodes_query))
names(nodes.highlight) <- nodes_query
## induce the shortest paths (one for each significant term) to the ontology root
subg <- dDAGinduce(g, nodes_query, path.mode="shortest_paths")
## color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
visDAG(g=subg, data=-1*log10(eTerm$adjp[V(subg)$name]), node.info="both", node.attrs=list(color=nodes.highlight))

## 9c) MP enrichment analysis
eTerm <- dEnricher(data, identity="symbol", genome="Mm", ontology="MP")
## write into the file called 'enrichment_MP.txt'
output <- dEnricherView(eTerm, top_num=NULL, sortBy="adjp", details=TRUE)
write.table(output, file="enrichment_MP.txt", quote=F, row.names=F,col.names=T,sep="\t")
## visualise the top significant terms in the MP heirarchy
## first, load the MP ontology
ig.MP <- dRDataLoader(RData='ig.MP')
g <- ig.MP
## select the top most significant 10 terms
top <- dEnricherView(eTerm, top_num=10, details=TRUE)
top
nodes_query <- rownames(top)
nodes.highlight <- rep("red", length(nodes_query))
names(nodes.highlight) <- nodes_query
## induce all possible paths to the ontology root
subg <- dDAGinduce(g, nodes_query)
## color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
visDAG(g=subg, data=-1*log10(eTerm$adjp[V(subg)$name]), node.info=c("none","term_id","term_name","both","full_term_name")[5], layout.orientation=c("left_right","top_bottom","bottom_top","right_left")[1], node.attrs=list(color=nodes.highlight))

## 9d) DO enrichment analysis
eTerm <- dEnricher(data, identity="symbol", genome="Mm", ontology="DO")
## write into the file called 'enrichment_DO.txt'
output <- dEnricherView(eTerm, top_num=NULL, sortBy="adjp", details=TRUE)
write.table(output, file="enrichment_DO.txt", quote=F, row.names=F,col.names=T,sep="\t")
## visualise the top significant terms in the DO heirarchy
## first, load the DO ontology
ig.DO <- dRDataLoader(RData='ig.DO')
g <- ig.DO
## select the top most significant 10 terms
top <- dEnricherView(eTerm, top_num=10, details=TRUE)
top
nodes_query <- rownames(top)
nodes.highlight <- rep("red", length(nodes_query))
## induce all possible shortest paths to the ontology root
subg <- dDAGinduce(g, nodes_query)
## color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
visDAG(g=subg, data=-1*log10(eTerm$adjp[V(subg)$name]), node.info="both", zlim=c(0,4), node.attrs=list(color=nodes.highlight))

## 9e) PS enrichment analysis
## use all common ancestors
eTerm <- dEnricher(data, identity="symbol", genome="Mm", ontology="PS", sizeRange=c(10,20000), min.overlap=0)
output <- dEnricherView(eTerm, top_num=NULL, sortBy="none", details=TRUE)
write.table(output, file="enrichment_PS.txt", quote=F, row.names=F,col.names=T,sep="\t")
output
## use common ancestors collapsed onto the known NCBI taxonomy
eTerm <- dEnricher(data, identity="symbol", genome="Mm", ontology="PS2", sizeRange=c(10,20000), min.overlap=0)
output <- dEnricherView(eTerm, top_num=NULL, sortBy="none", details=TRUE)
write.table(output, file="enrichment_PS2.txt", quote=F, row.names=F,col.names=T,sep="\t")
output
