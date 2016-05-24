# This is a demo for multilayer omics dataset from Hiratani et al
# 
# This multilay omics dataset (available from <a href="http://www.ncbi.nlm.nih.gov/pubmed/19952138" target="19952138">http://www.ncbi.nlm.nih.gov/pubmed/19952138</a>) involves genome-wide replication-timing profiles of 22 cell lines from early mouse embryogenesis. These cell lines can be categorised into: 1) pluripotent cells, including ESCs (ESC_46C, ESC_D3 and ESC_TT2) and iPSCs (iPSC, iPSC_1D4 and iPSC_2D4); 2) partially-reprogrammed iPSCs (piPSC_1A2, piPSC_1B3 and piPSC_V3); 3) early epiblast (EPL and EMB3_D3); 4) late epiblast (EpiSC5 and EpiSC7); 5) Ectoderm (EBM6_D3, EBM9_D3, NPC_46C and NPC_TT2); 6) Mesoderm and Endoderm; and 7) late Mesoderm (Myoblast, MEF_female and MEF_male).
#
# The dataset is extracted for RefSeq gene TSS locations, including:
## RT: a replication timing matrix of 17,292 genes X 22 samples;
## CpG: a matrix of 17,292 genes X 1 containing gene additional information on promoter CpG classification (see <a href="http://www.ncbi.nlm.nih.gov/pubmed/17603471" target="17603471">http://www.ncbi.nlm.nih.gov/pubmed/17603471</a>), with '1' for HCP (high CpG density promoters), '-1' for LCP (low CpG density promoters), '0' for ICP (intermediate CpG density promoters), and 'NA' for unclassified;
## EX: an expression matrix of 17,292 genes X 8 samples, and samples include pluripotent cells (ESC_D3); early epiblast (EMB3_D3); late epiblast (EpiSC7); Ectoderm (EBM6_D3 and EBM9_D3); Mesoderm and Endoderm.
###############################################################################
library(dnet)

# Load or install packages specifically used in this demo
if(!require(affy)){
    install.packages("affy",repos="http://www.stats.bris.ac.uk/R")
    library(affy)
}
if(!require(limma)){
    install.packages("limma",repos="http://www.stats.bris.ac.uk/R")
    library(limma)
}

# This multilayer omics dataset involves the information on DNA replication timing, promoter CpG classification and gene expression. It consists of digitised replication timing, promoter CpG status and expression levels of 17,292 genes in a variety of samples.
load(url("http://dnet.r-forge.r-project.org/data/Datasets/Hiratani_TableS1.RData"))
ls() # you should see three variables: 'RT', 'CpG' and 'EX'

# Entrez Gene information for the mouse
load(url("http://dnet.r-forge.r-project.org/data/Mm/org.Mm.eg.RData"))
gene_info <- org.Mm.eg$gene_info
gene_info[1:2,]

# create the 'eset' object for RT
group <- c(rep("ESC",3), rep("iPSC",3), rep("eEpiblast",2), rep("lEpiblast",2), rep("Ectoderm",4), rep("Mesoderm",1), rep("Endoderm",1), rep("piPSC",3), rep("Myoblast",3))
pdata <- data.frame(group=group, row.names=colnames(RT))
esetGene <- new("ExpressionSet", exprs=as.matrix(RT), phenoData=as(pdata,"AnnotatedDataFrame"))
esetGene

# An igraph object that contains a functional protein association network in mouse. The network is extracted from the STRING database (version 9.0.5). Only those associations with medium confidence (score>=0.4) are retained.
load(url("http://dnet.r-forge.r-project.org/data/Mm/org.Mm.string.RData"))
org.Mm.string

# extract network that only contains genes in esetGene
ind <- match(V(org.Mm.string)$symbol, rownames(esetGene))
## for extracted expression
esetGeneSub <- esetGene[ind[!is.na(ind)],]
esetGeneSub
## for extracted graph
nodes_mapped <- V(org.Mm.string)$name[!is.na(ind)]
network <- dNetInduce(g=org.Mm.string, nodes_query=nodes_mapped, knn=0, remove.loops=T, largest.comp=T)
V(network)$name <- V(network)$symbol
network

# 1) preparation of node p-values
## define the design matrix in a order manner
all <- as.vector(pData(esetGene)$group)
level <- levels(factor(all))
index_level <- sapply(level, function(x) which(all==x)[1])
level_sorted <- all[sort(index_level, decreasing=F)]
design <- sapply(level_sorted, function(x) as.numeric(all==x)) # Convert a factor column to multiple boolean columns
## a linear model is fitted for every gene by the function lmFit
fit <- lmFit(exprs(esetGene), design)
## define a contrast matrix
contrasts <- dContrast(level_sorted, contrast.type="pairwise")
contrast.matrix <- makeContrasts(contrasts=contrasts$each, levels=design)
colnames(contrast.matrix) <- contrasts$name
## computes moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
## for p-value
pvals <- as.matrix(fit2$p.value)
# for adjusted p-value
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
#pval <- dPvalAggregate(pvals, method="orderStatistic", order=ncol(pvals))
#pval <- dPvalAggregate(pvals, method="fishers", order=ncol(pvals))

# 2) identification of module
g <- dNetPipeline(g=network, pval=pval, nsize=40)
glayout <- layout.fruchterman.reingold(g)

# 3) color nodes according to communities identified via a spin-glass model and simulated annealing
#com <- spinglass.community(g, spins=3)
com <- walktrap.community(g, modularity=T)
com$csize <- sapply(1:length(com),function(x) sum(com$membership==x))
vgroups <- com$membership
colormap <- "yellow-darkorange"
palette.name <- visColormap(colormap=colormap)
mcolors <- palette.name(length(com))
vcolors <- mcolors[vgroups]
com$significance <- dCommSignif(g, com)

# 4) size nodes according to degrees
vdegrees <- igraph::degree(g)

# 5) sort nodes: first by communities and then degrees
tmp <- data.frame(ind=1:vcount(g), vgroups, vdegrees)
ordering <- tmp[order(vgroups,vdegrees),]$ind

# 6) visualise graph using 1-dimensional arc diagram
visNetArc(g, ordering=ordering, labels=V(g)$geneSymbol, vertex.label.color=vcolors, vertex.color=vcolors, vertex.frame.color=vcolors, vertex.size=log(vdegrees)+0.1, vertex.label.cex=0.4)

# 7) visualise graph using circle diagram
# 7a) drawn into a single circle 
visNetCircle(g=g, com=com, ordering=ordering, colormap=colormap, vertex.label=V(g)$symbol, vertex.size=igraph::degree(g)+5, vertex.label.color="black", vertex.label.cex=0.6, vertex.label.dist=0.75, vertex.shape="sphere", edge.color.within="grey", edge.color.crossing="black", edge.width=1, edge.lty=1, mark.shape=1, mark.expand=10)
# 7b) drawn into multiple circles
visNetCircle(g=g, com=com, circles="multiple", ordering=ordering, colormap=colormap, vertex.label=V(g)$symbol, vertex.size=igraph::degree(g)+5, vertex.label.color="black", vertex.label.cex=0.6, vertex.label.dist=0.25, vertex.shape="sphere", edge.color.within="grey", edge.color.crossing="black", edge.width=1, edge.lty=1, mark.shape=1, mark.expand=10)

# 8) as comparison, also visualise graph on 2-dimensional layout 
mark.groups <- communities(com)
mark.col <- visColoralpha(mcolors, alpha=0.2)
mark.border <- visColoralpha(mcolors, alpha=0.2)
edge.color <- c("grey", "black")[crossing(com,g)+1]
visNet(g, glayout=glayout, vertex.label=V(g)$geneSymbol, vertex.color=vcolors, vertex.frame.color=vcolors, vertex.shape="sphere", mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

legend_name <- paste("C",1:length(mcolors)," (n=",com$csize,", pval=",signif(com$significance,digits=2),")",sep='')
legend("bottomleft", legend=legend_name, fill=mcolors, bty="n", cex=0.6)

# 9) color by score and FC
# colored by score
visNet(g, glayout=glayout, pattern=V(g)$score, zlim=c(-1*ceiling(max(abs(V(g)$score))),ceiling(max(abs(V(g)$score)))), vertex.shape="circle", mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)
# colored by FC
colormap <- "darkgreen-lightgreen-lightpink-darkred"
logFC <- fit2$coefficients[V(g)$name,my_contrast]
visNet(g, glayout=glayout, pattern=logFC, colormap=colormap, vertex.shape="circle", mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

# 10) color by additional data
data <- exprs(esetGene)[V(g)$name,]
visNetMul(g=g, data=data, height=ceiling(sqrt(ncol(data)))*2, newpage=T,glayout=glayout,colormap=colormap,vertex.label=NA,vertex.shape="sphere", vertex.size=16,mtext.cex=0.7,border.color="888888", mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

# 11) color by additional data (be reordered)
sReorder <- dNetReorder(g, data, feature="edge", node.normalise="degree", amplifier=3, metric="none")
visNetReorder(g=g, data=data, sReorder=sReorder, height=ceiling(sqrt(ncol(data)))*2, newpage=T, glayout=glayout, colormap=colormap, vertex.label=NA,vertex.shape="sphere", vertex.size=16,mtext.cex=0.4,border.color="888888", mark.groups=mark.groups, mark.col=mark.col, mark.border=NA, mark.shape=1, mark.expand=10, edge.color=edge.color)

# 12) heatmap of subnetwork
visHeatmapAdv(data, colormap=colormap)
hmap <- data.frame(Symbol=rownames(data), data)
write.table(hmap, file=paste(my_contrast,".txt", sep=""), quote=F, row.names=F,col.names=T,sep="\t")

# 13) Write the subnetwork into a SIF-formatted file (Simple Interaction File)
sif <- data.frame(source=get.edgelist(g)[,1], type="interaction", target=get.edgelist(g)[,2])
write.table(sif, file=paste(my_contrast,".sif", sep=""), quote=F, row.names=F,col.names=F,sep="\t")
