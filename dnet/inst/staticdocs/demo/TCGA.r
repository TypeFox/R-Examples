# This is a demo for TCGA mutational profile dataset from Kandoth et al
# 
# This dataset is available from TCGA (see <a href="http://www.ncbi.nlm.nih.gov/pubmed/24132290" target="24132290">http://www.ncbi.nlm.nih.gov/pubmed/24132290</a>), containing somatic mutational profiles for 3096 cancer samples with survival data. These cancer samples belong to one of 12 major cancer types, including breast adenocarcinoma (BRCA), lung adenocarcinoma (LUAD), lung squamous cell carcinoma (LUSC), uterine corpus endometrial carcinoma (UCEC), glioblastoma multiforme (GBM), head and neck squamous cell carcinoma (HNSC), colon and rectal carcinoma (COAD/READ), bladder urothelial carcinoma (BLCA), kidney renal clear cell carcinoma (KIRC), ovarian serous carcinoma (OV) and acute myeloid leukaemia (LAML). For each patient sample, somatic mutations are represented as a profile of  states on genes, where non-zero entry indicates a gene for which how many mutations have occurred in the tumor relative to germ line. The dataset is provided as an 'ExpressionSet' object.
## assayData: exprs(TCGA_mutations), a matrix of 19171 genes X 3096 samples;
## phenoData: pData(TCGA_mutations), variables describing sample phenotypes (i.e. columns in assayData), including clinical/survival information about samples: "time" (i.e. survival time in days), "status" (i.e., survival status: 0=alive; 1=dead), "Age" (the patient age in years), "Gender" (the patient gender: male/female), "TCGA_tumor_type", "Tumor_stage", "Tumor_grade"
## featureData: fData(TCGA_mutations), variables describing features (i.e. rows in assayData), including information about features/genes: "EntrezID" for gene EntrezID, "Symbol" for gene symbol, "Desc" for gene description, "Synonyms" for gene symbol alias
###############################################################################
library(dnet)

# Load or install packages specifically used in this demo
for(pkg in c("Biobase","survival")){
    if(!require(pkg, character.only=T)){
        source("http://bioconductor.org/biocLite.R")
        biocLite(pkg)
        lapply(pkg, library, character.only=T)
    }
}

# load an "ExpressionSet" object
TCGA_mutations <- dRDataLoader(RData='TCGA_mutations')
eset <- TCGA_mutations
# extract information about phenotype data
pd <- pData(eset)
pd[1:3,]
# extract information about feature/gene data
fd <- fData(eset)
fd[1:3,]
# extract information about mutational data
md <- exprs(eset)
md[1:3,1:3]
# number of samples for each cancer type
tumor_type <- sort(unique(pData(eset)$TCGA_tumor_type))
table(pData(eset)$TCGA_tumor_type)

# Survival analysis across tumor types using Cox proportional hazards model
# Cox regression yields an equation for the hazard/risk as a function of several explanatory variables

## fit a Cox proportional hazards model for age, gender, tumor type
data <- pd
fit <- survival::coxph(Surv(time,status) ~ Age + Gender + TCGA_tumor_type, data=data)
res <- anova(fit, test="Chisq")

## Now with gene mutational data in subject, adjust for other explanatory variables (or called covariates) include: age, gender, and tumor type
## only those genes with mutations at least 1% of samples will be analysed
flag <- sapply(1:nrow(md), function(i) ifelse(sum(md[i,]!=0)>=0.01*ncol(md), T, F))
esetGene <- eset[flag, ]
md_selected <- md[flag,]
## survival analysis to obtain hazard ratio (HR) and pvaules
HR <- rep(1, nrow(md_selected))
pvals <- rep(1, nrow(md_selected))
for(i in 1:nrow(md_selected)){
    ## fit a Cox proportional hazards model
    data <- cbind(pd, gene=md_selected[i,])
    fit <- survival::coxph(formula=Surv(time,status) ~ Age + Gender + TCGA_tumor_type + gene, data=data)
    ## ANOVA (Chisq test)
    res <- as.matrix(anova(fit, test="Chisq"))
    HR[i] <- res[5,2]
    pvals[i] <- res[5,4]
}
names(HR) <- rownames(md_selected)
names(pvals) <- rownames(md_selected)

# An igraph object that contains a functional protein association network in human. The network is extracted from the STRING database (version 10). Only those associations with medium confidence (score>=400) are retained.
org.Hs.string <- dRDataLoader(RData='org.Hs.string')
# restrict to those edges with high confidence (score>=700)
network <- subgraph.edges(org.Hs.string, eids=E(org.Hs.string)[combined_score>=700])
network

# extract network that only contains genes in pvals
ind <- match(V(network)$symbol, names(pvals))
## for extracted graph
nodes_mapped <- V(network)$name[!is.na(ind)]
network <- dNetInduce(g=network, nodes_query=nodes_mapped, knn=0, remove.loops=F, largest.comp=T)
V(network)$name <- V(network)$symbol
network

# Identification of gene-active network
net <- dNetPipeline(g=network, pval=pvals, method="customised", significance.threshold=2.5e-02)
net

# visualisation of the gene-active network itself
## the layout of the network visualisation (fixed in different visuals) 
glayout <- layout.fruchterman.reingold(net)
## color nodes according to communities (identified via a spin-glass model and simulated annealing)
com <- spinglass.community(net, spins=25)
com$csize <- sapply(1:length(com),function(x) sum(com$membership==x))
vgroups <- com$membership
colormap <- "yellow-darkorange"
palette.name <- visColormap(colormap=colormap)
mcolors <- palette.name(length(com))
vcolors <- mcolors[vgroups]
com$significance <- dCommSignif(net, com)
## node sizes according to degrees
vdegrees <- igraph::degree(net)
## highlight different communities
mark.groups <- communities(com)
mark.col <- visColoralpha(mcolors, alpha=0.2)
mark.border <- visColoralpha(mcolors, alpha=0.2)
edge.color <- c("#C0C0C0", "#000000")[crossing(com,net)+1]
edge.color <- visColoralpha(edge.color, alpha=0.5)
## visualise the subnetwrok
visNet(g=net, glayout=glayout, vertex.label=V(net)$geneSymbol, vertex.color=vcolors, vertex.frame.color=vcolors, vertex.shape="sphere", mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color, newpage=F, vertex.label.color="blue", vertex.label.dist=0.4, vertex.label.font=2)
legend_name <- paste("C",1:length(mcolors)," (n=",com$csize,", pval=",signif(com$significance,digits=2),")",sep='')
legend("topleft", legend=legend_name, fill=mcolors, bty="n", cex=0.6)

# fit a Cox proportional hazards model using the subnetwork
## for the whole network
data_g <- t(md[V(net)$name,])
data_g <- apply(data_g!=0, 1, sum)
data <- cbind(pd, net=data_g)
fit <- coxph(formula=Surv(time,status) ~ Age + Gender + TCGA_tumor_type + net, data=data)
res <- as.matrix(anova(fit))
HR_g <- res[5,2]
pvals_g <- res[5,4]
## for the cumulative nodes from the network
cg_names <- names(sort(HR[V(net)$name], decreasing=T))
cg_signif <- matrix(1, nrow=length(cg_names), ncol=2)
rownames(cg_signif) <- cg_names
colnames(cg_signif) <- c("HR", "pvalue")
for(i in 1:length(cg_names)){
    data_g <- t(md[cg_names[1:i],])
    if(i!=1){
        data_g <- apply(data_g!=0, 1, sum)
    }else{
        data_g <- as.vector(data_g!=0)
    }
    data <- cbind(pd, cnet=data_g)
    fit <- coxph(formula=Surv(time,status) ~ Age + Gender + TCGA_tumor_type + cnet, data=data)
    res <- as.matrix(anova(fit))
    cg_signif[i,] <- res[5,c(2,4)]
}
cg_signif[cg_signif[,2]==0,2] <- min(cg_signif[cg_signif[,2]!=0,2])
naive <- sample(HR, length(cg_names))
bp.HR.list <- list(All=naive, Neti=HR[cg_names], Netc=cg_signif[2:nrow(cg_signif),1])
par(las=2, mar=c(10,8,4,2)) # all axis labels horizontal
boxplot(bp.HR.list, outline=F, horizontal=F, names=c("naive\n(genes in random)", "dnet\n(genes individually)", "dnet \n(genes in combination)"), col=c("red","green","blue"), ylab="Cox hazard ratio (HR)", log="y", ylim=c(0.1,100))
# Two-sample Kolmogorov-Smirnov test
## Genes randomly choosen versus genes in the network (used individually)
stats::ks.test(x=naive, y=HR[cg_names], alternative="two.sided", exact=NULL)
## Genes in the network (used individually) versuse genes in the network (used in combination)
stats::ks.test(x=HR[cg_names], y=cg_signif[2:nrow(cg_signif),1], alternative="two.sided", exact=NULL)

# Network-based sample classifications and visualisations on 2D sample landscape
# it uses the gene-active subnetwork overlaid by mutation frequency data
frac_mutated <- sapply(tumor_type, function(x) {
    e <- eset[, which(pData(eset)$TCGA_tumor_type==x)]
    apply(exprs(e)!=0,1,sum)/ncol(e)
})
rownames(frac_mutated) <- fData(eset)$Symbol
data <- frac_mutated[V(net)$name,]
sReorder <- dNetReorder(g=net, data, feature="edge", node.normalise="degree", amplifier=3, metric="none")
visNetReorder(g=net, data=data, sReorder=sReorder, height=ceiling(sqrt(ncol(data)))*3, newpage=T, glayout=glayout, colormap="darkgreen-lightgreen-lightpink-darkred", vertex.label=NA,vertex.shape="sphere", vertex.size=16,mtext.cex=0.8,border.color="888888", zlim=c(0,0.1), mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

# output the subnetwork and their mutation frequency data
## Write the subnetwork into a SIF-formatted file (Simple Interaction File)
sif <- data.frame(source=get.edgelist(net)[,1], type="interaction", target=get.edgelist(net)[,2])
write.table(sif, file=paste("Survival_TCGA.sif", sep=""), quote=F, row.names=F,col.names=F,sep="\t")
## Output the corresponding mutation frequency data
hmap <- data.frame(Symbol=rownames(data), data)
write.table(hmap, file=paste("Survival_TCGA.txt", sep=""), quote=F, row.names=F,col.names=T,sep="\t")

# define the "mutation ubiquity" of genes in terms of a vector which stores the fraction of samples (within a tumor type) having the mutated gene 
# sparseness for a vector is: 1) one if the vector contains only a single non-zero value; 2) zero if and only if all elements are equal; 3) otherwise, the value interpolates smoothly between the two extremes
sparseness <- sapply(1:nrow(frac_mutated), function(i){
    v <- frac_mutated[i,]
    n <- length(v)
    norm1 <- sum(abs(v))
    norm2 <- sqrt(sum(v^2))
    (sqrt(n)-norm1/norm2) / (sqrt(n)-1)
})
sparseness <- matrix(sparseness, ncol=1)
rownames(sparseness) <- rownames(frac_mutated)
# derive the "mutation ubiquity" of genes: mutational fraction with the same type, and fraction consistent across different types
ubiquity <- 1- sparseness
dev.new()
hist(ubiquity,20, xlab="Cross-tumor mutation ubiquity", xlim=c(0,1))

# GSEA using the network as a gene set against the cross-tumor mutation ubiquity and each tumor type
customised.genesets <- list(net=V(net)$name)
eTerm <- dGSEA(data=cbind(ubiquity, frac_mutated), identity="symbol", genome="Hs", ontology="Customised", customised.genesets=customised.genesets, weight=0, nperm=2000)
## Comparing normalised enrichement score (NES)
frac_pvalue <- as.vector(eTerm$pvalue)
frac_fdr <- stats::p.adjust(frac_pvalue, method="BH")
frac_nes <- as.vector(eTerm$nes)
frac_es <- as.vector(eTerm$es)
df <- cbind(frac_es, frac_nes, frac_pvalue, frac_fdr)
rownames(df) <- colnames(eTerm$es)
rownames(df)[nrow(df)] <- "Mutation\nubiquity"
ind <- sort.int(frac_es, index.return=T)$ix
data <- df[ind,]
par(las=1) # make label text perpendicular to axis
par(mar=c(5,8,4,2)) # increase y-axis margin.
z <- data[,2]
barY <- barplot(z, xlab="Normalised enrichment score (NES)", horiz=TRUE, names.arg=rownames(data), cex.names=0.7, cex.lab=0.7, cex.axis=0.7, col="transparent")
## GSEA plot for ubiquity
visGSEA(eTerm, which_sample=1)

# Evolutionary analysis for genes in the subnetwork
## get a list of genes in the subnetwork
eTerm <- dEnricher(data=V(net)$name, identity="symbol", genome="Hs", ontology="PS2", sizeRange=c(10,20000), min.overlap=0)
## Look at the evolution relevance along the path to the eukaryotic common ancestor
output <- dEnricherView(eTerm, top_num=NULL, sortBy="none", details=TRUE)
write.table(output, file="enrichment_PS2.txt", quote=F, row.names=F,col.names=T,sep="\t")
output
## load Entrezgene info
org.Hs.eg <- dRDataLoader(RData='org.Hs.eg')
gene_info <- org.Hs.eg$gene_info
entrez <- unlist(eTerm$overlap[6], use.names=F)
## build neighbor-joining tree
data <- frac_mutated[V(net)$name,]
tree_bs <- visTreeBootstrap(t(data), nodelabels.arg=list(cex=0.7,bg="white-pink-violet"), metric=c("euclidean","pearson","spearman","cos","manhattan","kendall","mi")[3], num.bootstrap=2000, plot.phylo.arg=list(cex=1, edge.width=1.2))
flag <- match(tree_bs$tip.label, colnames(data))
base <- sapply(eTerm$overlap, function(x){
    as.character(gene_info[match(x,rownames(gene_info)),2])
})
## reordering via hierarchical clustering
if(1){
    cluster_order <- matrix(1, nrow=length(base))
    base_order <- matrix(1, nrow=length(base))
    for(i in 1:length(base)){
        tmp <- base[[i]]
        ind <- match(tmp, rownames(data))
        if(length(ind)>0){
            base_order[ind] <- i
            tmpD <- data[ind,]
            if(length(tmp) != 1){
                distance <- as.dist(sDistance(tmpD, metric="pearson"))
                cluster <- hclust(distance, method="average")
                cluster_order[ind] <- cluster$order
            }else if(length(tmp) == 1){
                cluster_order[ind] <- 1
            }
        }
    }
    ## contruct data frame including 1st column for temporary index, 2nd for cluster order, 3rd for base/cluster ID
    df <- data.frame(ind=1:nrow(data), cluster_order, base_order)
    # order by: first base, then hexagon
    ordering <- df[order(base_order,cluster_order),]$ind
}
RowSideColors <- sapply(1:length(base), function(x) base_order==x)
RowSideColors <- t(RowSideColors)
rslab <- ifelse(eTerm$adjp<0.05," (FDR<0.05)","")
rslab <- paste(gsub(".*:","",eTerm$set_info$name), rslab, sep="")
rownames(RowSideColors) <- rslab
colnames(RowSideColors) <- rownames(data)
RowSideColors <- ifelse(RowSideColors==T, "gray","white")
RowSideColors <- RowSideColors[, ordering]
base_order1 <- base_order[ordering]
basesep_index <- sapply(unique(base_order1), function(x) which(base_order1[length(base_order1):1]==x)[1])
basesep_index <- basesep_index[1:length(basesep_index)-1]
labRow <- sapply(pvals[match(V(net)$name, names(pvals))], function(x){
    if(x < 0.005){
        " ***"
    }else if(x < 0.01){
        " **"
    }else if(x<0.05){
        " *"
    }else{
        ""
    }
})
labRow <- paste(rownames(data), labRow, sep="")
visHeatmapAdv(data=data[ordering,flag], Rowv=F, Colv=F, colormap="lightyellow-orange", zlim=c(0,0.12), keysize=1.5, RowSideColors=RowSideColors, RowSideWidth=2, RowSideLabelLocation="top", add.expr=abline(h=(basesep_index-0.5), lty=2,lwd=1,col="black"), offsetRow=-0.5, labRow=labRow[ordering], KeyValueName="Frequency", margins=c(6,6))

# Cross-tumor mutation ubiquity versus common ancestors
ind <- match(V(net)$name, rownames(ubiquity))
net_ubiquity <- ubiquity[ind]
net_ubiquity <- net_ubiquity[ordering]
names(net_ubiquity) <- rownames(data)[ordering]
data <- cbind(net_ubiquity, base_order1)
par(las=2, mar=c(12,8,4,2)) # all axis labels horizontal
lbls <- eTerm$set_info$name[unique(base_order1)]
lbls <- gsub(".*:","",lbls)
visBoxplotAdv(formula=net_ubiquity ~ base_order1, data=data, method=c("center","hex","square","swarm")[4], pch=19, xlab="", ylab="Cross-tumor mutation ubiquity", ylim=c(0,1), labels=lbls)
## Deuterostomia versus all ancestors
stats::ks.test(x=net_ubiquity[base_order1==6], y=net_ubiquity, alternative="two.sided", exact=NULL)
## Deuterostomia versus ancestors before Deuterostomia
stats::ks.test(x=net_ubiquity[base_order1==6], y=net_ubiquity[base_order1<6], alternative="two.sided", exact=NULL)
## Deuterostomia versus ancestors after Deuterostomia
stats::ks.test(x=net_ubiquity[base_order1==6], y=net_ubiquity[base_order1>6], alternative="two.sided", exact=NULL)

# GOBP enrichment analysis
eTerm <- dEnricher(data=V(net)$name, identity="symbol", genome="Hs", ontology="GOBP")
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

# GOMF enrichment analysis
eTerm <- dEnricher(data=V(net)$name, identity="symbol", genome="Hs", ontology="GOMF")
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

# MP enrichment analysis
#eTerm <- dEnricher(data=V(net)$name, identity="symbol", genome="Hs", ontology="MP", ontology.algorithm="elim")
eTerm <- dEnricher(data=V(net)$name, identity="symbol", genome="Hs", ontology="MP", min.overlap=4)
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

# DO enrichment analysis
#eTerm <- dEnricher(data=V(net)$name, identity="symbol", genome="Hs", ontology="DO", ontology.algorithm="pc")
eTerm <- dEnricher(data=V(net)$name, identity="symbol", genome="Hs", ontology="DO")
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
names(nodes.highlight) <- nodes_query
## induce all possible shortest paths to the ontology root
subg <- dDAGinduce(g, nodes_query)
## color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
visDAG(g=subg, data=-1*log10(eTerm$adjp[V(subg)$name]), node.info="both", zlim=c(0,4), node.attrs=list(color=nodes.highlight))

# Pathway enrichment analysis
## using CP
eTerm <- dEnricher(data=V(net)$name, identity="symbol", genome="Hs", ontology="MsigdbC2CP")
dEnricherView(eTerm, top_num=10, sortBy="adjp", details=TRUE)
## using KEGG
eTerm <- dEnricher(data=V(net)$name, identity="symbol", genome="Hs", ontology="MsigdbC2KEGG")
dEnricherView(eTerm, top_num=10, sortBy="adjp", details=TRUE)
## uisng REACTOME
eTerm <- dEnricher(data=V(net)$name, identity="symbol", genome="Hs", ontology="MsigdbC2REACTOME")
dEnricherView(eTerm, top_num=10, sortBy="adjp", details=TRUE)
## uisng BIOCARTA
eTerm <- dEnricher(data=V(net)$name, identity="symbol", genome="Hs", ontology="MsigdbC2BIOCARTA", min.overlap=2)
dEnricherView(eTerm, top_num=10, sortBy="adjp", details=TRUE)

# SCOP superfamily domain enrichment analysis
eTerm <- dEnricher(data=V(net)$name, identity="symbol", genome="Hs", ontology="SF")
dEnricherView(eTerm, top_num=10, sortBy="adjp", details=TRUE)

# DGIdb druggable gene category enrichment analysis
eTerm <- dEnricher(data=V(net)$name, identity="symbol", genome="Hs", ontology="DGIdb", sizeRange=c(10,5000), min.overlap=2)
dEnricherView(eTerm, top_num=10, sortBy="adjp", details=F)

# TFBS enrichment analysis
eTerm <- dEnricher(data=V(net)$name, identity="symbol", genome="Hs", ontology="MsigdbC3TFT")
dEnricherView(eTerm, top_num=10, sortBy="adjp", details=F)

# miRNA target enrichment analysis
eTerm <- dEnricher(data=V(net)$name, identity="symbol", genome="Hs", ontology="MsigdbC3MIR")
dEnricherView(eTerm, top_num=10, sortBy="adjp", details=F)
