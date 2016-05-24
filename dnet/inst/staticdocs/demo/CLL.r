# This is a demo for chronic lymphocytic leukemia patient expression data and network data
# 
###############################################################################
library(dnet)

# Load or install packages specifically used in this demo
for(pkg in c("Biobase","limma")){
    if(!require(pkg, character.only=T)){
        source("http://bioconductor.org/biocLite.R")
        biocLite(pkg)
        lapply(pkg, library, character.only=T)
    }
}

# This dataset involves 130 patients with chronic lymphocytic leukemia (CLL). When enrolled in the study, these CLL patients had not received prior therapy for CLL. Additional covariate about sampling time to first treatment (years) is available. The dataset has been normalised and log2-transformed, and provided as an "ExpressionSet" object.
CLL <- dRDataLoader(RData='CLL')
CLL

# Non-specific probesets filtering: the signals below log2(30) are considered as being technically unreliable (i.e. an empirically determined value of minimum sensitivity). Also, those probesets with >70% of samples having technically unriliable signals were excluded from further consideration
sensVal <- log2(30)
filter_flag <- apply(exprs(CLL)<sensVal,1,sum) > 0.7*ncol(exprs(CLL))
eset <- CLL[!filter_flag,]

# Create esetNew after replacing all those less than "sensVal" with "sensVal"
new.matrix <- exprs(eset)
new.matrix[new.matrix <= sensVal] <- sensVal
esetNew <- new("ExpressionSet",exprs=new.matrix,phenoData=as(pData(eset),"AnnotatedDataFrame"),featureData=as(fData(eset),"AnnotatedDataFrame"))

# A function to convert probeset-centric to entrezgene-centric expression levels
prob2gene <- function(eset){
    fdat <- fData(eset)
    tmp <- as.matrix(unique(fdat[c("EntrezID", "Symbol", "Desc")]))
    forder <- tmp[order(as.numeric(tmp[,1])),]
    forder <- forder[!is.na(forder[,1]),]
    rownames(forder) <- forder[,2]
    system.time({
        dat <- exprs(eset)
        edat <- matrix(data=NA, nrow=nrow(forder), ncol=ncol(dat))
        for (i in 1:nrow(forder)){
            ind <- which(fdat$EntrezID==forder[i,1])
            if (length(ind) == 1){
                edat[i,] <- dat[ind,]
            }else{
                edat[i,] <- apply(dat[ind,],2,mean)
            }
        }
    })
    rownames(edat) <- rownames(forder) # as gene symbols
    colnames(edat) <- rownames(pData(eset))
    esetGene <- new("ExpressionSet",exprs=data.frame(edat),phenoData=as(pData(eset),"AnnotatedDataFrame"),featureData=as(data.frame(forder),"AnnotatedDataFrame"))
    return(esetGene)
}
esetGene <- prob2gene(esetNew)
esetGene

# An igraph object that contains a functional protein association network in human. The network is extracted from the STRING database (version 10). Only those associations with medium confidence (score>=400) are retained.
org.Hs.string <- dRDataLoader(RData='org.Hs.string')
org.Hs.string

# extract network that only contains genes in esetGene
ind <- match(V(org.Hs.string)$symbol, rownames(esetGene))
## for extracted expression
esetGeneSub <- esetGene[ind[!is.na(ind)],]
esetGeneSub
## for extracted graph
nodes_mapped <- V(org.Hs.string)$name[!is.na(ind)]
network <- dNetInduce(g=org.Hs.string, nodes_query=nodes_mapped, knn=0, remove.loops=T, largest.comp=T)
V(network)$name <- V(network)$symbol
network

# according to sampling time to first treatment (years), patient samples are categorised into 3 groups: early-phase disease (E) if they were collected more than 4 years before treatment, intermediate phase (I) if collected 4 or less, but 1 or more, years before treatment (yellow bars), or late phase (L) if collected less than 1 year before treatment.
EIL <- sapply(pData(esetGeneSub)$Time, function(x) {
    if(x>4){
        "E"
    }else if(x<1){
        "L"
    }else{
        "I"
    }
})

# 1) preparation of node p-values
design <- model.matrix(~ -1 + factor(EIL))
colnames(design)<- c("E", "I", "L")
contrast.matrix <- makeContrasts(L-E, I-E, L-I, levels=design)
colnames(contrast.matrix) <- c("L_E", "I_E", "L_I")
contrast.matrix 
fit <- lmFit(exprs(esetGene), design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
# for p-value
pvals<-as.matrix(fit2$p.value)
# for adjusted p-value
adjpvals <- sapply(1:ncol(pvals),function(x) {
    p.adjust(pvals[,x], method="BH")
})
colnames(adjpvals) <- colnames(pvals)
# num of differentially expressed genes
apply(adjpvals<1e-1, 2, sum)
# only for the comparisons from late phase (L) against the early stage (E) 
my_contrast <- "L_E"
# get the p-values and calculate the scores thereupon
pval <- pvals[,my_contrast]

# 2) identification of module
g <- dNetPipeline(g=network, pval=pval, nsize=20)
glayout <- layout.fruchterman.reingold(g)

# 3) color nodes according to communities identified via a spin-glass model and simulated annealing
#com <- walktrap.community(g, modularity=T)
com <- spinglass.community(g, spins=25)
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
ind <- union(which(EIL=="L"), which(EIL=="E"))
ind <- sample(ind, 10) # sampling randomly 10 without replacement.
data <- exprs(esetGene)[V(g)$name,ind]
colnames(data) <- EIL[ind]
visNetMul(g=g, data=data, height=ceiling(sqrt(ncol(data)))*2, newpage=T,glayout=glayout,colormap=colormap,vertex.label=NA,vertex.shape="sphere", vertex.size=18,mtext.cex=0.8,border.color="888888", mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

# 11) color by additional data (be reordered)
sReorder <- dNetReorder(g, data, feature="edge", node.normalise="degree", amplifier=2, metric="none")
visNetReorder(g=g, data=data, sReorder=sReorder, height=ceiling(sqrt(ncol(data)))*2, newpage=T, glayout=glayout, colormap=colormap, vertex.label=NA,vertex.shape="sphere", vertex.size=18,mtext.cex=0.8,border.color="888888", mark.groups=mark.groups, mark.col=mark.col, mark.border=NA, mark.shape=1, mark.expand=10, edge.color=edge.color)
