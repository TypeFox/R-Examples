### Function for Normalization of Counts in the interactiontable:
# Input: - original interaction matrix   (interaction table including zeros, raw discrete counts!)
#        - baittable = classification of samples and controls
#        - different normalization methods can be chosen in "norm":
#         "sumtotal", "upperquartile", "DESeq", "TMM", "quantile"
# Output: normed interaction-matrix (normed count table), scaling factors

norm.inttable <- function( inttab.mat, baittab, norm = c("sumtotal", "upperquartile", "DESeq", "TMM", "quantile")) {

#require(limma)
#require(edgeR)
#require(DESeq)
#require(aroma.light)
norm <- match.arg(norm)

baittab$V3 <- as.character(baittab$V3)
baittab$V1 <- as.character(baittab$V1)

Cpos <- match(baittab$V1[grep("C",baittab$V3)], colnames(inttab.mat))     # columns in inttab.mat for "C"controls
Tpos <- match(baittab$V1[grep("T",baittab$V3)], colnames(inttab.mat))     # columns in inttab.mat for "T"samples


                                  ##########    Normalization     ############
                                  
# "sumtotal" :
if(norm=="sumtotal") {
sumcount <- apply(inttab.mat, 2, function(x){sum(x)})    # sum total in samples
scal.fac <- sumcount
scal.fac[Cpos] <- scal.fac[Cpos] / median(sumcount[Cpos])
scal.fac[Tpos] <- scal.fac[Tpos] / median(sumcount[Tpos])
inttab.norm <- apply(inttab.mat, 2, function(x){x/sum(x)} )

inttab.norm[ ,Cpos] <- inttab.norm[ ,Cpos] * median(sumcount[Cpos])   # within Ctrl replicates
inttab.norm[ ,Tpos] <- inttab.norm[ ,Tpos] * median(sumcount[Tpos])   # within Bait replicates
}


# "upper-quartile" (proposed in Bullard et al. 2010):
if(norm=="upperquartile") {
upper.quartiles <- apply(inttab.mat, 2, function(x){quantile(x, probs=0.75)})
zerocount <- which(upper.quartiles==0)
scal.fac <- apply(inttab.mat, 2, function(x){qu<-quantile(x,probs=0.75); if(qu!=0)qu else 1})
scal.fac[setdiff(Cpos,zerocount)] <- scal.fac[setdiff(Cpos,zerocount)] / median(scal.fac[Cpos])
scal.fac[setdiff(Tpos,zerocount)] <- scal.fac[setdiff(Tpos,zerocount)] / median(scal.fac[Tpos])         # scaling factors

inttab.norm <- apply(inttab.mat, 2, function(x){qu<-quantile(x,probs=0.75); if(qu!=0)x/qu else x})             # normed by 75% quantile of the samples
scaling <-  ifelse(upper.quartiles==0,1,upper.quartiles)
inttab.norm[,setdiff(Cpos,zerocount)] <- inttab.norm[,setdiff(Cpos,zerocount)] * median(scaling[Cpos]) # transfer back to the count scales
inttab.norm[,setdiff(Tpos,zerocount)] <- inttab.norm[,setdiff(Tpos,zerocount)] * median(scaling[Tpos]) # normalized count matrix
}                                                                              


# "DESeq" - normalization method in the DESeq package (Anders et al. 2010):
if(norm=="DESeq") {
cds <- newCountDataSet( inttab.mat, baittab$V3)       # build countDataSet

cds1 <- estimateSizeFactors( cds[ ,Cpos] )           
cds2 <- estimateSizeFactors( cds[ ,Tpos] )
sizeFactors( cds ) [Cpos] <- sizeFactors( cds1 )
sizeFactors( cds ) [Tpos] <- sizeFactors( cds2 )
scal.fac <- sizeFactors( cds )                                               # scaling factors

inttab.norm <- scale(counts(cds), center=FALSE, scale=sizeFactors(cds))
colnames(inttab.norm) <-colnames(counts(cds))                                # normalized count matrix
}


# "TMM" - normalization method in the edgeR package (Robinson et al. 2010):
if(norm=="TMM") {
d1 <- DGEList(inttab.mat[, Cpos])
d1 <- calcNormFactors(d1, method="TMM")
eff.libsize1 <- d1$samples$lib.size * d1$samples$norm.factors

d2 <- DGEList(inttab.mat[, Tpos])
d2 <- calcNormFactors(d2, method="TMM")
eff.libsize2 <- d2$samples$lib.size * d2$samples$norm.factors

eff.libsize <- rep(NA, times=dim(inttab.mat)[2])
eff.libsize[Cpos] <- eff.libsize1
eff.libsize[Tpos] <- eff.libsize2
scal.fac <- eff.libsize                                               # scaling factors
scal.fac[Cpos] <- scal.fac[Cpos] / median(eff.libsize1)
scal.fac[Tpos] <- scal.fac[Tpos] / median(eff.libsize2)

inttab.norm <- scale(inttab.mat, center=FALSE, scale=eff.libsize)
inttab.norm[ ,Cpos] <- inttab.norm[ ,Cpos] * median(eff.libsize1)   
inttab.norm[ ,Tpos] <- inttab.norm[ ,Tpos] * median(eff.libsize2)     # normalized count matrix
}


# "quantile" normalization:
if(norm=="quantile") {
inttab.norm <- inttab.mat
inttab.norm[ ,Cpos] <- normalizeQuantileRank (inttab.mat[ ,Cpos], robust=TRUE)  # quantile normalization within replicates
inttab.norm[ ,Tpos] <- normalizeQuantileRank (inttab.mat[ ,Tpos], robust=TRUE)
scal.fac <- NA                                                                  # no scaling factors available for this method
}


return(list(inttab.norm, scal.fac))         # output

}               # function end
