# load the package
library(momr)

#' all the data in the package
# data(package="momr")

#' load the raw and frequency test dataset
data("hs_3.3_metahit_sample_dat_raw")
data("hs_3.3_metahit_sample_dat_freq")


#' NORMALIZATION
#' This should be performed with the whole dataset (complete catalogue). 
#' But here is an exemple with the subset of the data for illustration purposes
data(hs_3.3_metahit_genesize)
norm.data <- normFreqRPKM(dat=hs_3.3_metahit_sample_dat_raw, cat=hs_3.3_metahit_genesize)


#' CLUSTERING OF SAMPLES
#pdf(file="visual_output.pdf") # open the pdf_visu

hc.data <- hierClust(data=hs_3.3_metahit_sample_dat_freq[,1:5], side="col")
clust.order <- hc.data$mat.hclust$order
#' order samples followin the hierarchical clustering
ordered.samples <- colnames(hs_3.3_metahit_sample_dat_freq[,1:5])[clust.order]
#' how close are the two first samples (spearman, rho)
hc.data$mat.rho[ordered.samples[1], ordered.samples[2]]
# select the samples closely related together
close.samples <- filt.hierClust(hc.data$mat.rho, hclust.method = "ward", plot = TRUE, filt = 0.37)

#' CLUSTER GENES ON THE MGS CATALOG
#' load the curated mgs data for the hs_3.3_metahit catalog
data("mgs_hs_3.3_metahit_sup500")
#' project a list of genes onto the mgs
genebag <- rownames(hs_3.3_metahit_sample_dat_freq)
mgs <- projectOntoMGS(genebag=genebag, list.mgs=mgs_hs_3.3_metahit_sup500)
#' extract the profile of a list of genes from the whole dataset
mgs.dat <- extractProfiles(mgs, hs_3.3_metahit_sample_dat_freq, silent=FALSE)
#' plot the barcodes
par(mfrow=c(5,1), mar=c(1,0,0,0))
for(i in 1:5){ 
  plotBarcode(mgs.dat[[i]])
}
#' compute the filtered vectors
mgs.mean.vect <- computeFilteredVectors(profile=mgs.dat, type="mean")


#' TEST RELATIONS
#' for the first 1000 genes
res.test <- testRelations(data=hs_3.3_metahit_sample_dat_freq[1:500,],
                          trait=c(rep(1,150),rep(2,142)),type="wilcoxon")
head(res.test)
print(paste("There are",sum(res.test$p<0.05, na.rm=TRUE),"significant genes and",
            sum(res.test$q<0.05, na.rm=TRUE), "after adjustment for multiple testing"))

res.test.mgs <- testRelations(data=mgs.mean.vect, trait=c(rep(1,150),rep(2,142)),type="wilcoxon")


#' DOWNSIZING & UPSIZING
#' downsize the matrix
data.downsized <- downsizeMatrix(data=hs_3.3_metahit_sample_dat_raw[,1:5],level=600,repetitions=1)
colSums(data.downsized, na.rm=TRUE)

#' downsize the genecount
data.genenb <- downsizeGC(data=hs_3.3_metahit_sample_dat_raw[,1:5], level=600, repetitions=3)
par(mfrow=c(1,1), mar=c(4,4,4,4))
plot(density(colMeans(data.genenb, na.rm=TRUE)), main="density of downsized gene richness")

#dev.off()
#' End of test file
