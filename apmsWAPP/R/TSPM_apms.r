######  Analysis of AP-MS Data: Detection of protein interaction partners and separation from contaminants
#####   on the base of spectral count Data
###### Method: Application of the TSPM model + pre- and postprocessing steps

#Steps:   - Normalization
#         - Filtering
#         - TSPM model
#         - p-value Adjustment for multiple testing by:
#             * Benjamini-Hochberg procedure: FDR controlled p-values
#             * permutation approach + Westfall&Young algorithm: FWER controlled p-values for each protein interaction-candidate

# Assumptions:
#     - Input of only one bait protein; controls required; - replicates for each group required
#     - minimum number of replicates for bait and control experiment preferred to be 3 

# Output:
#     - TSPM output (id, LFC, rawp, BH.padj, LRT dispersion)
#     - case WY-adjustment: WY adj. pvalues + counter 
#     - (filtered) (normalized) count matrix


tspm_apms  <- function(counts, baittab, norm = c("none", "sumtotal", "upperquartile", "DESeq", "TMM", "quantile"),Filter=TRUE, filter.method= c("IQR", "overallVar", "noVar"), var.cutoff=NA, limit=0, adj.method=c("BH","WY")) {
#require(multtest)
#require(gtools)

norm <- match.arg(norm)
adj.method <- match.arg(adj.method)

baittab <- read.table(baittab) 
baittab$V3 <- as.character(baittab$V3)
baittab$V2 <- as.character(baittab$V2)
baittab$V1 <- as.character(baittab$V1)
if(all(baittab$V3%in%c("C","T"))==FALSE) stop("Sample and Controls need to be defined by T respectively C")
if(setequal(baittab$V1, colnames(counts))==FALSE) stop("Names of samples need to be consistent in Baittable and Count-matrix")

baittab <- baittab[order(baittab$V3),]                    # baittab sorted: controls first, followed by bait samples
counts <- counts[ ,match(baittab$V1, colnames(counts))]   # counts in the same order as baittable

counts.org <- counts

                        
                              ####  NORMALIZATION  ####
if (norm!="none")  {
  if(norm=="quantile")  {
  norm.out <- norm.inttable (counts, baittab, norm)
  counts <- round(norm.out[[1]])
  lib.size <- rep(1, times=dim(counts)[2])      }
  else  {
  norm.out <- norm.inttable (counts, baittab, norm)
  lib.size <- norm.out[[2]]
  counts <- norm.out[[1]]          }
}
else { lib.size <- rep(1, times=dim(counts)[2]) }


                             ####    FILTERING   ####
if(Filter==TRUE) {
filter.method <- match.arg(filter.method)
counts.nf <- varFilter (counts, baittab, func=filter.method, var.cutoff, limit)  
pos.f <- match(rownames(counts.nf), rownames(counts))
if (norm!="quantile") { counts.f <- counts.org[pos.f,]   }    # filtered unnormalized counts
else   {counts.f <- counts[pos.f,] }                          # case quantile-norm
}

                            #####  TSPM model  #####
# TSPM parameter:
x1 <- factor(baittab$V3[match(colnames(counts),baittab$V1)], levels=c("C","T")) # important: order C -> T in levels
x0 <- rep(1, times=dim(counts)[2])                            
                                                             
# TSPM Modell for diff.expr.
if (Filter==TRUE) new.counts <- counts.f
else {if(norm=="quantile") new.counts <- counts else new.counts <- counts.org }

tspm <- TSPM (new.counts, x1, x0, lib.size, alpha.wh=0.05)

if(adj.method=="BH") {
if (Filter==TRUE) return(list(tspm, counts.nf)) else return(list(tspm, counts))    # output1: filtered (normalized) counts, output2: (normalized) counts
}


if (adj.method=="WY")   {
org.LRT <- tspm$LRT
classlabel <- ifelse(x1=="C",0,1)
perms <- mt.sample.label(classlabel,test="t",fixed.seed.sampling="n",B=0)
perms2 <- apply(perms, 2, function(x){ifelse(x==0,"C","T")} )
# permutation table:
proteins <- rownames(new.counts)
permmat <- as.data.frame(matrix(data=NA, nrow=length(proteins), ncol=(dim(perms)-1) ) )
rownames(permmat) <- proteins

for ( j in 2:nrow(perms2)) {
x1 <- factor(perms2[j,], levels=c("C","T"))
tspm.perm <- TSPM (new.counts, x1, x0, lib.size, alpha.wh=0.05)    # permuted TSPM model
permmat[,(j-1)] <- tspm.perm$LRT
}

set.seed(12345)                                   # Westfall&Young Algorithm
pp <- permute(c(1:length(proteins)))
WY.padj <- WY.permalg(org.LRT[pp], permmat[pp,])

if (Filter==TRUE) return(list(tspm, WY.padj, counts.nf, permmat)) else return(list(tspm, WY.padj, counts, permmat))     # output1: filtered (normalized) counts, output2: (normalized) counts
}        # wy end


}        # end
