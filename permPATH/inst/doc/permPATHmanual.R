## ----echo=TRUE-----------------------------------------------------------
# Generate toy phenotype and gene expression data sets
# This example consists of 40 genes grouped into 5 pathways and 100 patients
# grp is a binary trait (e.g., case vs control)
# bp is a continuous trait (e.g., blood pressure)
set.seed(1234)
library(permPATH)

n = 100
K = 40
grp = rep(1:0,each=n/2)
bp = rnorm(n)

pdat = data.frame(grp, bp)
rm(grp, bp)
expdat = matrix(rnorm(K*n),K,n)

## Assign marker names g1,...,gK to the expression data set and
## patient ids id1,...,idn to the expression and phenotype data
gnames = paste("g",1:K,sep="")
rownames(expdat) = gnames
patid = paste("id",1:n,sep="")
rownames(pdat) = patid
colnames(expdat) = patid

# Group the K genes into M pathways of sizes n1,...,nM
M = 5
p = runif(M)
p = p/sum(p)
nM = rmultinom(1, size=K, prob=p)
gset = lapply(nM, function(x){gnames[sample(x)]})
names(gset) = paste("pathway",1:M,sep="")
names(gset)

# Carry out permutation analysis with grp as the outcome
# using the two-sample Wilcoxon test with B=100 random permutations.
# The score is the maxmean test statistic
res = perm.path(expdat, y=pdat[["grp"]], local.test="wilcoxon", 
            global.test="maxmean", B=100, gset=gset, min.num=2, 
            max.num=50, sort="score")

# Output results for top pathways
head(res[["res"]])

# Output individual test statistics
res[["stats"]]

# Carry out permutation analysis with bp as the outcome
# using the Spearman test with B=100 random permutations.
# The score is the maxmean test statistic
res = perm.path(expdat, y=pdat[["bp"]], local.test="spearman", 
            global.test="maxmean", B=100, gset=gset, min.num=2, 
            max.num=50, sort="score")

# Output results for top pathways
head(res[["res"]])

# Output individual test statistics
res[["stats"]]

## ----echo=TRUE-----------------------------------------------------------
# Generate gene symbols
set.seed(1234)
library(permPATH)

gnames = c("CCL13", "CCL19", "CCL2", "CCL3", "CCL3L1", "CCL4",    
          "CCL5", "CCL7", "CCL8", "CCR1", "CCR2", "CCR3", "CCR5",
          "CD14", "CD180", "CD2", "CD209", "CD40", "CD44", "CD80",
          "CD86", "CD8A", "CDC42", "CEBPA", "CSF2", "CXCL1", "CXCL10",
          "CXCR4", "EIF2AK2", "ELK1", "ERBB2", "FCAR", "HLAA", 
          "HLADQA1", "HLADQB1", "HSPA1A", "IFIT3", "IFNA1", "IFNB1",
          "IFNG", "IL10", "IL12A", "IL12B", "IL16", "IL1A", "IL1B",
          "IL2", "IL6", "IL8", "INHBA", "IRF1", "IRF3", "ITGAM", 
          "LTA", "LYN", "MAP3K7", "MAP4K4", "MAPK8", "MAPK8IP3", 
          "MYD88", "NFKB1", "NFKBIA", "NFKBIL1", "NFRKB", "PELI1",    
          "PTGS2", "REL", "RELA", "RIPK2", "SARM1", "STK4", "TAP2",    
          "TGFB1", "TIRAP", "TLR1", "TLR10", "TLR2", "TLR3", "TLR4",
          "TLR5", "TLR6", "TLR7", "TLR8", "TLR9", "TNF", "UBE2N", "B2M",
          "RPL13A", "ACTB", "HGDC", "RTC1", "RTC2", "RTC3", "PPC1", "PPC2", "PPC3")

# extract publicly available pre-defined pathways
xx = readLines("http://software.broadinstitute.org/gsea/resources/msigdb/4.0/c2.cp.reactome.v4.0.symbols.gmt")
pnames = as.character(sapply(xx, function(x){unlist(strsplit(x, "\t", fixed=TRUE))[1]}))
anno = as.character(sapply(xx, function(x){unlist(strsplit(x, "\t", fixed=TRUE))[2]}))
gset = lapply(xx, function(x){unlist(strsplit(x, "\t", fixed=TRUE))[-c(1,2)]})
names(gset) = pnames
gset = list(gset, pnames, anno)

#intersect gene nsymbols with gene symbols from pathways
ind = unlist(lapply(gset[[1]], function(x){ifelse(length(intersect(x,gnames))>1, TRUE, FALSE)}))
gset[[1]] = gset[[1]][ind]
gset[[2]] = gset[[2]][ind]
gset[[3]] = gset[[3]][ind]
gset[[1]] = lapply(gset[[1]], function(x){intersect(x, gnames)})
names(gset[[1]]) = gset[[2]]
names(gset[[3]]) = gset[[2]]

#create gene expression data
n = 220
K = length(gnames)
expdat = matrix(abs(rnorm(K*n)), K, n)
rownames(expdat) = gnames
patid = paste("id",1:n,sep="")
colnames(expdat) = patid

grp = rep(1:0,each=n/2)
bp = abs(rnorm(n))
pdat = data.frame(grp, bp)
rm(grp, bp)

# Carry out permutation analysis with grp as the outcome
# using the two-sample Wilcoxon test with B=10000 random permutations.
# The score is the maxmean test statistic
res = perm.path(expdat, y=pdat[["grp"]], local.test="wilcoxon", 
            global.test="maxmean", B=10^4, gset=gset[[1]], min.num=2, 
            max.num=50, sort="score", anno=gset[[3]])

# Output results for top pathways
head(res[["res"]])

# Carry out permutation analysis with bp as the outcome
# using the Spearman test with B=10000 random permutations. 
# The score is the maxmean test statistic
res = perm.path(expdat, y=pdat[["grp"]], local.test="spearman", 
            global.test="maxmean", B=10^4, gset=gset[[1]], min.num=2, 
            max.num=50, sort="score", anno=gset[[3]])

# Output results for top pathways
head(res[["res"]])

## ----echo=TRUE-----------------------------------------------------------
library(permPATH)
set.seed(1234)
n = 100
K = 40
grp = rep(1:0,each=n/2)
bp = rnorm(n)

pdat = data.frame(grp, bp)
rm(grp, bp)
expdat = matrix(rnorm(K*n),K,n)

## Assign marker names g1,...,gK to the expression data set and
## patient ids id1,...,idn to the expression and phenotype data
gnames = paste("g",1:K,sep="")
rownames(expdat) = gnames
patid = paste("id",1:n,sep="")
rownames(pdat) = patid
colnames(expdat) = patid

#Group the K genes into M pathways of sizes n1,...,nM
M = 5
p = runif(M)
p = p/sum(p)
nM = rmultinom(1, size=K, prob=p)
gset = lapply(nM, function(x){gnames[sample(x)]})
names(gset) = paste("pathway",1:M,sep="")

## Carry out permutation analysis with grp as the outcome
## using the two-sample Wilcoxon with B=100 random permutations
res = perm.path(expdat, y=pdat[["grp"]], local.test="wilcoxon", global.test="maxmean", B=100, gset=gset, 
           min.num=2, max.num=50, sort="score")

# create an html file
#epermPATH2HTML(rstab, dir="/dir/", fname="tophits")

sessionInfo()

