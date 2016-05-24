# last modified: 3 June 2008


# This reproduces Table 2 of Wichert et al. (2004)


# Note that for time series with even sample size the numbers
# below differ from the ones published in Wichert et al (2004)
# because in the original publication the g-statistic was erroneously
# computed with the intensity at frequency \pi included. 
# This is fixed in version 1.0.4 of GeneCycle.


# Note: in the following multiple testing is performed in four ways:
#      - as in Wichert et al. using the original Benjamini-Hochberg (1995) 
#         algorithm (this yields the "conservative" estimate)
#      - using the FDR approach of Storey and Tibshirani (2003)
#         (this yields the "less conservative" estimate)
#      - using semiparametric Fdr and local fdr, as estimated by the
#         algorithm in fdrtool
#
# Other FDR controlling approaches may lead to different sets of
# significant genes.


library("GeneCycle")
library("qvalue")


# load data sets

data(caulobacter)

# download these files from http://www.strimmerlab.org/data.html 
load("spellmann-yeast.rda.gz")
load("fibroblasts.rda.gz")
load("humanhela.rda.gz")



# calculate p-values and determine number of significant genes
find.periodic.genes <- function(dataset)
{
  cat("Computing p-values ...\n")
  pval = fisher.g.test(dataset)
   
  n1 = sum(qvalue(pval, lambda=0)$qvalues < 0.05)
  n2 = sum(qvalue(pval)$qvalues < 0.05)

  fdr.out <- fdrtool(pval, statistic="pvalue", plot=FALSE)
  n3 <- sum(fdr.out$qval < 0.05) 
  n4 <- sum(fdr.out$lfdr < 0.2) 

  cat("Conservative estimate (Wichert et al.) =", n1, "\n")
  cat("Less conservative estimate =", n2, "\n")
  cat("Semiparametric Fdr < 0.05 (fdrtool) =", n3, "\n")
  cat("Semiparametric fdr < 0.2 (fdrtool) =", n4, "\n")
}


### data analysis ###

dim(cdc15)       # 24 x 4289
find.periodic.genes(cdc15)
# Conservative estimate (Wichert et al.) = 221
# Less conservative estimate = 294
# Semiparametric Fdr < 0.05 (fdrtool) = 275 
# Semiparametric fdr < 0.2 (fdrtool) = 319 

dim(cdc28)       # 17 x 1365
find.periodic.genes(cdc28)
# Conservative estimate (Wichert et al.) = 27   (!! note: typographic error in ms. !!)
# Less conservative estimate = 56
# Semiparametric Fdr < 0.05 (fdrtool) = 63 
# Semiparametric fdr < 0.2 (fdrtool) = 179 

dim(alpha)       # 18 x 4415
find.periodic.genes(alpha)
# Conservative estimate (Wichert et al.) = 170
# Less conservative estimate = 203
# Semiparametric Fdr < 0.05 (fdrtool) = 213 
# Semiparametric fdr < 0.2 (fdrtool) = 284 

dim(elution)     # 14 x 5695
find.periodic.genes(elution)
# Conservative estimate (Wichert et al.) = 72
# Less conservative estimate = 101
# Semiparametric Fdr < 0.05 (fdrtool) = 125 
# Semiparametric fdr < 0.2 (fdrtool) = 490 

dim(caulobacter)     # 11 x 1444
find.periodic.genes(caulobacter)
# Conservative estimate (Wichert et al.) = 45
# Less conservative estimate = 45
# Semiparametric Fdr < 0.05 (fdrtool) = 52 
# Semiparametric fdr < 0.2 (fdrtool) = 94 

dim(humanN2)     # 13 x 4574
find.periodic.genes(humanN2)
# Conservative estimate (Wichert et al.) = 0
# Less conservative estimate = 0
# Semiparametric Fdr < 0.05 (fdrtool) = 0 
# Semiparametric fdr < 0.2 (fdrtool) = 0 

dim(humanN3)     # 12 x 5079
find.periodic.genes(humanN3)
# Conservative estimate (Wichert et al.) = 0
# Less conservative estimate = 0
# Semiparametric Fdr < 0.05 (fdrtool) = 0 
# Semiparametric fdr < 0.2 (fdrtool) = 0 

dim(score1)      # 12 x 14728
find.periodic.genes(score1)
# Conservative estimate (Wichert et al.) = 2
# Less conservative estimate = 2
# Semiparametric Fdr < 0.05 (fdrtool) = 2 
# Semiparametric fdr < 0.2 (fdrtool) = 1 

dim(score2)      # 26 x 15472
find.periodic.genes(score2)
# Conservative estimate (Wichert et al.) = 72
# Less conservative estimate = 211
# Semiparametric Fdr < 0.05 (fdrtool) = 185 
# Semiparametric fdr < 0.2 (fdrtool) = 430 

dim(score3)      # 48 x 39724
find.periodic.genes(score3)
# Conservative estimate (Wichert et al.) = 4250
# Less conservative estimate = 4505
# Semiparametric Fdr < 0.05 (fdrtool) = 4549 
# Semiparametric fdr < 0.2 (fdrtool) = 4447 

dim(score4)      # 19 x 39192
find.periodic.genes(score4)
# Conservative estimate (Wichert et al.) = 57
# Less conservative estimate = 60
# Semiparametric Fdr < 0.05 (fdrtool) = 60 
# Semiparametric fdr < 0.2 (fdrtool) = 122 

dim(score5)      #  9 x 34890
find.periodic.genes(score5)
# Conservative estimate (Wichert et al.) = 0
# Less conservative estimate = 0
# Semiparametric Fdr < 0.05 (fdrtool) = 0 
# Semiparametric fdr < 0.2 (fdrtool) = 0 

