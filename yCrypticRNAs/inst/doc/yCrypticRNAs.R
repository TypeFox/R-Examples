## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = FALSE, comment = "", prompt = TRUE)

## ----eval = FALSE--------------------------------------------------------
#  install.packages("yCrypticRNAs")

## ------------------------------------------------------------------------
library("yCrypticRNAs")

## ---- results="hide", width = 20-----------------------------------------
#create A coverageDataSet
# Suppose an RNA-seq experiment that was done in duplicates in wild-type cells and mutant cells.
# Data were sequenced in strand-specific manner and we wish to use fragments for paired-end reads.
types <- c("wt", "wt", "mut", "mut")
bamfiles <- system.file("extdata", paste0(types, "_rep", 1:2,".bam"),
                        package = "yCrypticRNAs")
scaling_factors <- c(0.069872847, 0.081113079, 0.088520251, 0.069911116)

data("annotations")
data <- coverageDataSet(bamfiles = bamfiles, annotations = annotations,
                        types = types, sf = scaling_factors,
                        paired_end = TRUE, as_fragments = TRUE)

## ------------------------------------------------------------------------
# load introns annotations
data("introns")

# create geneCoverage dataset for FLO8 gene
flo8 <- gene_coverage(coverageDataSet = data, name = "YER109C", introns = introns)

## ---- width=60-----------------------------------------------------------
# using the 3'/5' ratio method
ratio_score(geneCoverage = flo8)

# using the 3' enrichemnt method
enrichment_score(geneCoverage = flo8)

# using the probabilistic method
zscore_score(geneCoverage = flo8)


## ---- eval=FALSE---------------------------------------------------------
#  genome_wide_scores(coverageDataSet = data, method = "ratio",
#                     outfile = "/tmp/ratio_scores.txt", introns = introns)
#  genome_wide_scores(coverageDataSet = data, method = "enrichment",
#                     outfile = "/tmp/enrichment_scores.txt", introns = introns)
#  genome_wide_scores(coverageDataSet = data, method = "probabilistic",
#                     outfile = "/tmp/zscores.txt", introns = introns)

## ---- fig.width=7, fig.height=6------------------------------------------
#calculating the cTSS for FLO8 gene
flo8_cTSS <- initiation_sites(
  name = "YER109C", bamfiles = bamfiles,
  types = types, annotations = annotations,
  introns = introns, sf = scaling_factors, replicates = 100
)

#visualize results
par(mfrow = c(3,2))
plot(flo8, cTSS = flo8_cTSS, method = "methodC_gaussian")
plot(flo8, cTSS = flo8_cTSS, method = "methodA")
plot(flo8, cTSS = flo8_cTSS, method = "methodB")
plot(flo8, cTSS = flo8_cTSS, method = "methodC")
plot(flo8, cTSS = flo8_cTSS, method = "methodD")

