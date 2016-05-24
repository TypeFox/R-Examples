## ------------------------------------------------------------------------
library(vcfR)

vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = "pinfsc50")
dna_file <- system.file("extdata", "pinf_sc50.fasta", package = "pinfsc50")
gff_file <- system.file("extdata", "pinf_sc50.gff", package = "pinfsc50")

vcf <- read.vcfR(vcf_file, verbose = FALSE)
dna <- ape::read.dna(dna_file, format = "fasta")
gff <- read.table(gff_file, sep="\t", quote="")

chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
chrom <- masker(chrom, min_DP = 300, max_DP = 700)
chrom <- proc.chromR(chrom, verbose = FALSE)


## ------------------------------------------------------------------------
head(chrom)

## ------------------------------------------------------------------------
dp <- extract.gt(chrom, element="DP", as.numeric=TRUE)
rownames(dp) <- 1:nrow(dp)
head(dp)

## ---- fig.height=7, fig.width=7------------------------------------------
heatmap.bp(dp[1001:1500,])

## ------------------------------------------------------------------------
is.na(dp[na.omit(dp == 0)]) <- TRUE

## ---- fig.height=7, fig.width=7------------------------------------------
heatmap.bp(dp[1001:1500,])

## ---- fig.height=4, fig.width=7------------------------------------------
par(mar=c(8,4,4,2))
barplot(apply(dp, MARGIN=2, mean, na.rm=TRUE), las=3)
par(mar=c(5,4,4,2))

