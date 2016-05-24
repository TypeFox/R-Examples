## ---- eval=FALSE---------------------------------------------------------
#  vignette("intro_to_vcfR", package="vcfR")

## ------------------------------------------------------------------------
library(vcfR)

# Find the files.
vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = "pinfsc50")
dna_file <- system.file("extdata", "pinf_sc50.fasta", package = "pinfsc50")
gff_file <- system.file("extdata", "pinf_sc50.gff", package = "pinfsc50")

# Input the files.
vcf <- read.vcfR(vcf_file, verbose = FALSE)
dna <- ape::read.dna(dna_file, format = "fasta")
gff <- read.table(gff_file, sep="\t", quote="")

# Create a chromR object.
chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=TRUE)
chrom <- masker(chrom, min_QUAL=0, min_DP=350, max_DP=650, min_MQ=59.5, max_MQ=60.5)
chrom <- proc.chromR(chrom, verbose = TRUE)

## ---- fig.height=7, fig.width=7------------------------------------------
chromoqc(chrom, dp.alpha = 22)

## ---- fig.height=7, fig.width=7------------------------------------------
chrom <- proc.chromR(chrom, verbose=FALSE, win.size=1e4)
chromoqc(chrom, dp.alpha = 22)

## ------------------------------------------------------------------------
chrom <- proc.chromR(chrom, verbose=FALSE, win.size=1e3)

## ------------------------------------------------------------------------
head(chrom)

## ---- fig.height=7, fig.width=7------------------------------------------
plot(chrom)

