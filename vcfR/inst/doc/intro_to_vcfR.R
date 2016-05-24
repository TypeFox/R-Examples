## ------------------------------------------------------------------------
pkg <- "pinfsc50"
vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = pkg)
dna_file <- system.file("extdata", "pinf_sc50.fasta", package = pkg)
gff_file <- system.file("extdata", "pinf_sc50.gff", package = pkg)

## ----read.vcfR-----------------------------------------------------------
library(vcfR)
vcf <- read.vcfR( vcf_file, verbose = FALSE )

## ----read.dna------------------------------------------------------------
dna <- ape::read.dna(dna_file, format = "fasta")

## ----gff-----------------------------------------------------------------
gff <- read.table(gff_file, sep="\t", quote="")

## ----create.chromR-------------------------------------------------------
library(vcfR)
chrom <- create.chromR(name='Supercontig', vcf=vcf, seq=dna, ann=gff)

## ----plot chrom, fig.align='center', fig.height=7, fig.width=7-----------
plot(chrom)

## ----masker, fig.align='center', fig.height=7, fig.width=7---------------
chrom <- masker(chrom, min_QUAL = 1, min_DP = 300, max_DP = 700, min_MQ = 59.9,  max_MQ = 60.1)
plot(chrom)

## ----proc.chromR, fig.align='center', fig.height=7, fig.width=7----------
chrom <- proc.chromR(chrom, verbose=TRUE)
plot(chrom)

## ----chromoqc1, fig.align='center', fig.height=7, fig.width=7------------
chromoqc(chrom, dp.alpha=20)

## ----masker proc.chromR, fig.align='center', fig.height=7, fig.width=7----
chrom@var.info$mask <- TRUE
chrom <- proc.chromR(chrom, verbose=FALSE)
chromoqc(chrom, dp.alpha=20)
chrom <- masker(chrom, min_QUAL = 1, min_DP = 300, max_DP = 700, min_MQ = 59.9,  max_MQ = 60.1)
chrom <- proc.chromR(chrom, verbose=FALSE)

## ----chromoqc2, fig.align='center', fig.height=7, fig.width=7------------
chromoqc(chrom, xlim=c(5e+05, 6e+05))

## ----chromoqc3, fig.align='center', fig.height=7, fig.width=7------------
record <- 130
buffer=1e3
chromoqc(chrom, xlim=c( gff[record,4]-buffer, gff[record,5]+buffer ))

