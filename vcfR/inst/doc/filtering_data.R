## ------------------------------------------------------------------------
library(vcfR)

vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = "pinfsc50")
dna_file <- system.file("extdata", "pinf_sc50.fasta", package = "pinfsc50")
gff_file <- system.file("extdata", "pinf_sc50.gff", package = "pinfsc50")

vcf <- read.vcfR(vcf_file, verbose = FALSE)
dna <- ape::read.dna(dna_file, format = "fasta")
gff <- read.table(gff_file, sep="\t", quote="")

chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
chrom <- proc.chromR(chrom, verbose = TRUE)


## ------------------------------------------------------------------------
head(chrom)

## ---- tidy=TRUE----------------------------------------------------------
strwrap(grep("ID=MQ,", chrom@vcf@meta, value=T))
strwrap(grep("ID=DP,", chrom@vcf@meta, value=T))

## ---- fig.height=7, fig.width=7------------------------------------------
plot(chrom)

## ---- fig.height=7, fig.width=7------------------------------------------
chromoqc(chrom, dp.alpha = 22)

## ---- fig.height=7, fig.width=7------------------------------------------
chrom <- masker(chrom, min_QUAL=0, min_DP=350, max_DP=650, min_MQ=59.5, max_MQ=60.5)
chrom <- proc.chromR(chrom, verbose = FALSE)
plot(chrom)

## ---- fig.height=7, fig.width=7------------------------------------------
chromoqc(chrom, dp.alpha = 22)

## ---- eval=FALSE---------------------------------------------------------
#  write.vcf(chrom, file="good_variants.vcf.gz", mask=TRUE)

