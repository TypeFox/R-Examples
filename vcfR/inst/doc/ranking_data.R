## ------------------------------------------------------------------------
library(vcfR)

vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = "pinfsc50")
dna_file <- system.file("extdata", "pinf_sc50.fasta", package = "pinfsc50")
gff_file <- system.file("extdata", "pinf_sc50.gff", package = "pinfsc50")

vcf <- read.vcfR(vcf_file, verbose = FALSE)
dna <- ape::read.dna(dna_file, format = "fasta")
gff <- read.table(gff_file, sep="\t", quote = "")

chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
#chrom <- masker(chrom, min_DP = 900, max_DP = 1500)
chrom <- proc.chromR(chrom, verbose = TRUE)


## ------------------------------------------------------------------------
head(chrom)

## ------------------------------------------------------------------------
gq <- extract.gt(chrom, element="GQ", as.numeric=TRUE)
dp <- extract.gt(chrom, element="DP", as.numeric=TRUE)

## ---- fig.align='center', fig.width=7, fig.height=5----------------------
#hist(gq[,1])
par( mar = c(8,4,4,2) )
boxplot(gq, las=2, col=2:5, main="Genotype Quality (GQ)")

dp2 <- dp
dp2[ dp2 == 0 ] <- NA
boxplot(dp2, las=2, col=2:5, main="Sequence Depth (DP)", log="y")
abline(h=10^c(0:4), lty=3, col="#808080")
par( mar = c(5,4,4,2) )

## ------------------------------------------------------------------------
mids <- apply(dp, MARGIN=2, median, na.rm=TRUE)
dp2 <- sweep(dp, MARGIN=2, mids, FUN="-")
dp2 <- abs(dp2)
dp2 <- -1 * dp2

## ---- fig.align='center', fig.width=7, fig.height=5----------------------
par( mar = c(8,4,4,2) )
boxplot(dp2, las=2, col=2:5, main="Sequence Depth (DP)")
par( mar = c(5,4,4,2) )

## ------------------------------------------------------------------------
gq2 <- gq/100
range(gq2, na.rm=TRUE)

amins <- abs(apply(dp2, MARGIN=2, min, na.rm = TRUE))
dp2 <- sweep(dp2, MARGIN=2, STATS = amins, FUN="+")
dp2 <- sweep(dp2, MARGIN=2, STATS = amins, FUN="/")
range(dp2, na.rm=TRUE)

## ---- fig.align='center'-------------------------------------------------
scores <- dp2 + gq2
scores <- rowSums(scores, na.rm = TRUE)

## ---- fig.align='center'-------------------------------------------------
hist(scores, col=4)

## ------------------------------------------------------------------------
chrom <- rank.variants.chromR(chrom, scores)
head(chrom@var.info)

## ------------------------------------------------------------------------
head(chrom@var.info[,c('POS', 'MQ', 'DP', 'window_number', 'rank')])

## ------------------------------------------------------------------------
chrom@var.info$mask[chrom@var.info$rank > 1] <- FALSE

## ---- fig.height=7, fig.width=7------------------------------------------
chromoqc(chrom, dp.alpha='66')

## ------------------------------------------------------------------------
chrom <- masker( chrom, min_QUAL=0, min_DP=350, max_DP=650, min_MQ=59.5, max_MQ=60.5 )
chrom@var.info$mask[ chrom@var.info$rank > 1 ] <- FALSE

## ---- fig.height=7, fig.width=7------------------------------------------
chromoqc(chrom, dp.alpha='66')

