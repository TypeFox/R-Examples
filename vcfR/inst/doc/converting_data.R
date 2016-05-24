## ------------------------------------------------------------------------
library(vcfR)
data(vcfR_example)

## ----write.vcf, eval=FALSE-----------------------------------------------
#  write.vcf(vcf, "test.vcf.gz")
#  unlink("test.vcf.gz") # Clean up after our example is done.

## ----genind, eval=TRUE---------------------------------------------------
my_genind <- vcfR2genind(vcf)
class(my_genind)
my_genind

## ----genclone, eval=TRUE-------------------------------------------------
my_genclone <- poppr::as.genclone(my_genind)
class(my_genclone)
my_genclone

## ----genlight, eval=FALSE------------------------------------------------
#  vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = "pinfsc50")
#  vcf <- read.vcfR(vcf_file, verbose = FALSE)
#  x <- vcfR2genlight(vcf)
#  x

## ----load vcf dna gff----------------------------------------------------
# Find the files.
vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = "pinfsc50")
dna_file <- system.file("extdata", "pinf_sc50.fasta", package = "pinfsc50")
gff_file <- system.file("extdata", "pinf_sc50.gff", package = "pinfsc50")
# Read in data.
vcf <- read.vcfR(vcf_file, verbose = FALSE)
dna <- ape::read.dna(dna_file, format="fasta")
gff <- read.table(gff_file, sep="\t", quote = "")

## ----vcfR2DNAbin, tidy=TRUE----------------------------------------------
record <- 130
my_dnabin1 <- vcfR2DNAbin(vcf, consensus = TRUE, extract.haps = FALSE, gt.split="|", ref.seq=dna[,gff[record,4]:gff[record,5]], start.pos=gff[record,4], verbose=FALSE)
my_dnabin1

## ----image_DNAbin1, fig.align='center', fig.width=7, fig.height=7--------
par(mar=c(5,8,4,2))
ape::image.DNAbin(my_dnabin1[,ape::seg.sites(my_dnabin1)])
par(mar=c(5,4,4,2))

## ----vcfR2DNAbin_2, tidy=TRUE--------------------------------------------
my_dnabin1 <- vcfR2DNAbin(vcf, consensus=FALSE, extract.haps=TRUE, gt.split="|", ref.seq=dna[,gff[record,4]:gff[record,5]], start.pos=gff[record,4], verbose=FALSE)

## ----image_DNAbin_2, fig.align='center', fig.width=7, fig.height=7-------
par(mar=c(5,8,4,2))
ape::image.DNAbin(my_dnabin1[,ape::seg.sites(my_dnabin1)])
par(mar=c(5,4,4,2))

## ---- eval=FALSE---------------------------------------------------------
#  write.dna( my_dnabin1, file = 'my_gene.fasta', format = 'fasta' )
#  unlink('my_gene.fasta') # Clean up after we're done with the example.

## ----vcfR2loci, eval=FALSE-----------------------------------------------
#  system.time( my_loci <- vcfR2loci(vcf) )
#  class(my_loci)

