exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
files <- c("ENr334-100k.maf", "ENr334-100k.fa", "gencode.ENr334-100k.gff")
unzip(exampleArchive, files)

# Read a fasta file, ENr334-100k.fa
# this file represents a 4-way alignment of the encode region
# ENr334 starting from hg18 chr6 position 41405894
idx.offset <- 41405894
m1 <- read.msa("ENr334-100k.fa", offset=idx.offset)
m1

# Now read in only a subset represented in a feature file
f <- read.feat("gencode.ENr334-100k.gff")
f$seqname <- "hg18"  # need to tweak source name to match name in alignment
m1 <- read.msa("ENr334-100k.fa", features=f, offset=idx.offset)

# Can also subset on certain features
do.cats <- c("CDS", "5'flank", "3'flank")
m1 <- read.msa("ENr334-100k.fa", features=f, offset=idx.offset,
               do.cats=do.cats)

# Can read MAFs similarly, but don't need offset because
# MAF file is annotated with coordinates
m2 <- read.msa("ENr334-100k.maf", features=f, do.cats=do.cats)
# Also, note that when features is given and the file is
# in MAF format, the first sequence is automatically
# stripped of gaps
ncol.msa(m1)
ncol.msa(m2)
ncol.msa(m1, "hg18")

unlink(files) # clean up
