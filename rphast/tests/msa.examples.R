library(rphast)

# ncol.msa
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
ncol.msa(m)
ncol.msa(m, names.msa(m))
print(m, print.seq=TRUE)

# nrow.msa
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
nrow.msa(m)

# offset.msa
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
offset.msa(m)
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"),
         offset=500000, pointer.only=TRUE)
offset.msa(m)

# alphabet.msa
m <- msa(seqs=c("a--acgtaa", "NN-nnnTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
alphabet.msa(m)

# is.ordered.msa
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
is.ordered.msa(m)
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), is.ordered=FALSE,
         pointer.only=TRUE)
is.ordered.msa(m)

# from.pointer.msa
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"))
m
m <- from.pointer.msa(m)
m

# as.pointer.msa
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
m
m <- as.pointer.msa(m)
m


# sub.msa
print(sub.msa(m, c("human", "rat")), print.seq=TRUE)
print(sub.msa(m, c("human", "rat"), end.col=6), print.seq=TRUE)
print(sub.msa(m, c("human", "rat"), start.col=3), print.seq=TRUE)
              
m <- msa(seqs=c("ACGT---AT", "AGGTAGTAA", "AGGAAGTAG"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
print(sub.msa(m, c("human", "rat"), start.col=3, end.col=6),
      print.seq=TRUE)
print(sub.msa(m, c("mouse"), keep=FALSE, refseq="human",
              start.col=3, end.col=4),
      print.seq=TRUE)

# names.msa
m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
         names=c("human", "mouse", "rat"),
         pointer.only=TRUE)
names.msa(m)
m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"), pointer.only=TRUE)
names.msa(m)

# strip.gaps.msa
m <- msa(seqs=c("A--ACGTAT-", "AG-AGGTAA-", "AGGAGGTA--"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
print(strip.gaps.msa(m, c("human", "mouse")), print.seq=TRUE)
print(strip.gaps.msa(m, strip.mode="any.gaps"), print.seq=TRUE)
print(strip.gaps.msa(m, strip.mode="all.gaps"), print.seq=TRUE)
print(m, print.seq=TRUE)

# read.msa
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
files <- c("ENr334-100k.maf", "ENr334-100k.fa", "gencode.ENr334-100k.gff")
unzip(exampleArchive, files)

# Read a fasta file, ENr334-100k.fa
# this file represents a 4-way alignment of the first 100k
# bp of the encode region ENr334 starting from hg18 chr6
# position 41405894
idx.offset <- 41405894
m1 <- read.msa("ENr334-100k.fa", offset=idx.offset, pointer.only=TRUE)
m1

# Now read in only a subset represented in a feature file
f <- read.feat("gencode.ENr334-100k.gff")
f$seqname <- "hg18"  # need to tweak source name to match name in alignment
m1 <- read.msa("ENr334-100k.fa", features=f, offset=idx.offset, pointer.only=TRUE)

# Can also subset on certain features
do.cats <- c("CDS", "5'flank", "3'flank")
m1 <- read.msa("ENr334-100k.fa", features=f, offset=idx.offset,
               do.cats=do.cats, pointer.only=TRUE)

# Can read MAFs similarly, but don't need offset because
# MAF file is annotated with coordinates
m2 <- read.msa("ENr334-100k.maf", features=f, do.cats=do.cats, pointer.only=TRUE)
# Also, note that when features is given and the file is
# in MAF format, the first sequence is automatically
# stripped of gaps
ncol.msa(m1)
ncol.msa(m2)
ncol.msa(m1, "hg18")

unlink(files) # clean up

rm(list = ls())
gc()
