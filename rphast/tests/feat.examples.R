library(rphast)
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
featFile <- "gencode.ENr334.gp"

# read.feat
unzip(exampleArchive, featFile)
f <- read.feat(featFile, pointer.only=TRUE)
dim.feat(f)
unlink(featFile)

# as.pointer.feat
f <- as.pointer.feat(f)
summary(f)

# write.feat
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
f <- feat(seq, src, feature, start, end)
f <- feat(seq, src, feature, start, end)
write.feat(f, "test.gff")
unlink("test.gff")

# nrow.feat
nrow.feat(as.pointer.feat(f))

rm(list = ls())
gc()
