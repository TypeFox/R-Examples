exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
featFile <- "gencode.ENr334-100k.gff"
unzip(exampleArchive, featFile)
f <- read.feat(featFile)
dim(f)
f[1:10,]
unlink(featFile)
