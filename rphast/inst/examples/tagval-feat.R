exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
featFile <- "sol1.gp"
unzip(exampleArchive, featFile)
f <- read.feat(featFile)
geneName <- tagval.feat(f, "transcript_id")
geneName[1:10]
length(unique(geneName)) # number of unique genes
unlink(featFile)
rm(f, geneName)
