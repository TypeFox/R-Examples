exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
featFile <- "sol1.gp"
unzip(exampleArchive, featFile)
f <- read.feat(featFile)
plot.gene(f)
plot.gene(f, xlim=c(0, 10000))  #zoom in

unlink(featFile)
