exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
featFile <- "sol1.gp"
unzip(exampleArchive, featFile)
f <- read.feat(featFile)
featTrack <- as.track.feat(f, "basic feature track")
f <- add.introns.feat(f)
geneTrack <- as.track.feat(f, "gene track", is.gene=TRUE)
plot.track(list(featTrack, geneTrack))
plot.track(list(featTrack, geneTrack, geneTrack, geneTrack, geneTrack),
           xlim=c(14800, 16000))
unlink(featFile)
