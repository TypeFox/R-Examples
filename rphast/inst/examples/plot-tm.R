exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
filename <- "rev.mod"
unzip(exampleArchive, filename)
tm <- read.tm(filename)
plot(tm)
plot(tm, show.eq.freq=FALSE)
plot(tm, max.cex=20, eq.freq.max.cex=1,
     col=matrix(1:16, nrow=4),
     eq.freq.col=c("red", "green"),
     filled=TRUE, add=TRUE)
plot.rate.matrix(tm[["rate.matrix"]],
                 eq.freq=tm[["backgd"]],
                 filled=FALSE)
plot.rate.matrix(tm[["rate.matrix"]],
                 eq.freq=tm[["backgd"]],
                 filled=TRUE, add=TRUE)
unlink(filename)
