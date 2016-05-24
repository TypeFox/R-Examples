exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
filename <- "rev.mod"
unzip(exampleArchive, filename)
tm <- read.tm(filename)
tm <- add.ls.mod(tm, branch="mm9", subst.mod="HKY85")
plot.lsmodel.tm(tm, 1)
tm$ls.model$backgd <- c(0.9, 0.05, 0.03, 0.02)
plot.lsmodel.tm(tm, 1)
plot.rate.matrix(tm[["rate.matrix"]],
                 eq.freq=tm[["backgd"]],
                 filled=FALSE,
                 alphabet=tm[["alphabet"]])
unlink(filename)
