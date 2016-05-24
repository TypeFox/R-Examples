exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
filename <- "rev.mod"
unzip(exampleArchive, filename)
tm <- read.tm(filename)

# change background frequencies to new value, adjusting rate matrix
mod.backgd.tm(tm, c(0.25, 0.25, 0.25, 0.25))

# change background frequencies so that GC content is 0.6
mod.backgd.tm(tm, gc=0.6)

unlink(filename)
