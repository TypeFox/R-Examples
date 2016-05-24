exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
file <- "coding.hmm"
unzip(exampleArchive, file)
# this is a 5-state hmm with states representing
# intergenic, intron, first, second, and third codon positions.
h <- read.hmm(file)
h
unlink(file)
