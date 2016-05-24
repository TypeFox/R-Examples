exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
files <- c("ENr334-100k.maf", "ENr334-100k.fa", "ENr334-100k.ss", "sol1.maf", "sol1.gp")
unzip(exampleArchive, files=files)
guess.format.msa(files)
# the last file is not an alignment, which is why it returns UNKNOWN
unlink(files)
