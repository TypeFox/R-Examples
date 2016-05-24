exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, "ENr334-100k.maf")
m <- read.msa("ENr334-100k.maf")
mod <- phyloFit(m, tree="((hg18,(mm9,rn4)),canFam2)")
x <- total.expected.subs.msa(sub.msa(m, start.col=41447839, end.col=41448033, refseq="hg18"), mod)
dim(x)
dimnames(x)
x["mm9-rn4",,]
unlink("ENr334-100k.maf")
