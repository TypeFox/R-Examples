files <- c("rev.mod", "ENr334-100k.maf", "ENr334-100k.fa", "small.gff")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, files)
msa <- read.msa("ENr334-100k.fa")
mod <- read.tm("rev.mod")
likelihood.msa(msa, mod)
like1 <- likelihood.msa(msa, mod, by.column=TRUE)
length(like1)==ncol.msa(msa)
sum(like1)
msa <- read.msa("ENr334-100k.maf")
likelihood.msa(msa, mod)
like2 <- likelihood.msa(msa, mod, by.column=TRUE)
sum(like2)
mod$subst.mod <- "JC69"
likelihood.msa(msa, mod)
#'
# can also get likelihood by feature
features <- read.feat("small.gff")
features$seqname <- names(msa)[1]
likelihood.msa(msa, mod, features=features)
unlink(files)
