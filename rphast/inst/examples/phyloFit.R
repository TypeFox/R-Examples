exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
files <- c("ENr334-100k.maf", "ENr334-100k.fa", "gencode.ENr334-100k.gff", "rev.mod")
unzip(exampleArchive, files)
m <- read.msa("ENr334-100k.maf")
mod <- phyloFit(m, tree="((hg18, (mm9, rn4)), canFam2)")
mod
phyloFit(m, init.mod=mod)
likelihood.msa(m, mod)
mod$likelihood
print(mod$likelihood, digits=10)
f <- read.feat("gencode.ENr334-100k.gff")
mod <- phyloFit(m, tree="((hg18, (mm9, rn4)), canFam2)",
                features=f, quiet=TRUE)
names(mod)
mod$other
mod[["5'flank"]]
phyloFit(m, init.mod=mod$AR, nrates=3, alpha=4.0)
phyloFit(m, init.mod=mod$AR, rate.constants=c(10, 5, 1))
# background frequencies options

# this should use the background frequencies from the initial mod
phyloFit(m, init.mod=mod$AR, quiet=TRUE)$backgd
mod$AR$backgd

# this should use the background frequencies from the data
phyloFit(m, init.mod=mod$AR, init.backgd.from.data=TRUE, quiet=TRUE)$backgd
mod$AR$backgd

# this should optimize the background frequencies
phyloFit(m, init.mod=mod$AR, no.opt=NULL, quiet=TRUE)$backgd
mod$AR$backgd

unlink(files)
