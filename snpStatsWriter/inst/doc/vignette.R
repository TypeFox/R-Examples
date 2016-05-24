
## ------------------------------------------------------------------------
library(snpStatsWriter)
data(testdata,package="snpStats")
A.small <- Autosomes[1:6,1:10]
nsnps <- ncol(A.small)


## ------------------------------------------------------------------------
f <- tempfile()
write.snphap(A.small, file=f)
head(read.table(f,sep="\t"))
unlink(f)


## ------------------------------------------------------------------------
pf <- tempfile() ## pedigree file
mf <- tempfile() ## marker file
write.mach(A.small, a1=rep("1",nsnps), a2=rep("2",nsnps), pedfile=pf, mfile=mf)
head(read.table(mf))
head(read.table(pf))
unlink(pf)
unlink(mf)


## ------------------------------------------------------------------------
pf <- tempfile()
write.impute(A.small, a1=rep("1",nsnps), a2=rep("2",nsnps), bp=1:nsnps, pedfile=pf)
head(read.table(pf))
unlink(pf)


## ------------------------------------------------------------------------
gf <- tempfile() ## genotype file
mf <- tempfile() ## marker file
write.beagle(A.small, a1=rep("1",nsnps), a2=rep("2",nsnps), bp=1:nsnps, gfile=gf, mfile=mf)
head(read.table(gf,header=TRUE))
head(read.table(mf))
unlink(gf)
unlink(mf)


## ------------------------------------------------------------------------
f <- tempfile()
write.phase(A.small, file=f)
head(scan(f,what=""))
unlink(f)


