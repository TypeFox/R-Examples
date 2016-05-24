## ---- fig.cap="Cartoon representation of VCF file organization", echo=FALSE, fig.height=4, fig.width=4, fig.align='center'----
par(mar=c(0.1,0.1,0.1,0.1))
plot(c(0,5), c(0,5), type="n", frame.plot=FALSE, axes=FALSE, xlab="", ylab="")
rect(xleft=0, ybottom=4, xright=3, ytop=5)
rect(xleft=0, ybottom=0, xright=2, ytop=4)
rect(xleft=2, ybottom=0, xright=5, ytop=4)
text(1.5, 4.5, "meta")
text(1.0, 2.5, "fix")
text(3.5, 2.5, "gt")
par(mar=c(5,4,4,2))

## ------------------------------------------------------------------------
library(vcfR)
data(vcfR_example)

## ---- echo=FALSE---------------------------------------------------------
strwrap(vcf@meta[1:7])

## ---- echo=FALSE---------------------------------------------------------
vcf@fix[1:6, 1:7]

## ------------------------------------------------------------------------
strwrap(grep('DP', vcf@meta, value = TRUE))

## ---- echo=TRUE----------------------------------------------------------
unlist(strsplit(as.character(vcf@fix[1, 8]), split=";"))

## ---- echo=FALSE---------------------------------------------------------
vcf@gt[1:6, 1:8]

## ------------------------------------------------------------------------
head(vcf)

