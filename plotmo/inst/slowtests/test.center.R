# test.center.R: test plotmo's center and ndiscrete args
# Stephen Milborrow, Berea Apr 2011

library(rpart.plot)
library(plotmo)
library(earth)
data(etitanic)
options(warn=1) # print warnings as they occur
if(!interactive())
    postscript(paper="letter")
et <- etitanic[, c("survived", "pclass", "sex", "age")]
et$pclassn <- as.numeric(et$pclass)
et <- et[c(30:80,330:380,630:680), ]

par(mfrow=c(3,3))
par(mar=c(3, 3.5, 3, 0.5))
par(mgp=c(1.5, .5, 0))

ndiscrete <- 0

#--- row 1

set.seed(844)
a1 <- lm(survived~pclassn+sex, data=et)
plotmo(a1, all2=T, do.par=F, degree1=NA, degree2=1, center=TRUE, clip=F,
       main="a1: survived~pclassn+sex\n(default ndiscrete)",
       pt.col=ifelse(et$survived, "black", "red"),
       pt.pch=".", pt.cex=2.5, lab=c(1,1,1))

set.seed(844)
plotmo(a1, degree1=1, all2=T, degree2=0, do.par=F, xflip=T, center=TRUE, clip=F,
       grid.levels=list(sex="f"), ndiscrete=ndiscrete,
       main="pclassn with sex=\"female\"",
       smooth.col="lightblue", smooth.lwd=2,
       pt.col=ifelse(et$survived, "black", "red"),
       pt.pch=".", pt.cex=2.5)

set.seed(844)
plotmo(a1, degree1=1, all2=T, degree2=0, do.par=F, xflip=T, center=TRUE, clip=F,
       grid.levels=list(sex="m"), ndiscrete=ndiscrete,
       main="pclassn with sex=\"male\"",
       smooth.col="lightblue", smooth.lwd=2,
       pt.col=ifelse(et$survived, "black", "red"),
       pt.pch=".", pt.cex=2.5)

#--- row 2

a2 <- lm(survived~pclassn*sex, data=et)
set.seed(844)
plotmo(a2, all2=T, do.par=F, degree2=1, degree1=0, center=TRUE, clip=F,
       main="a2: survived~pclassn*sex\n(default ndiscrete)")

set.seed(844)
plotmo(a2, degree1=1, all2=T, degree2=0, do.par=F, xflip=T, center=TRUE, clip=F,
       grid.levels=list(sex="f"), ndiscrete=ndiscrete,
       main="pclassn with sex=\"female\"",
       smooth.col="lightblue",  smooth.lwd=2,
       pt.col=ifelse(et$survived, "black", "red"),
       pt.pch=".", pt.cex=2.5)

set.seed(844)
plotmo(a2, degree1=1, all2=T, degree2=0, do.par=F, xflip=T, center=TRUE, clip=F,
       grid.levels=list(sex="m"), ndiscrete=ndiscrete,
       main="pclassn with sex=\"male\"",
       smooth.col="lightblue",  smooth.lwd=2,
       pt.col=ifelse(et$survived, "black", "red"),
       pt.pch=".", pt.cex=2.5)

#--- row 3

par(mfg=c(3,2))
a3 <- lm(survived~pclassn, data=et)
set.seed(844)
plotmo(a3, do.par=F, xflip=T, center=TRUE, clip=F, ndiscrete=ndiscrete,
       main="a3: survived~pclassn", degree1.col=1,
       smooth.col="lightblue",  smooth.lwd=2,
       pt.col=ifelse(et$survived, "black", "red"),
       pt.pch=".", pt.cex=2.5)

plot(0, 0, type="n", axes=FALSE, xlab="", ylab="")

#--- row 1

a4 <- earth(survived~pclassn+age, data=et, degree=2)

set.seed(844)
plotmo(a4, do.par=F, center=TRUE, clip=F, ylim=c(-.6,.7),
       main="earth: survived~pclassn+age\n(default ndiscrete)", degree1=0, all2=T)

set.seed(844)
plotmo(a4, do.par=F, xflip=F, all1=T, center=TRUE, clip=F, ylim=c(-.6,.7),
       main="a4, age with pclassn=1st", ndiscrete=ndiscrete,
       degree2=0, degree1=2,
       # grid.levels=list(pclassn="1st"),
       grid.levels=list(pclassn=1),
       smooth.col="lightblue",  smooth.lwd=2,
       pt.col=ifelse(et$survived, "black", "red"),
       pt.pch=".", pt.cex=2.5)

set.seed(844)
plotmo(a4, do.par=F, xflip=F, all1=T, center=TRUE, clip=F, ylim=c(-.6,.7),
       main="age with pclassn=3rd", ndiscrete=ndiscrete,
       degree2=0, degree1=2,
       grid.levels=list(pclassn=3),
       smooth.col="lightblue",  smooth.lwd=2,
       pt.col=ifelse(et$survived, "black", "red"),
       pt.pch=".", pt.cex=2.5)

#--- row 2

set.seed(844)
plotmo(a4, do.par=F, center=TRUE, clip=F, type2="im",
       main="a4 earth: survived~pclassn+age\n(default ndiscrete)", degree1=0, all2=T, yflip=T,
       pt.col=ifelse(et$survived, 1, "red"),
       image.col=gray(seq(6, 10, length=10) / 10), xflip=T,
       pt.pch=".", pt.cex=2)

set.seed(844)
plotmo(a4, do.par=F, xflip=F, all1=T, center=TRUE, clip=F,
       main="pclassn with age=10", ndiscrete=ndiscrete,
       degree2=0, degree1=1,
       grid.levels=list(age=10),
       smooth.col="lightblue",  smooth.lwd=2,
       pt.col=ifelse(et$survived, "black", "red"),
       pt.pch=".", pt.cex=2.5)

set.seed(844)
plotmo(a4, do.par=F, xflip=F, all1=T, center=TRUE, clip=F,
       main="pclassn with age=40",  ndiscrete=ndiscrete,
       degree2=0, degree1=1,
       grid.levels=list(age=40),
       smooth.col="lightblue",  smooth.lwd=2,
       pt.col=ifelse(et$survived, "black", "red"),
       pt.pch=".", pt.cex=2.5)

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}
