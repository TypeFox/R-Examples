### R code from vignette source 'bcsft.Rnw'

###################################################
### code chunk number 1: bcsft.Rnw:70-74
###################################################
library(qtl)
listeria.bc2s3<-read.cross(format="csv",
  file=system.file(file.path("sampledata", "listeria.csv"), package = "qtl"),
      BC.gen=2, F.gen=3)


###################################################
### code chunk number 2: bcsft.Rnw:78-80
###################################################
data(hyper)
hyper3 <- convert2bcsft(hyper, BC.gen = 3)


###################################################
### code chunk number 3: bcsft.Rnw:88-92
###################################################
listeria.f2<-read.cross(format="csv",
   file=system.file(file.path("sampledata", "listeria.csv"), package = "qtl"))
map.bc2s3 <- est.map(listeria.bc2s3)
map.f2<-est.map(listeria.f2)


###################################################
### code chunk number 4: bcsft.Rnw:97-99
###################################################
map.bc1 <- est.map(hyper)
map.bc3<-est.map(hyper3)


###################################################
### code chunk number 5: bcsft.Rnw:103-104
###################################################
plot(map.f2, map.bc2s3, label=FALSE, main="")


###################################################
### code chunk number 6: bcsft.Rnw:111-112
###################################################
plot(map.bc1, map.bc3, label=FALSE, main="")


###################################################
### code chunk number 7: bcsft.Rnw:122-130
###################################################
listeria.bc2s3<-replace.map(listeria.bc2s3, map.f2)
listeria.f2<-replace.map(listeria.f2, map.f2)

listeria.f2<-calc.genoprob(listeria.f2, step=1 )
one.f2<-scanone(listeria.f2, method="em",pheno.col=1)

listeria.bc2s3<-calc.genoprob(listeria.bc2s3, step=1 )
one.bc2s3<-scanone(listeria.bc2s3, method="em",pheno.col=1)


###################################################
### code chunk number 8: bcsft.Rnw:134-135
###################################################
plot(one.f2, one.bc2s3, col=c("red", "purple"))


###################################################
### code chunk number 9: bcsft.Rnw:141-149
###################################################
hyper3<-replace.map(hyper3, map.bc1)
hyper<-replace.map(hyper, map.bc1)

hyper<-calc.genoprob(hyper, step=1 )
one.hyp<-scanone(hyper, method="em",pheno.col=1)

hyper3<-calc.genoprob(hyper3, step=1 )
one.hyp3<-scanone(hyper3, method="em",pheno.col=1)


###################################################
### code chunk number 10: bcsft.Rnw:153-154
###################################################
plot(one.hyp, one.hyp3, col=c("red", "purple"))


