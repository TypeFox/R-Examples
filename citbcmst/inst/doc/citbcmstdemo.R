### R code from vignette source 'citbcmstdemo.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: foo
###################################################
foo <- packageDescription("citbcmst")


###################################################
### code chunk number 2: centroid
###################################################
library(citbcmst)
data(citbcmst)
summary(citbcmst)


###################################################
### code chunk number 3: affy
###################################################
load(list.files(system.file("extdata", package="citbcmst"), full.names=TRUE)[1])# load exp.norm.bertheau07 object stored in /inst/extdata
exp.annot.bertheau07 <- data.frame(id=rownames(exp.norm.bertheau07), stringsAsFactors=F, row.names=rownames(exp.norm.bertheau07) )


###################################################
### code chunk number 4: fig1plot
###################################################
citbcmst.bertheau07 <- cit.assignBcmst(   data=exp.norm.bertheau07,
                                          data.annot=exp.annot.bertheau07,
                                          data.colId="id",
                                          data.colMap="id" ,
                                          citbcmst.colMap="Probe.Set.ID",
                                          dist.method="dlda",
                                          plot=TRUE
                                      )


###################################################
### code chunk number 5: affy2
###################################################
str(citbcmst.bertheau07)
table(citbcmst.bertheau07$citbcmst)
table(citbcmst.bertheau07$citbcmst.mixte)
table(citbcmst.bertheau07$citbcmst.core)
table(citbcmst.bertheau07$citbcmst.confidence)


###################################################
### code chunk number 6: fig1
###################################################
citbcmst.bertheau07 <- cit.assignBcmst(   data=exp.norm.bertheau07,
                                          data.annot=exp.annot.bertheau07,
                                          data.colId="id",
                                          data.colMap="id" ,
                                          citbcmst.colMap="Probe.Set.ID",
                                          dist.method="dlda",
                                          plot=TRUE
                                      )


###################################################
### code chunk number 7: oligo
###################################################
load(list.files(system.file("extdata", package="citbcmst"), full.names=TRUE)[2]) # load exp.norm.chanrion08 object stored in /inst/extdata
exp.annot.chanrion08 <- data.frame(id=rownames(exp.norm.chanrion08), gs=rownames(exp.norm.chanrion08), stringsAsFactors=F, row.names=rownames(exp.norm.chanrion08) )


###################################################
### code chunk number 8: fig2plot
###################################################
citbcmst.chanrion08 <- cit.assignBcmst(   data=exp.norm.chanrion08,
                                          data.annot=exp.annot.chanrion08,
                                          data.colId="id",
                                          data.colMap="gs" ,
                                          citbcmst.colMap="Gene.Symbol",
                                          dist.method="pearson",
                                          plot=TRUE
                                      )


###################################################
### code chunk number 9: oligo2
###################################################
str(citbcmst.chanrion08)
table(citbcmst.chanrion08$citbcmst)
table(citbcmst.chanrion08$citbcmst.mixte)
table(citbcmst.chanrion08$citbcmst.core)
table(citbcmst.chanrion08$citbcmst.confidence)


###################################################
### code chunk number 10: fig2
###################################################
citbcmst.chanrion08 <- cit.assignBcmst(   data=exp.norm.chanrion08,
                                          data.annot=exp.annot.chanrion08,
                                          data.colId="id",
                                          data.colMap="gs" ,
                                          citbcmst.colMap="Gene.Symbol",
                                          dist.method="pearson",
                                          plot=TRUE
                                      )


