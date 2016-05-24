### R code from vignette source 'citccmstdemo.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: foo
###################################################
foo <- packageDescription("citccmst")


###################################################
### code chunk number 2: centroid
###################################################
library(citccmst)
summary(citccmst)


###################################################
### code chunk number 3: affy
###################################################
load(list.files(system.file("extdata", package="citccmst"), full.names=TRUE)) # load citvalid.exp.norm object stored in /inst/extdata

citvalid.exp.annot <- data.frame(id=rownames(citvalid.exp.norm), stringsAsFactors=FALSE, row.names=rownames(citvalid.exp.norm) )


###################################################
### code chunk number 4: fig1plot
###################################################
citvalid.citccmst <- cit.assignCcmst(   data=citvalid.exp.norm,
                                        data.annot=citvalid.exp.annot,
                                        data.colId="id",
                                        data.colMap="id" ,
                                        citccmst.colMap="Probe.Set.ID",
                                        dist.method="dqda",
                                        plot=T
                )


###################################################
### code chunk number 5: affy2
###################################################
str(citvalid.citccmst)
table(citvalid.citccmst$citccmst)
table(citvalid.citccmst$citccmst.mixed)
table(citvalid.citccmst$citccmst.core)
table(citvalid.citccmst$citccmst.confidence)


###################################################
### code chunk number 6: fig1
###################################################
citvalid.citccmst <- cit.assignCcmst(   data=citvalid.exp.norm,
                                        data.annot=citvalid.exp.annot,
                                        data.colId="id",
                                        data.colMap="id" ,
                                        citccmst.colMap="Probe.Set.ID",
                                        dist.method="dqda",
                                        plot=T
                )


