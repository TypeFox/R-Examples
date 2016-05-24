### R code from vignette source 'jetset.Rnw'

###################################################
### code chunk number 1: jetset.Rnw:39-40
###################################################
  library(jetset)


###################################################
### code chunk number 2: jetset.Rnw:49-50
###################################################
  ls("package:jetset")


###################################################
### code chunk number 3: jetset.Rnw:71-72
###################################################
  head(scores.hgu95av2)


###################################################
### code chunk number 4: jetset.Rnw:109-113
###################################################
  jmap('hgu95av2', eg = "2099")
  jmap('hgu133a', eg = "2099")
  jmap('hgu133plus2', eg = "2099")
  jmap('u133x3p', eg = "2099")


###################################################
### code chunk number 5: jetset.Rnw:123-125
###################################################
  jmap('hgu95av2', ensembl = "ENSG00000091831")
  jmap('hgu133a', ensembl = "ENSG00000091831")


###################################################
### code chunk number 6: jetset.Rnw:133-134
###################################################
  jmap('hgu133a', symbol = c("ESR1", "ERBB2", "AURKA"))


###################################################
### code chunk number 7: jetset.Rnw:144-145
###################################################
  jmap('u133x3p', alias = c("P53", "HER-2", "K-RAS"))


###################################################
### code chunk number 8: jetset.Rnw:154-155
###################################################
  jscores('hgu95av2', symbol = 'STAT1')


###################################################
### code chunk number 9: jetset.Rnw:162-163
###################################################
  jmap('hgu95av2', symbol = 'STAT1')


###################################################
### code chunk number 10: jetset.Rnw:168-170
###################################################
  allscores <- jscores('hgu95av2')
  str(allscores)


###################################################
### code chunk number 11: jetset.Rnw:179-183
###################################################
table(!is.na(scores.hgu95av2$EntrezID))
table(!is.na(scores.hgu133a$EntrezID))
table(!is.na(scores.hgu133plus2$EntrezID))
table(!is.na(scores.u133x3p$EntrezID))


###################################################
### code chunk number 12: jetset.Rnw:188-192
###################################################
length(na.omit(unique(scores.hgu95av2$EntrezID)))
length(na.omit(unique(scores.hgu133a$EntrezID)))
length(na.omit(unique(scores.hgu133plus2$EntrezID)))
length(na.omit(unique(scores.u133x3p$EntrezID)))


###################################################
### code chunk number 13: sessionInfo
###################################################
sessionInfo()


