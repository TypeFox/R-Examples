### R code from vignette source 'DEU_PISA2012_description.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: DEU_PISA2012_description.Rnw:25-26
###################################################
library(pairwise)


###################################################
### code chunk number 2: DEU_PISA2012_description.Rnw:31-32
###################################################
data(DEU_PISA2012)


###################################################
### code chunk number 3: DEU_PISA2012_description.Rnw:41-42
###################################################
names(DEU_PISA2012)


###################################################
### code chunk number 4: DEU_PISA2012_description.Rnw:47-48
###################################################
lapply(DEU_PISA2012,class)


###################################################
### code chunk number 5: DEU_PISA2012_description.Rnw:54-55
###################################################
names(DEU_PISA2012$id)


###################################################
### code chunk number 6: DEU_PISA2012_description.Rnw:62-63
###################################################
table(DEU_PISA2012$id$BOOKID)


###################################################
### code chunk number 7: DEU_PISA2012_description.Rnw:67-68
###################################################
table(DEU_PISA2012$id$QuestID)


###################################################
### code chunk number 8: DEU_PISA2012_description.Rnw:75-76
###################################################
names(DEU_PISA2012$covariate)


###################################################
### code chunk number 9: DEU_PISA2012_description.Rnw:91-92
###################################################
names(DEU_PISA2012$cog)


###################################################
### code chunk number 10: DEU_PISA2012_description.Rnw:98-99
###################################################
names(DEU_PISA2012$cog$pv)


###################################################
### code chunk number 11: DEU_PISA2012_description.Rnw:104-105
###################################################
names(DEU_PISA2012$cog$pv$MATH)


###################################################
### code chunk number 12: DEU_PISA2012_description.Rnw:108-109
###################################################
names(DEU_PISA2012$cog$pv$READ)


###################################################
### code chunk number 13: DEU_PISA2012_description.Rnw:112-113
###################################################
names(DEU_PISA2012$cog$pv$SCIE)


###################################################
### code chunk number 14: DEU_PISA2012_description.Rnw:117-118
###################################################
length(DEU_PISA2012$cog$pv$MATH$PV1MATH)


###################################################
### code chunk number 15: DEU_PISA2012_description.Rnw:123-124
###################################################
names(DEU_PISA2012$cog$dat)


###################################################
### code chunk number 16: DEU_PISA2012_description.Rnw:129-130
###################################################
rapply(DEU_PISA2012$cog$dat,names,classes = "list",how="list")


###################################################
### code chunk number 17: DEU_PISA2012_description.Rnw:137-138
###################################################
names(DEU_PISA2012$cog$dat$MATH)


###################################################
### code chunk number 18: DEU_PISA2012_description.Rnw:143-146
###################################################
dim(DEU_PISA2012$cog$dat$MATH$resp)
dim(DEU_PISA2012$cog$dat$MATH$inc7)
dim(DEU_PISA2012$cog$dat$MATH$inc8)


###################################################
### code chunk number 19: DEU_PISA2012_description.Rnw:161-162
###################################################
names(DEU_PISA2012$ncog)


###################################################
### code chunk number 20: DEU_PISA2012_description.Rnw:168-169
###################################################
names(DEU_PISA2012$ncog$wle)


###################################################
### code chunk number 21: DEU_PISA2012_description.Rnw:175-176
###################################################
names(DEU_PISA2012$ncog$dat)


###################################################
### code chunk number 22: DEU_PISA2012_description.Rnw:183-184
###################################################
names(DEU_PISA2012$ncog$dat$CLSMAN)


###################################################
### code chunk number 23: DEU_PISA2012_description.Rnw:189-190
###################################################
lapply(DEU_PISA2012$ncog$dat$CLSMAN,dim)


###################################################
### code chunk number 24: DEU_PISA2012_description.Rnw:195-196
###################################################
colnames(DEU_PISA2012$ncog$dat$CLSMAN$resp)


###################################################
### code chunk number 25: DEU_PISA2012_description.Rnw:207-208
###################################################
names(DEU_PISA2012$weights)


###################################################
### code chunk number 26: DEU_PISA2012_description.Rnw:217-218
###################################################
(sum(DEU_PISA2012$cog$dat$MATH$inc7))/(prod(dim(DEU_PISA2012$cog$dat$MATH$inc7)))*100


###################################################
### code chunk number 27: DEU_PISA2012_description.Rnw:223-224
###################################################
(sum(DEU_PISA2012$cog$dat$READ$inc7))/(prod(dim(DEU_PISA2012$cog$dat$READ$inc7)))*100


###################################################
### code chunk number 28: DEU_PISA2012_description.Rnw:229-230
###################################################
(sum(DEU_PISA2012$cog$dat$SCIE$inc7))/(prod(dim(DEU_PISA2012$cog$dat$SCIE$inc7)))*100


