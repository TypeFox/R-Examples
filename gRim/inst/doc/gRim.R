### R code from vignette source 'gRim.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: gRim.Rnw:72-76
###################################################
dir.create("figures")
oopt <- options()
options("digits"=4, "width"=80, "prompt"=" ", "continue"="  ")
options(useFancyQuotes="UTF-8")


###################################################
### code chunk number 2: gRim.Rnw:83-87
###################################################
options("width"=85)
library(gRim)
library(Rgraphviz)
ps.options(family="serif")


###################################################
### code chunk number 3: gRim.Rnw:119-122
###################################################
args(dmod)
args(cmod)
args(mmod)


###################################################
### code chunk number 4: gRim.Rnw:140-142
###################################################
data(reinis)
str(reinis)


###################################################
### code chunk number 5: gRim.Rnw:154-158
###################################################
data(reinis)
dm1<-dmod(list(c("smoke","systol"),c("smoke","mental","phys")), data=reinis)
dm1<-dmod(~smoke:systol + smoke:mental:phys, data=reinis)
dm1


###################################################
### code chunk number 6: gRim.Rnw:178-180
###################################################
formula(dm1)
str(terms(dm1))


###################################################
### code chunk number 7: gRim.Rnw:186-187
###################################################
summary(dm1)


###################################################
### code chunk number 8: gRim.Rnw:218-221
###################################################
dm2 <- dmod(~.^2, margin=c("smo","men","phy","sys"),
            data=reinis)
formula(dm2)


###################################################
### code chunk number 9: gRim.Rnw:225-228
###################################################
dm3 <- dmod(list(c("smoke", "systol"), c("smoke", "mental", "phys")),
            data=reinis, interactions=2)
formula(dm3)


###################################################
### code chunk number 10: gRim.Rnw:245-246
###################################################
iplot(dm1)


###################################################
### code chunk number 11: gRim.Rnw:260-264
###################################################
data(carcass)
cm1 <- cmod(~Fat11:Fat12:Fat13, data=carcass)
cm1 <- cmod(~Fat11:Fat12 + Fat12:Fat13 + Fat11:Fat13, data=carcass)
cm1


###################################################
### code chunk number 12: gRim.Rnw:270-271
###################################################
iplot(cm1)


###################################################
### code chunk number 13: gRim.Rnw:278-281
###################################################
data(milkcomp1)
mm1 <- mmod(~.^., data=milkcomp1)
mm1


###################################################
### code chunk number 14: gRim.Rnw:287-288
###################################################
iplot(mm1)


###################################################
### code chunk number 15: gRim.Rnw:305-307
###################################################
ms <- dmod(~.^., marginal=c("phys","mental","systol","family"), data=reinis)
formula(ms)


###################################################
### code chunk number 16: gRim.Rnw:313-315
###################################################
ms1 <- update(ms, list(dedge=~phys:mental))
formula(ms1)


###################################################
### code chunk number 17: gRim.Rnw:321-323
###################################################
ms2<- update(ms, list(dedge=~phys:mental+systol:family))
formula(ms2)


###################################################
### code chunk number 18: gRim.Rnw:329-331
###################################################
ms3 <- update(ms, list(dedge=~phys:mental:systol))
formula(ms3)


###################################################
### code chunk number 19: gRim.Rnw:337-339
###################################################
ms4 <- update(ms, list(dterm=~phys:mental:systol) )
formula(ms4)


###################################################
### code chunk number 20: gRim.Rnw:345-347
###################################################
ms5 <- update(ms, list(aterm=~phys:mental+phys:systol+mental:systol) )
formula(ms5)


###################################################
### code chunk number 21: gRim.Rnw:353-355
###################################################
ms6 <- update(ms, list(aedge=~phys:mental+systol:family))
formula(ms6)


###################################################
### code chunk number 22: gRim.Rnw:376-377
###################################################
cit <- ciTest(reinis, set=c("systol","smoke","family","phys"))


###################################################
### code chunk number 23: gRim.Rnw:410-411
###################################################
cit$slice


###################################################
### code chunk number 24: gRim.Rnw:443-444
###################################################
ciTest(reinis, set=c("systol","smoke","family","phys"), method='MC')


###################################################
### code chunk number 25: gRim.Rnw:458-460
###################################################
dm5 <- dmod(~ment:phys:systol+ment:systol:family+phys:systol:smoke,
            data=reinis)


###################################################
### code chunk number 26: fundamentalfig1
###################################################
iplot(dm5)


###################################################
### code chunk number 27: gRim.Rnw:485-487
###################################################
testdelete(dm5, ~smoke:systol)
testdelete(dm5, ~family:systol)


###################################################
### code chunk number 28: gRim.Rnw:506-507
###################################################
testadd(dm5, ~smoke:mental)


###################################################
### code chunk number 29: gRim.Rnw:533-534
###################################################
ed.in <- getInEdges(ugList(dm5$glist), type="decomposable")


###################################################
### code chunk number 30: gRim.Rnw:547-548
###################################################
ed.out <- getOutEdges(ugList(dm5$glist), type="decomposable")


###################################################
### code chunk number 31: gRim.Rnw:556-558
###################################################
args(testInEdges)
args(testOutEdges)


###################################################
### code chunk number 32: gRim.Rnw:572-574
###################################################
testInEdges(dm5, getInEdges(ugList(dm5$glist), type="decomposable"),
             k=log(sum(reinis)))


###################################################
### code chunk number 33: gRim.Rnw:596-599
###################################################
dm.sat <- dmod(~.^., data=reinis)
dm.back <- backward(dm.sat)
iplot(dm.back)


###################################################
### code chunk number 34: gRim.Rnw:614-617
###################################################
dm.i   <- dmod(~.^1, data=reinis)
dm.forw <- forward(dm.i)
iplot(dm.forw)


###################################################
### code chunk number 35: gRim.Rnw:681-685
###################################################
fix <- list(c("smoke","phys","systol"), c("systol","protein"))
fix <- do.call(rbind, unlist(lapply(fix, names2pairs),recursive=FALSE))
fix
dm.s3 <- backward(dm.sat, fixin=fix, details=1)


###################################################
### code chunk number 36: gRim.Rnw:696-697
###################################################
dm.i3 <- forward(dm.i, fixout=fix, details=1)


###################################################
### code chunk number 37: gRim.Rnw:903-904
###################################################
loglinGenDim(dm2$glist, reinis)


###################################################
### code chunk number 38: gRim.Rnw:958-961
###################################################
dm3 <- dmod(list(c("smoke", "systol"), c("smoke", "mental", "phys")),
            data=reinis)
names(dm3)


###################################################
### code chunk number 39: gRim.Rnw:969-970
###################################################
str(dm3$glist)


###################################################
### code chunk number 40: gRim.Rnw:974-975
###################################################
str(dm3$glistNUM)


###################################################
### code chunk number 41: gRim.Rnw:981-982
###################################################
dm3$varNames


###################################################
### code chunk number 42: gRim.Rnw:990-991
###################################################
str(dm3[c("varNames","conNames","conLevels")])


