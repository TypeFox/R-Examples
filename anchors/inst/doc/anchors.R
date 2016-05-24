### R code from vignette source 'anchors.Rnw'

###################################################
### code chunk number 1: R options
###################################################
options(width = 60)
options(SweaveHooks = list(fig = function() par(mar=c(3,3,1,0.5),mgp = c(2,1,0))))


###################################################
### code chunk number 2: anchors.Rnw:135-136
###################################################
library(anchors)


###################################################
### code chunk number 3: anchors.Rnw:140-141 (eval = FALSE)
###################################################
## help(package="anchors")


###################################################
### code chunk number 4: anchors.Rnw:144-145 (eval = FALSE)
###################################################
##   data(package="anchors")


###################################################
### code chunk number 5: anchors.Rnw:149-150 (eval = FALSE)
###################################################
## demo(package="anchors")


###################################################
### code chunk number 6: anchors.Rnw:189-190 (eval = FALSE)
###################################################
## data(freedom)


###################################################
### code chunk number 7: anchors.Rnw:214-215 (eval = FALSE)
###################################################
## demo(anchors.freedom)


###################################################
### code chunk number 8: anchors.Rnw:257-258
###################################################
library(anchors)


###################################################
### code chunk number 9: anchors.Rnw:263-264 (eval = FALSE)
###################################################
## help(package="anchors")


###################################################
### code chunk number 10: anchors.Rnw:268-269 (eval = FALSE)
###################################################
## help(anchors)


###################################################
### code chunk number 11: anchors.Rnw:426-430
###################################################
library(anchors)
data(freedom)
a1 <- anchors(self ~ vign2+vign3+vign4+vign5+vign6, freedom, method="C")
summary(a1)


###################################################
### code chunk number 12: vignette-order
###################################################
getOption("SweaveHooks")[["fig"]]()
vo1<-anchors.order(~vign2+vign3+vign4+vign5+vign6, freedom)
summary(vo1,top=10,digits=3)
barplot(vo1)


###################################################
### code chunk number 13: anchors.Rnw:465-467
###################################################
a2 <- anchors(self ~ vign2+vign3+vign5+vign4+vign6, freedom, method="C")
summary(a2)


###################################################
### code chunk number 14: anchors.Rnw:489-494
###################################################
data(freedom)
fo <- list(self = self ~ 1,
           vign = cbind(vign1,vign3,vign6) ~ 1, 
           cpolr= ~ as.factor(country) + sex + age + educ)
ent <- anchors(fo, data = freedom, method="C", combn=TRUE)


###################################################
### code chunk number 15: anchors.Rnw:496-497
###################################################
summary(ent,digits=3)


###################################################
### code chunk number 16: entropy1
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(ent)


###################################################
### code chunk number 17: anchors.Rnw:721-724
###################################################
fo <- list(self =  self ~ sex + age + educ + factor(country) ,
             vign = cbind(vign1,vign2,vign3,vign4,vign5,vign6) ~ 1 ,
             tau  =       ~ sex + age + educ + factor(country) )


###################################################
### code chunk number 18: anchors.Rnw:744-745
###################################################
out  <- chopit( fo, data=freedom)


###################################################
### code chunk number 19: anchors.Rnw:748-749
###################################################
summary(out)


