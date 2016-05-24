### R code from vignette source 'PigeonExample.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: PigeonExample.Rnw:21-22
###################################################
options(width=67)


###################################################
### code chunk number 2: PigeonExample.Rnw:26-28
###################################################
library(indicspecies)
data(pigeons)


###################################################
### code chunk number 3: PigeonExample.Rnw:31-32
###################################################
ls()


###################################################
### code chunk number 4: PigeonExample.Rnw:36-37
###################################################
dfood


###################################################
### code chunk number 5: PigeonExample.Rnw:40-42
###################################################
plot(hclust(dfood, method="average"), h=-1, xlab="", 
     ylab="Distance", main="", sub="", ylim=c(0.2,1))


###################################################
### code chunk number 6: PigeonExample.Rnw:49-53
###################################################
diet.pop.barcelona = colSums(diet.barcelona)
round(diet.pop.barcelona/sum(diet.pop.barcelona), dig=3)
diet.pop.moia = colSums(diet.moia)
round(diet.pop.moia/sum(diet.pop.moia), dig=3)


###################################################
### code chunk number 7: PigeonExample.Rnw:64-66
###################################################
nichevar(P=diet.barcelona, mode="single")
nichevar(P=diet.moia, mode="single")


###################################################
### code chunk number 8: PigeonExample.Rnw:69-74
###################################################
popvar.barcelona = nichevar(P=diet.barcelona, D=dfood, 
                            mode="single")
popvar.barcelona
popvar.moia= nichevar(P=diet.moia, D=dfood, mode="single")
popvar.moia


###################################################
### code chunk number 9: PigeonExample.Rnw:84-86
###################################################
nicheoverlap(P1=diet.barcelona, P2=diet.moia, mode="single")
nicheoverlap(P1=diet.barcelona, P2=diet.moia, mode="single", D = dfood)


###################################################
### code chunk number 10: PigeonExample.Rnw:94-96
###################################################
round(sweep(diet.barcelona, 1, FUN="/", 
            rowSums(diet.barcelona)), dig=3)


###################################################
### code chunk number 11: PigeonExample.Rnw:99-101
###################################################
round(sweep(diet.moia, 1, FUN="/", 
            rowSums(diet.moia)), dig=3)


###################################################
### code chunk number 12: PigeonExample.Rnw:111-115
###################################################
indvar.barcelona<-nichevar(P=diet.barcelona, D=dfood)
summary(indvar.barcelona)
indvar.moia<-nichevar(P=diet.moia, D=dfood)
summary(indvar.moia)


###################################################
### code chunk number 13: PigeonExample.Rnw:120-121
###################################################
wilcox.test(indvar.barcelona$B, indvar.moia$B)


###################################################
### code chunk number 14: PigeonExample.Rnw:130-134
###################################################
Spec.barcelona = mean(indvar.barcelona$B)/popvar.barcelona$B
Spec.barcelona
Spec.moia = mean(indvar.moia$B)/popvar.moia$B
Spec.moia


###################################################
### code chunk number 15: PigeonExample.Rnw:141-143
###################################################
Spec.ind.barcelona = indvar.barcelona$B/popvar.barcelona$B
Spec.ind.moia = indvar.moia$B/popvar.moia$B


###################################################
### code chunk number 16: PigeonExample.Rnw:146-147
###################################################
wilcox.test(Spec.ind.barcelona, Spec.ind.moia)


###################################################
### code chunk number 17: PigeonExample.Rnw:153-155
###################################################
O.barcelona = nicheoverlap(diet.barcelona,D=dfood, mode="pairwise")
O.moia = nicheoverlap(diet.moia,D=dfood, mode="pairwise")


###################################################
### code chunk number 18: PigeonExample.Rnw:158-160
###################################################
mean(O.barcelona[lower.tri(O.barcelona)])
mean(O.moia[lower.tri(O.moia)])


###################################################
### code chunk number 19: PigeonExample.Rnw:163-167
###################################################
O.barcelona.ind = (rowSums(O.barcelona)-1)/(nrow(O.barcelona)-1)
summary(O.barcelona.ind)
O.moia.ind = (rowSums(O.moia)-1)/(nrow(O.moia)-1)
summary(O.moia.ind)


###################################################
### code chunk number 20: PigeonExample.Rnw:170-171
###################################################
wilcox.test(O.barcelona.ind, O.moia.ind)


