### R code from vignette source 'kst.Rnw'

###################################################
### code chunk number 1: kst.Rnw:56-58
###################################################
options(width = 80)
library("kst")


###################################################
### code chunk number 2: kstructure
###################################################
# An endorelation representing a surmise relation
kst <- endorelation(graph=set(tuple(1,1), tuple(2,2), tuple(3,3),
  tuple(4,4), tuple(2,1), tuple(3,1), tuple(4,1), tuple(3,2), tuple(4,2)))
kstructure(kst)
# A set of sets representing knowledge states (e.g., clauses of a surmise system)
kst <- kstructure(set(set("a"), set("a","b"), set("a","c"), set("d","e"), 
   set("a","b","d","e"), set("a","c","d","e"), set("a","b","c","d","e")))
kst


###################################################
### code chunk number 3: set_options
###################################################
sets_options("quote",FALSE)
kst


###################################################
### code chunk number 4: kdomain
###################################################
kdomain(kst)


###################################################
### code chunk number 5: knotions
###################################################
knotions(kst)


###################################################
### code chunk number 6: katoms
###################################################
katoms(kst, items=set("a","b","c"))


###################################################
### code chunk number 7: ktrace
###################################################
ktrace(kst, items=set("c","d","e"))


###################################################
### code chunk number 8: kfringe
###################################################
kfringe(kst, operation="inner")
kfringe(kst, state=set("a", "b"), operation="inner")
kfringe(kst, operation="outer")
kfringe(kst, state=set("a", "b"), operation="outer")


###################################################
### code chunk number 9: kstructure_is_wellgraded
###################################################
kstructure_is_wellgraded(kst)


###################################################
### code chunk number 10: plot
###################################################
if(require("Rgraphviz")) {plot(kst)}


###################################################
### code chunk number 11: kst.Rnw:164-165
###################################################
if(require("Rgraphviz")) {plot(kst)}


###################################################
### code chunk number 12: as.relation
###################################################
as.relation(kst)


###################################################
### code chunk number 13: kassess
###################################################
rp <- data.frame(a=c(1,1,0,1,1,1,1,0,0,0),b=c(0,1,0,1,0,1,0,1,0,0),
   c=c(0,0,0,0,1,1,1,0,1,0),d=c(0,0,1,1,1,1,0,0,0,1), e=c(0,0,1,1,1,1,0,0,0,0))
kassess(kst, rpatterns=rp)


###################################################
### code chunk number 14: kvalidate
###################################################
kvalidate(kst, rpatterns=rp, method="gamma")
kvalidate(kst, rpatterns=rp, method="percent")
kvalidate(kst, rpatterns=rp, method="VC")
kvalidate(kst, rpatterns=rp, method="DA")


###################################################
### code chunk number 15: closure
###################################################
closure(kst, operation="union")


###################################################
### code chunk number 16: reduction
###################################################
reduction(kst, operation="discrimination")


###################################################
### code chunk number 17: kspace
###################################################
ksp <- kspace(kst)
ksp


###################################################
### code chunk number 18: kstructure_is_space
###################################################
kstructure_is_kspace(ksp)


###################################################
### code chunk number 19: kbase
###################################################
kbase(ksp)


###################################################
### code chunk number 20: lpath
###################################################
lp <- lpath(ksp)
lp


###################################################
### code chunk number 21: lpath_is_gradation
###################################################
lpath_is_gradation(lp)


