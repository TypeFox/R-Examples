### R code from vignette source 'Trees.Rnw'

###################################################
### code chunk number 1: Trees.Rnw:49-52
###################################################
options(width=70)
foo <- packageDescription("phangorn")
options("show.signif.stars" = FALSE)


###################################################
### code chunk number 2: Trees.Rnw:68-71
###################################################
library(phangorn)
primates <- read.phyDat("primates.dna", format="phylip", 
    type="DNA")


###################################################
### code chunk number 3: Trees.Rnw:77-80
###################################################
dm  <- dist.ml(primates)
treeUPGMA  <- upgma(dm)
treeNJ  <- NJ(dm)


###################################################
### code chunk number 4: plotNJ
###################################################
layout(matrix(c(1,2), 2, 1), height=c(1,2))
par(mar = c(0,0,2,0)+ 0.1)
plot(treeUPGMA, main="UPGMA")
plot(treeNJ, "unrooted", main="NJ")


###################################################
### code chunk number 5: figNJ
###################################################
getOption("SweaveHooks")[["fig"]]()
layout(matrix(c(1,2), 2, 1), height=c(1,2))
par(mar = c(0,0,2,0)+ 0.1)
plot(treeUPGMA, main="UPGMA")
plot(treeNJ, "unrooted", main="NJ")


###################################################
### code chunk number 6: Trees.Rnw:102-104
###################################################
parsimony(treeUPGMA, primates)
parsimony(treeNJ, primates)


###################################################
### code chunk number 7: Trees.Rnw:107-110
###################################################
treePars  <- optim.parsimony(treeUPGMA, primates)
treeRatchet  <- pratchet(primates, trace = 0)
parsimony(c(treePars, treeRatchet), primates)


###################################################
### code chunk number 8: Trees.Rnw:113-114 (eval = FALSE)
###################################################
## (trees <- bab(subset(primates,1:10)))


###################################################
### code chunk number 9: Trees.Rnw:120-122
###################################################
fit = pml(treeNJ, data=primates)
fit


###################################################
### code chunk number 10: Trees.Rnw:125-126
###################################################
methods(class="pml")


###################################################
### code chunk number 11: Trees.Rnw:129-131
###################################################
fitJC  <- optim.pml(fit, TRUE)
logLik(fitJC)


###################################################
### code chunk number 12: Trees.Rnw:134-138
###################################################
fitGTR <- update(fit, k=4, inv=0.2) 
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, 
    rearrangement = "NNI", control = pml.control(trace = 0))
fitGTR 


###################################################
### code chunk number 13: Trees.Rnw:143-146
###################################################
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, 
    rearrangement = "stochastic", control = pml.control(trace = 0))
fitGTR 


###################################################
### code chunk number 14: Trees.Rnw:150-151
###################################################
anova(fitJC, fitGTR) 


###################################################
### code chunk number 15: Trees.Rnw:154-155
###################################################
SH.test(fitGTR, fitJC) 


###################################################
### code chunk number 16: Trees.Rnw:158-162
###################################################
AIC(fitJC)
AIC(fitGTR) 
AICc(fitGTR) 
BIC(fitGTR) 


###################################################
### code chunk number 17: Trees.Rnw:165-166
###################################################
load("Trees.RData")


###################################################
### code chunk number 18: Trees.Rnw:168-169 (eval = FALSE)
###################################################
## mt = modelTest(primates)


###################################################
### code chunk number 19: Trees.Rnw:173-175
###################################################
library(xtable)
print(xtable(mt, caption="Summary table of modelTest", label="tab:modelTest"), include.rownames=FALSE)


###################################################
### code chunk number 20: Trees.Rnw:179-182
###################################################
env <- attr(mt, "env")
ls(envir=env)
(fit <- eval(get("HKY+G+I", env), env))


###################################################
### code chunk number 21: Trees.Rnw:185-187 (eval = FALSE)
###################################################
## bs = bootstrap.pml(fitJC, bs=100, optNni=TRUE, 
##     control = pml.control(trace = 0))


###################################################
### code chunk number 22: plotBS
###################################################
par(mfrow=c(2,1))
par(mar=c(1,1,3,1))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
title("a)")
cnet <- consensusNet(bs, p=0.2)
plot(cnet, "2D", show.edge.label=TRUE)
title("b)")


###################################################
### code chunk number 23: figBS
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(2,1))
par(mar=c(1,1,3,1))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
title("a)")
cnet <- consensusNet(bs, p=0.2)
plot(cnet, "2D", show.edge.label=TRUE)
title("b)")


###################################################
### code chunk number 24: Trees.Rnw:218-220
###################################################
options(prompt=" ")
options(continue="  ")


###################################################
### code chunk number 25: Trees.Rnw:222-251 (eval = FALSE)
###################################################
## library(phangorn)
## file="myfile"
## dat = read.phyDat(file)
## dm = dist.ml(dat, "F81")
## tree = NJ(dm)
## # as alternative for a starting tree:
## tree <- pratchet(dat)          # parsimony tree 
## tree <- nnls.phylo(tree, dm)   # need edge weights
## 
## 
## # 1. alternative: quick and dirty: GTR + G 
## fitStart = pml(tree, dat, k=4)
## fit = optim.pml(fitStart, model="GTR", optGamma=TRUE, rearrangement="stochastic") 
##  
## # 2. alternative: preper with modelTest  
## mt <- modelTest(dat, tree=tree, multicore=TRUE)
## mt[order(mt$AICc),]
## # choose best model from the table according to AICc
## bestmodel <- mt$Model[which.min(mt$AICc)]
## 
## env = attr(mt, "env")
## fitStart = eval(get("GTR+G+I", env), env) 
## 
## # or let R search the table
## fitStart = eval(get(bestmodel, env), env) 
## # equivalent to:   fitStart = eval(get("GTR+G+I", env), env) 
## fit = optim.pml(fitStart, rearrangement = "stochastic", 
##     optGamma=TRUE, optInv=TRUE, model="GTR")
## bs = bootstrap.pml(fit, bs=100, optNni=TRUE, multicore=TRUE)


###################################################
### code chunk number 26: Trees.Rnw:257-278 (eval = FALSE)
###################################################
## library(phangorn)
## file="myfile"
## dat = read.phyDat(file, type = "AA")
## dm = dist.ml(dat, model="JTT")
## tree = NJ(dm)
## 
## # parallel will only work safely from command line 
## # and not at all windows
## (mt <- modelTest(dat, model=c("JTT", "LG", "WAG"), 
##     multicore=TRUE)) 
## # run all available amino acid models
## (mt <- modelTest(dat, model="all", multicore=TRUE))
## 
## fitStart = eval(get(mt$Model[which.min(mt$BIC)], env), env) 
## 
## fitNJ = pml(tree, dat, model="JTT", k=4, inv=.2)
## fit = optim.pml(fitNJ, rearrangement = "stochastic", 
##     optInv=TRUE, optGamma=TRUE)
## fit
## 
## bs = bootstrap.pml(fit, bs=100, optNni=TRUE, multicore=TRUE)


###################################################
### code chunk number 27: Trees.Rnw:290-291
###################################################
toLatex(sessionInfo())


