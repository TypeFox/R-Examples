###################################################
### chunk number 1: 
###################################################
library(tgp)
library(maptree)
#options(width=65)
seed <- 0; set.seed(seed)


###################################################
### chunk number 2: 
###################################################
exp2d.data <- exp2d.rand(lh=0, dopt=10)
X <- exp2d.data$X
Z <- exp2d.data$Z
Xcand <- lhs(1000, rbind(c(-2,6),c(-2,6)))


###################################################
### chunk number 3: 
###################################################
exp1 <- btgpllm(X=X, Z=Z, pred.n=FALSE, corr="exp", verb=0)


###################################################
### chunk number 4: mapt
###################################################
tgp.trees(exp1)


###################################################
### chunk number 5: 
###################################################
rl <- readline("press RETURN to continue: ")
graphics.off()


###################################################
### chunk number 6: 
###################################################
XX <- tgp.design(200, Xcand, exp1)
XX <- rbind(XX, c(-sqrt(1/2),0))


###################################################
### chunk number 7: cands
###################################################
plot(exp1$X, pch=19, cex=0.5)
points(XX)
mapT(exp1, add=TRUE)


###################################################
### chunk number 8: 
###################################################
rl <- readline("press RETURN to continue: ")
graphics.off()


###################################################
### chunk number 9: 
###################################################
exp.as <- btgpllm(X=X, Z=Z, XX=XX, corr="exp", improv=TRUE, 
                        Ds2x=TRUE, verb=0)


###################################################
### chunk number 10: expas
###################################################
par(mfrow=c(1,3), bty="n")
plot(exp.as, main="tgpllm,", layout="as", as="alm")
plot(exp.as, main="tgpllm,", layout='as', as='alc')
plot(exp.as, main="tgpllm,", layout='as', as='improv')


###################################################
### chunk number 11: 
###################################################
rl <- readline("press RETURN to continue: ")
graphics.off()


