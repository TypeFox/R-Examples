### R code from vignette source 'BinomialInference.Rnw'

###################################################
### code chunk number 1: BinomialInference.Rnw:17-20
###################################################
library(LearnBayes)
beta.par <- beta.select(list(p=0.5, x=0.2), list(p=0.75, x=.28))
beta.par


###################################################
### code chunk number 2: BinomialInference.Rnw:32-33
###################################################
triplot(beta.par, c(6, 4))


###################################################
### code chunk number 3: BinomialInference.Rnw:40-43
###################################################
beta.post.par <- beta.par + c(6, 4)
post.sample <- rbeta(1000, beta.post.par[1], beta.post.par[2])
quantile(post.sample, c(0.05, 0.95))


###################################################
### code chunk number 4: BinomialInference.Rnw:50-51
###################################################
predplot(beta.par, 10, 6)


###################################################
### code chunk number 5: BinomialInference.Rnw:60-65
###################################################
n <- 20
s <- 0:n
pred.probs <- pbetap(beta.par, n, s)
plot(s, pred.probs, type="h")
discint(cbind(s, pred.probs), 0.90)


