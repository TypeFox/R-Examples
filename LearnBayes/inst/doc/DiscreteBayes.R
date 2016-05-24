### R code from vignette source 'DiscreteBayes.Rnw'

###################################################
### code chunk number 1: DiscreteBayes.Rnw:21-23
###################################################
p <- seq(0, 1, by = 0.01)
prior <- 1 / 101 + 0 * p


###################################################
### code chunk number 2: DiscreteBayes.Rnw:25-28
###################################################
plot(p, prior, 
     type="h",
     main="Prior Distribution")


###################################################
### code chunk number 3: DiscreteBayes.Rnw:35-37
###################################################
library(LearnBayes)
post <- pdisc(p, prior, c(20, 12))


###################################################
### code chunk number 4: DiscreteBayes.Rnw:39-42
###################################################
plot(p, post, 
     type="h",
     main="Posterior Distribution")


###################################################
### code chunk number 5: DiscreteBayes.Rnw:47-48
###################################################
discint(cbind(p, post), 0.90)


###################################################
### code chunk number 6: DiscreteBayes.Rnw:57-60
###################################################
n <- 20
s <- 0:20
pred.probs <- pdiscp(p, post, n, s)


###################################################
### code chunk number 7: DiscreteBayes.Rnw:63-66
###################################################
plot(s, pred.probs, 
     type="h",
     main="Predictive Distribution")


###################################################
### code chunk number 8: DiscreteBayes.Rnw:72-74
###################################################
prior <- rep(1/11, 11)
names(prior) <- 20:30


###################################################
### code chunk number 9: DiscreteBayes.Rnw:78-79
###################################################
y <- c(24, 25, 31, 31, 22, 21, 26, 20, 16, 22)


###################################################
### code chunk number 10: DiscreteBayes.Rnw:83-84
###################################################
post <- discrete.bayes(dpois, prior, y)


###################################################
### code chunk number 11: DiscreteBayes.Rnw:89-90
###################################################
print(post)


###################################################
### code chunk number 12: DiscreteBayes.Rnw:93-94
###################################################
plot(post)


###################################################
### code chunk number 13: DiscreteBayes.Rnw:97-98
###################################################
summary(post)


