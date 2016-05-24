### R code from vignette source 'random-knee2.Rnw'

###################################################
### code chunk number 1: random-knee2.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: random-knee2.Rnw:21-24 (eval = FALSE)
###################################################
## library(catdata)
## data(kneesequential)
## data(kneecumulative)


###################################################
### code chunk number 3: random-knee2.Rnw:29-34 (eval = FALSE)
###################################################
## kneesequential$Age <- kneesequential$Age - 30
## kneesequential$Age2 <- kneesequential$Age^2
## 
## kneecumulative$Age <- kneecumulative$Age - 30
## kneecumulative$Age2<-kneecumulative$Age^2


###################################################
### code chunk number 4: random-knee2.Rnw:39-40 (eval = FALSE)
###################################################
## library(lme4)


###################################################
### code chunk number 5: random-knee2.Rnw:45-48 (eval = FALSE)
###################################################
## seqGH<-glmer(y~-1+Icept1+Icept2+Icept3+Icept4+Th+Age+Age2+(1|Person), 
##              family=binomial(link=logit),data=kneesequential, nAGQ = 25)
## summary(seqGH)


###################################################
### code chunk number 6: random-knee2.Rnw:53-54 (eval = FALSE)
###################################################
## library(MASS)


###################################################
### code chunk number 7: random-knee2.Rnw:59-62 (eval = FALSE)
###################################################
## seqPQL<-glmmPQL(y ~-1+Icept1+Icept2+Icept3+Icept4+Th+Age+Age2, 
## random=list(Person=~1), family=binomial(link=logit), data=kneesequential, niter=30)
## summary(seqPQL)


###################################################
### code chunk number 8: random-knee2.Rnw:67-68 (eval = FALSE)
###################################################
## library(ordinal)


###################################################
### code chunk number 9: random-knee2.Rnw:74-77 (eval = FALSE)
###################################################
## cumGH<-clmm2(as.factor(y)~1+Th+Age+Age2, random = as.factor(Person), data = 
## kneecumulative, link = "logistic",nAGQ=25,start=c(-5,-3,3,5,rep(0.001,4)),Hess=TRUE)
## summary(cumGH)


###################################################
### code chunk number 10: random-knee2.Rnw:82-85 (eval = FALSE)
###################################################
## cumLP<-clmm2(as.factor(y)~1+Th+Age+Age2, random = as.factor(Person), data = 
## kneecumulative, link = "logistic",start=c(-5,-3,3,5,rep(0.001,4)), Hess = TRUE)
## summary(cumLP)


###################################################
### code chunk number 11: random-knee2.Rnw:88-89 (eval = FALSE)
###################################################
## detach(package:ordinal)


###################################################
### code chunk number 12: random-knee2.Rnw:92-93 (eval = FALSE)
###################################################
## detach(package:lme4)


