### R code from vignette source 'multivariate-knee.Rnw'

###################################################
### code chunk number 1: multivariate-knee.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: multivariate-knee.Rnw:19-22 (eval = FALSE)
###################################################
## library(catdata)
## data(knee)
## attach(knee)


###################################################
### code chunk number 3: multivariate-knee.Rnw:29-36 (eval = FALSE)
###################################################
## R2D <- rep(0, length(R2))
## R3D <- rep(0, length(R3))
## R4D <- rep(0, length(R3))
## 
## R2D[R2>2] <- 1
## R3D[R3>2] <- 1
## R4D[R4>2] <- 1


###################################################
### code chunk number 4: multivariate-knee.Rnw:42-46 (eval = FALSE)
###################################################
## N <- rep(knee$N, each=3)
## Th <- rep(knee$Th, each=3)
## Age <- rep(knee$Age, each=3)
## Sex <- rep(knee$Sex, each=3)


###################################################
### code chunk number 5: multivariate-knee.Rnw:51-53 (eval = FALSE)
###################################################
## Response <- c(rbind(R2D,R3D,R4D))
## Age2 <- Age^2


###################################################
### code chunk number 6: multivariate-knee.Rnw:58-60 (eval = FALSE)
###################################################
## Th <- as.factor(Th)
## Sex <- as.factor(Sex)


###################################################
### code chunk number 7: multivariate-knee.Rnw:65-66 (eval = FALSE)
###################################################
## library(gee)


###################################################
### code chunk number 8: multivariate-knee.Rnw:71-73 (eval = FALSE)
###################################################
## gee1a <- gee(Response ~ Th + Sex + Age + Age2, id=N,
## family=binomial(link=logit))


###################################################
### code chunk number 9: multivariate-knee.Rnw:76-77 (eval = FALSE)
###################################################
## summary(gee1a)


###################################################
### code chunk number 10: multivariate-knee.Rnw:82-84 (eval = FALSE)
###################################################
## gee2a <- gee(Response ~ Th + Sex + Age + Age2, id=N,
## family=binomial(link=logit), corstr="exchangeable")


###################################################
### code chunk number 11: multivariate-knee.Rnw:87-88 (eval = FALSE)
###################################################
## summary(gee2a)


###################################################
### code chunk number 12: multivariate-knee.Rnw:93-95 (eval = FALSE)
###################################################
## gee3a <- gee(Response ~ Th + Sex + Age + Age2, id=N,
## family=binomial(link=logit), corstr="AR-M", Mv=1)


###################################################
### code chunk number 13: multivariate-knee.Rnw:98-99 (eval = FALSE)
###################################################
## summary(gee3a)


###################################################
### code chunk number 14: multivariate-knee.Rnw:105-106 (eval = FALSE)
###################################################
## library(geepack)


###################################################
### code chunk number 15: multivariate-knee.Rnw:111-113 (eval = FALSE)
###################################################
## gee1b <- geeglm(Response ~ Th + Sex + Age + Age2, id=N,
## family=binomial(link=logit))


###################################################
### code chunk number 16: multivariate-knee.Rnw:116-117 (eval = FALSE)
###################################################
## summary(gee1b)


###################################################
### code chunk number 17: multivariate-knee.Rnw:122-124 (eval = FALSE)
###################################################
## gee2b <- geeglm(Response ~ Th + Sex + Age + Age2, id=N,
## family=binomial(link=logit), corstr="exchangeable")


###################################################
### code chunk number 18: multivariate-knee.Rnw:127-128 (eval = FALSE)
###################################################
## summary(gee2b)


###################################################
### code chunk number 19: multivariate-knee.Rnw:133-135 (eval = FALSE)
###################################################
## gee3b <- geeglm(Response ~ Th + Sex + Age + Age2, id=N,
## family=binomial(link=logit), corstr="ar1")


###################################################
### code chunk number 20: multivariate-knee.Rnw:138-139 (eval = FALSE)
###################################################
## summary(gee3b)


###################################################
### code chunk number 21: multivariate-knee.Rnw:146-149 (eval = FALSE)
###################################################
## glm1 <- glm(Response ~ Th + Sex + Age + Age2,
## family=binomial(link=logit))
## summary(glm1)


###################################################
### code chunk number 22: multivariate-knee.Rnw:155-157 (eval = FALSE)
###################################################
## Age <- Age-30
## Age2 <- Age^2


###################################################
### code chunk number 23: multivariate-knee.Rnw:166-168 (eval = FALSE)
###################################################
## gee1c <- gee(Response ~ Th + Sex + Age + Age2, id=N,
## family=binomial(link=logit))


###################################################
### code chunk number 24: multivariate-knee.Rnw:171-172 (eval = FALSE)
###################################################
## summary(gee1c)


###################################################
### code chunk number 25: multivariate-knee.Rnw:177-179 (eval = FALSE)
###################################################
## gee2c <- gee(Response ~ Th + Sex + Age + Age2, id=N,
## family=binomial(link=logit), corstr="exchangeable")


###################################################
### code chunk number 26: multivariate-knee.Rnw:182-183 (eval = FALSE)
###################################################
## summary(gee2c)


###################################################
### code chunk number 27: multivariate-knee.Rnw:188-190 (eval = FALSE)
###################################################
## gee3c <- gee(Response ~ Th + Sex + Age + Age2, id=N,
## family=binomial(link=logit), corstr="AR-M", Mv=1)


###################################################
### code chunk number 28: multivariate-knee.Rnw:193-194 (eval = FALSE)
###################################################
## summary(gee3c)


