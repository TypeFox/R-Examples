### R code from vignette source 'kcirtTutorial.Snw'

###################################################
### code chunk number 1: kcirtTutorial.Snw:111-114
###################################################
options(width=49)
options(prompt=" ")
options(continue="   ")


###################################################
### code chunk number 2: kcirtTutorial.Snw:126-138 (eval = FALSE)
###################################################
## options(stringsAsFactors=FALSE, width=140)
## 
## library(kcirt)
## 
## constructMap.ls <- list(
## c(1,2,3),
## c(1,2,3),
## c(1,2,3),
## c(1,1,2,3),
## c(1,2,2,3),
## c(1,2,3,3)
## )


###################################################
### code chunk number 3: kcirtTutorial.Snw:146-156 (eval = FALSE)
###################################################
## set.seed(99999)
## 
## mxLambda <- c(
## c(1,1,-1),
## c(1,-1,1),
## c(-1,1,1),
## c(1,-1,1,-1),
## c(1,-1,1,-1),
## c(1,-1,1,-1)
## )


###################################################
### code chunk number 4: kcirtTutorial.Snw:165-166 (eval = FALSE)
###################################################
## covEta <- diag(1, 3)


###################################################
### code chunk number 5: kcirtTutorial.Snw:176-177 (eval = FALSE)
###################################################
## qTypes <- c("R", "R", "R", "M", "M", "M")


###################################################
### code chunk number 6: kcirtTutorial.Snw:185-186 (eval = FALSE)
###################################################
## mu <- rep(0, length(mxLambda))


###################################################
### code chunk number 7: kcirtTutorial.Snw:194-198 (eval = FALSE)
###################################################
## mod1 <- kcirt.model(constructMap.ls=constructMap.ls, qTypes=qTypes, mxLambda=mxLambda, 
## mu=mu, covEta=covEta)
## 
## mod1


###################################################
### code chunk number 8: kcirtTutorial.Snw:206-207 (eval = FALSE)
###################################################
## kcirt.ystarinfo(mod1)


###################################################
### code chunk number 9: kcirtTutorial.Snw:219-225 (eval = FALSE)
###################################################
## N <- 200
## 
## set.seed(99999)
## mod2 <- kcirt.sim(mod1, N=N)
## 
## mod2$mxData


###################################################
### code chunk number 10: kcirtTutorial.Snw:234-236 (eval = FALSE)
###################################################
## mxHatLambda <- matrix(0, nrow(mod2$mxLambda), ncol(mod2$mxLambda))
## diag(mxHatLambda) <- diag(mod2$mxLambda) + rnorm(nrow(mod2$mxLambda), 0, 0.2)


###################################################
### code chunk number 11: kcirtTutorial.Snw:245-247 (eval = FALSE)
###################################################
## hatMu <- rep(0, length(mxLambda))
## mxHatEta <- matrix(0, N, 3)


###################################################
### code chunk number 12: kcirtTutorial.Snw:254-257 (eval = FALSE)
###################################################
## mod2$mxHatLambda <- mxHatLambda
## mod2$mxHatEta <- mxHatEta
## mod2$hatMu <- hatMu


###################################################
### code chunk number 13: kcirtTutorial.Snw:263-264 (eval = FALSE)
###################################################
## mod3 <- mod2


###################################################
### code chunk number 14: kcirtTutorial.Snw:273-286 (eval = FALSE)
###################################################
## mod3 <- kcirt.fitMSS(model=mod3, lambdaConstraint="self", kcpus=2,
## penalty="logit", usetruesigma=TRUE, mss.sd=c(0.03, 0.03, 0.5))
## 
## mod3$performance
## mean(mod3$performance)
## 
## par(mfrow=c(1, 4))
## for(i in 1:ncol(mod3$mxHatEta)) {
## plot(mod3$mxEta[ ,i], mod3$mxHatEta[ ,i],
## main = round(1000*cor(mod3$mxHatEta[ ,i], mod3$mxEta[ ,i])^2)/1000 )
## }
## plot(mod3$mxLambda, mod3$mxHatLambda)
## 


###################################################
### code chunk number 15: kcirtTutorial.Snw:318-327 (eval = FALSE)
###################################################
## constructMap.ls <- list(
## c(1,1,2),
## c(2,2,3),
## c(3,3,1)
## )
## 
## qTypes <- c("R", "R", "R")
## 
## mod1 <- kcirt.model(constructMap.ls=constructMap.ls, qTypes=qTypes)


###################################################
### code chunk number 16: kcirtTutorial.Snw:336-337 (eval = FALSE)
###################################################
## mod1$mxLambdaCTinfo


###################################################
### code chunk number 17: kcirtTutorial.Snw:365-397 (eval = FALSE)
###################################################
## constructMap.ls <- list(
## c(1,2,3),
## c(1,2,3),
## c(1,2,3),
## c(1,1,2,3),
## c(1,2,2,3),
## c(1,2,3,3)
## )
## 
## 
## mxLambda <- c(
## c(1,1,-1),
## c(1,-1,1),
## c(-1,1,1),
## c(1,-1,1,-1),
## c(1,-1,1,-1),
## c(1,-1,1,-1)
## )
## 
## 
## covEta <- diag(1, 3)
## 
## 
## qTypes <- c("R", "R", "R", "M", "M", "M")
## 
## 
## mu <- rep(0, length(mxLambda))
## 
## 
## mod1 <- kcirt.model(constructMap.ls=constructMap.ls, qTypes=qTypes, mxLambda=mxLambda,
## mu=mu, covEta=covEta)
## 


###################################################
### code chunk number 18: kcirtTutorial.Snw:405-408 (eval = FALSE)
###################################################
## set.seed(99999)
## mod1$mxLambda[ mod1$mxLambdaCTinfo == "WF" ] <-
## rnorm( sum(mod1$mxLambdaCTinfo == "WF" ), 0, 0.1 )


###################################################
### code chunk number 19: kcirtTutorial.Snw:416-420 (eval = FALSE)
###################################################
## set.seed(99999)
## N <- 1000
## modX1 <- kcirt.sim(mod1, N=N)
## modX2 <- kcirt.sim(mod1, N=N)


###################################################
### code chunk number 20: kcirtTutorial.Snw:428-467 (eval = FALSE)
###################################################
## 
## mxHatLambda <- matrix(0, nrow(modX1$mxLambda), ncol(modX1$mxLambda))
## diag(mxHatLambda) <- diag(modX1$mxLambda) + rnorm(nrow(modX1$mxLambda), 0, 0.2)
## 
## mxHatEta <- matrix(0, N, 3)
## hatMu <- rep(0, length(mxLambda))
## 
## modX1$mxHatLambda <- mxHatLambda
## modX1$mxHatEta <- mxHatEta
## modX1$hatMu <- hatMu
## 
## 
## mod3 <- modX1
## 
## 
## mod3 <- kcirt.fitMSS(model=mod3, lambdaConstraint="self", kcpus=2,
## penalty="logit", usetruesigma=TRUE, mss.sd=c(0.03, 0.03, 0.5)) ### first MSS pass
## mean(mod3$performance)
## 
## mod3 <- kcirt.fitMSS(model=mod3, lambdaConstraint="self", kcpus=2,
## penalty="logit", usetruesigma=TRUE, mss.sd=c(0.03, 0.03, 0.5)) ### second MSS pass
## mean(mod3$performance)
## 
## mod3 <- kcirt.fitMSS(model=mod3, lambdaConstraint="self", kcpus=2,
## penalty="logit", usetruesigma=TRUE, mss.sd=c(0.03, 0.03, 0.5)) ### third MSS pass
## mean(mod3$performance)
## 
## 
## 
## mod3$performance
## 
## 
## par(mfrow=c(1, 4))
## for(i in 1:ncol(mod3$mxHatEta)) {
## plot(mod3$mxEta[ ,i], mod3$mxHatEta[ ,i],
## main = round(1000*cor(mod3$mxHatEta[ ,i], mod3$mxEta[ ,i])^2)/1000 )
## }
## plot(mod3$mxLambda, mod3$mxHatLambda)
## 


###################################################
### code chunk number 21: kcirtTutorial.Snw:480-513 (eval = FALSE)
###################################################
## 
## 
## modX2$hatMu <- mod3$hatMu
## modX2$mxHatLambda <- mod3$mxHatLambda
## modX2$mxHatEta <- matrix(0, N, 3)
## 
## mod4 <- modX2
## 
## 
## mod4 <- kcirt.fitMSS(model=mod4, lambdaConstraint="self", kcpus=2,
## penalty="logit", usetruesigma=TRUE, mss.sd=c(0, 0, 0.5)) ### first MSS pass
## mean(mod4$performance)
## 
## mod4 <- kcirt.fitMSS(model=mod4, lambdaConstraint="self", kcpus=2,
## penalty="logit", usetruesigma=TRUE, mss.sd=c(0, 0, 0.5)) ### second MSS pass
## mean(mod4$performance)
## 
## mod4 <- kcirt.fitMSS(model=mod4, lambdaConstraint="self", kcpus=2,
## penalty="logit", usetruesigma=TRUE, mss.sd=c(0, 0, 0.5)) ### third MSS pass
## mean(mod4$performance)
## 
## 
## 
## mod4$performance
## 
## 
## par(mfrow=c(1, 4))
## for(i in 1:ncol(mod4$mxHatEta)) {
## plot(mod4$mxEta[ ,i], mod4$mxHatEta[ ,i],
## main = round(1000*cor(mod4$mxHatEta[ ,i], mod4$mxEta[ ,i])^2)/1000 )
## }
## plot(mod4$mxLambda, mod4$mxHatLambda)
## 


###################################################
### code chunk number 22: kcirtTutorial.Snw:533-565 (eval = FALSE)
###################################################
## constructMap.ls <- list(
## c(1,2,3),
## c(1,2,3),
## c(1,2,3),
## c(1,1,2,3),
## c(1,2,2,3),
## c(1,2,3,3)
## )
## 
## 
## mxLambda <- c(
## c(1,1,-1),
## c(1,-1,1),
## c(-1,1,1),
## c(1,-1,1,-1),
## c(1,-1,1,-1),
## c(1,-1,1,-1)
## )
## 
## 
## covEta <- diag(1, 3)
## 
## 
## qTypes <- c("R", "R", "R", "M", "M", "M")
## 
## 
## mu <- rep(0, length(mxLambda))
## 
## 
## mod1 <- kcirt.model(constructMap.ls=constructMap.ls, qTypes=qTypes, mxLambda=mxLambda,
## mu=mu, covEta=covEta)
## 


###################################################
### code chunk number 23: kcirtTutorial.Snw:575-579 (eval = FALSE)
###################################################
## set.seed(99999)
## N <- 1000
## modX1 <- kcirt.sim(mod1, N=N)
## modX2 <- kcirt.sim(mod1, N=N)


###################################################
### code chunk number 24: kcirtTutorial.Snw:587-626 (eval = FALSE)
###################################################
## 
## mxHatLambda <- matrix(0, nrow(modX1$mxLambda), ncol(modX1$mxLambda))
## diag(mxHatLambda) <- diag(modX1$mxLambda) + rnorm(nrow(modX1$mxLambda), 0, 0.2)
## 
## mxHatEta <- matrix(0, N, 3)
## hatMu <- rep(0, length(mxLambda))
## 
## modX1$mxHatLambda <- mxHatLambda
## modX1$mxHatEta <- mxHatEta
## modX1$hatMu <- hatMu
## 
## 
## mod3 <- modX1
## 
## 
## mod3 <- kcirt.fitMSS(model=mod3, lambdaConstraint="withinx", kcpus=2,
## penalty="logit", usetruesigma=TRUE, mss.sd=c(0.03, 0.03, 0.5)) ### first MSS pass
## mean(mod3$performance)
## 
## mod3 <- kcirt.fitMSS(model=mod3, lambdaConstraint="withinx", kcpus=2,
## penalty="logit", usetruesigma=TRUE, mss.sd=c(0.03, 0.03, 0.5)) ### second MSS pass
## mean(mod3$performance)
## 
## mod3 <- kcirt.fitMSS(model=mod3, lambdaConstraint="withinx", kcpus=2,
## penalty="logit", usetruesigma=TRUE, mss.sd=c(0.03, 0.03, 0.5)) ### third MSS pass
## mean(mod3$performance)
## 
## 
## 
## mod3$performance
## 
## 
## par(mfrow=c(1, 4))
## for(i in 1:ncol(mod3$mxHatEta)) {
##     plot(mod3$mxEta[ ,i], mod3$mxHatEta[ ,i],
##     main = round(1000*cor(mod3$mxHatEta[ ,i], mod3$mxEta[ ,i])^2)/1000 )
## }
## plot(mod3$mxLambda, mod3$mxHatLambda)
## 


###################################################
### code chunk number 25: kcirtTutorial.Snw:640-673 (eval = FALSE)
###################################################
## 
## 
## modX2$hatMu <- mod3$hatMu
## modX2$mxHatLambda <- mod3$mxHatLambda
## modX2$mxHatEta <- matrix(0, N, 3)
## 
## mod4 <- modX2
## 
## 
## mod4 <- kcirt.fitMSS(model=mod4, lambdaConstraint="self", kcpus=2,
## penalty="logit", usetruesigma=TRUE, mss.sd=c(0, 0, 0.5)) ### first MSS pass
## mean(mod4$performance)
## 
## mod4 <- kcirt.fitMSS(model=mod4, lambdaConstraint="self", kcpus=2,
## penalty="logit", usetruesigma=TRUE, mss.sd=c(0, 0, 0.5)) ### second MSS pass
## mean(mod4$performance)
## 
## mod4 <- kcirt.fitMSS(model=mod4, lambdaConstraint="self", kcpus=2,
## penalty="logit", usetruesigma=TRUE, mss.sd=c(0, 0, 0.5)) ### third MSS pass
## mean(mod4$performance)
## 
## 
## 
## mod4$performance
## 
## 
## par(mfrow=c(1, 4))
## for(i in 1:ncol(mod4$mxHatEta)) {
##     plot(mod4$mxEta[ ,i], mod4$mxHatEta[ ,i],
##     main = round(1000*cor(mod4$mxHatEta[ ,i], mod4$mxEta[ ,i])^2)/1000 )
## }
## plot(mod4$mxLambda, mod4$mxHatLambda)
## 


###################################################
### code chunk number 26: kcirtTutorial.Snw:720-726 (eval = FALSE)
###################################################
## 
## 
## setwd("/Users/dzes/KF_Files/KF_Creations/KF_R/kcirt/inst/doc")
## Sweave("kcirt.Snw")
## 
## 


