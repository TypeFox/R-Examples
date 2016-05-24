### R code from vignette source 'modbin-unemployment1.Rnw'

###################################################
### code chunk number 1: modbin-unemployment1.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=60)


###################################################
### code chunk number 2: modbin-unemployment1.Rnw:18-21 (eval = FALSE)
###################################################
## library(catdata)
## data(unemployment)
## attach(unemployment)


###################################################
### code chunk number 3: modbin-unemployment1.Rnw:26-34 (eval = FALSE)
###################################################
## durbin <- as.factor(durbin)
## table.durbin <- ftable(subset(unemployment, select=c("age", "durbin")), 
## col.vars="durbin")
## rels<-table.durbin[,1]/rowSums(table.durbin)
## age.new <- min(age):max(age)
## 
## model1 <- glm(table.durbin ~ age.new, family=binomial)
## summary(model1)


###################################################
### code chunk number 4: modbin-unemployment1.Rnw:39-42 (eval = FALSE)
###################################################
## plot(age.new, model1$fitted.values, xlab="Age", ylab="Observed/Fitted values", 
## type="l", ylim=c(0,1))
## points(age.new,table.durbin[,1]/rowSums(table.durbin))


###################################################
### code chunk number 5: modbin-unemployment1.Rnw:47-49 (eval = FALSE)
###################################################
## plot(model1$fitted.values,sqrt(rowSums(table.durbin))*rstandard(model1), 
## xlab="Predicted values", ylab="Residuals")


###################################################
### code chunk number 6: modbin-unemployment1.Rnw:51-55 (eval = FALSE)
###################################################
## qqnorm(sqrt(rowSums(table.durbin))*rstandard(model1), main="", 
##        ylab="Standardized deviance residuals")
## qqline(sqrt(rowSums(table.durbin))*rstandard(model1), lwd=2.5, 
##        lty="dashed", col="red")


