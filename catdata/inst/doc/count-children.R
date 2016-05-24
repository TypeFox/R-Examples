### R code from vignette source 'count-children.Rnw'

###################################################
### code chunk number 1: count-children.Rnw:11-12 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: count-children.Rnw:18-21 (eval = FALSE)
###################################################
## library(catdata)
## data(children)
## attach(children)


###################################################
### code chunk number 3: count-children.Rnw:25-28 (eval = FALSE)
###################################################
## pois <- glm(child ~ age+I(age^2)+I(age^3)+I(age^4)+dur+I(dur^2)+nation+god+univ, 
##             data = children,  family = poisson(link=log))
## summary(pois)


###################################################
### code chunk number 4: count-children.Rnw:32-38 (eval = FALSE)
###################################################
## x <- min(age):max(age)
## y <- exp(pois$coef[1]+pois$coef["age"]*x+pois$coef["I(age^2)"]*x^2+
##   pois$coef["I(age^3)"]*x^3+pois$coef["I(age^4)"]*x^4+pois$coef["dur"]*10+
##   pois$coef["I(dur^2)"]*100)
## par(cex=1.4)
## plot(x, y, ylab="Number of Children", xlab="Age")


###################################################
### code chunk number 5: count-children.Rnw:41-46 (eval = FALSE)
###################################################
## y <- (pois$coef[1]+pois$coef["age"]*x+pois$coef["I(age^2)"]*x^2+
##   pois$coef["I(age^3)"]*x^3+pois$coef["I(age^4)"]*x^4+pois$coef["dur"]*10+
##   pois$coef["I(dur^2)"]*100)
## par(cex=1.4)
## plot(x, y, ylab="Linear Predictor", xlab="Age")


###################################################
### code chunk number 6: count-children.Rnw:49-55 (eval = FALSE)
###################################################
## x <- min(dur):max(dur)
## y <- exp(pois$coef[1]+pois$coef["age"]*40+pois$coef["I(age^2)"]*40^2+
##   pois$coef["I(age^3)"]*40^3+pois$coef["I(age^4)"]*40^4+pois$coef["dur"]*x+
##   pois$coef["I(dur^2)"]*x^2)
## par(cex=1.4)
## plot(x, y, ylab="Number of Children", xlab="Duration of School Education")


###################################################
### code chunk number 7: count-children.Rnw:58-63 (eval = FALSE)
###################################################
## y <- (pois$coef[1]+pois$coef["age"]*40+pois$coef["I(age^2)"]*40^2+
##   pois$coef["I(age^3)"]*40^3+pois$coef["I(age^4)"]*40^4+pois$coef["dur"]*x+
##   pois$coef["I(dur^2)"]*x^2)
## par(cex=1.4)
## plot(x, y, ylab="Linear Predictor", xlab="Duration of School Education")


###################################################
### code chunk number 8: count-children.Rnw:69-70 (eval = FALSE)
###################################################
## anova(pois)


