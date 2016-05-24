### R code from vignette source 'multinomial-addiction2.Rnw'

###################################################
### code chunk number 1: multinomial-addiction2.Rnw:13-15 (eval = FALSE)
###################################################
## rm(list=ls(all=TRUE))
## options(width=60)


###################################################
### code chunk number 2: multinomial-addiction2.Rnw:20-23 (eval = FALSE)
###################################################
## library(catdata)
## data(addiction)
## attach(addiction)


###################################################
### code chunk number 3: multinomial-addiction2.Rnw:29-34 (eval = FALSE)
###################################################
## ill01 <- ill
## ill01[ill==0] <- 1
## ill01[ill==2] <- 0
## 
## age2 <- age^2


###################################################
### code chunk number 4: multinomial-addiction2.Rnw:39-42 (eval = FALSE)
###################################################
## m01vs2 <- glm(ill01 ~ as.factor(gender) + as.factor(university) + age + age2, 
## family=binomial())
## summary(m01vs2)


###################################################
### code chunk number 5: multinomial-addiction2.Rnw:47-51 (eval = FALSE)
###################################################
## detach(addiction)
## addiction2 <- addiction[addiction$ill!=2,]
## attach(addiction2)
## age2 <- age^2


###################################################
### code chunk number 6: multinomial-addiction2.Rnw:56-59 (eval = FALSE)
###################################################
## m0vs1 <- glm(ill ~ as.factor(gender) + as.factor(university) + age + age2, 
## family=binomial())
## summary(m0vs1)


