### R code from vignette source 'binary-unemployment.Rnw'

###################################################
### code chunk number 1: binary-unemployment.Rnw:11-12 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: binary-unemployment.Rnw:17-22 (eval = FALSE)
###################################################
## unemployment <- matrix(c(403, 238, 167, 175), nrow=2, ncol=2)
## rownames(unemployment) <- c("male","female")
## colnames(unemployment) <- c("<6 month",">6 month")
## unemployment
## rowSums(unemployment)


###################################################
### code chunk number 3: binary-unemployment.Rnw:25-29 (eval = FALSE)
###################################################
## ( odds_m <- 403/167 )
## ( odds_w <- 238/175 )
## ( log_odds_m <- log(403/167) )
## ( log_odds_w <- log(238/175) )


###################################################
### code chunk number 4: binary-unemployment.Rnw:33-35 (eval = FALSE)
###################################################
## gender <- c(rep(1, 403+167), rep(0,238+175))
## unemp <- c(rep(1, 403), rep(0, 167), rep(1, 238), rep(0, 175))


###################################################
### code chunk number 5: binary-unemployment.Rnw:38-39 (eval = FALSE)
###################################################
## table(gender, unemp)


###################################################
### code chunk number 6: binary-unemployment.Rnw:42-46 (eval = FALSE)
###################################################
## bin <- glm(unemp ~ gender, family=binomial)
## summary(bin)
## bin$coef
## exp(bin$coef)


###################################################
### code chunk number 7: binary-unemployment.Rnw:49-50 (eval = FALSE)
###################################################
## gender_effect <- c(rep(1, 403+167), rep(-1,238+175))


###################################################
### code chunk number 8: binary-unemployment.Rnw:53-54 (eval = FALSE)
###################################################
## table(gender_effect, unemp)


###################################################
### code chunk number 9: binary-unemployment.Rnw:57-61 (eval = FALSE)
###################################################
## bin_effect <- glm(unemp ~ gender_effect, family=binomial)
## summary(bin_effect)
## bin_effect$coef
## exp(bin_effect$coef)


###################################################
### code chunk number 10: binary-unemployment.Rnw:65-70 (eval = FALSE)
###################################################
## unemp_level <- matrix(c(202, 307, 87, 45,
##                         96, 162, 66, 18), nrow=4, ncol=2)
## colnames(unemp_level) <- c("Short term","Long term")
## unemp_level
## rowSums(unemp_level)


###################################################
### code chunk number 11: binary-unemployment.Rnw:74-77 (eval = FALSE)
###################################################
## level <- factor(c(rep(1, 202+96), rep(2,307+162), rep(3,87+66), rep(4,45+18))) 
## unemp_l <-  c(rep(1, 202), rep(0, 96), rep(1, 307), rep(0, 162),
##             rep(1, 87), rep(0, 66), rep(1, 45), rep(0, 18))


###################################################
### code chunk number 12: binary-unemployment.Rnw:80-81 (eval = FALSE)
###################################################
## table(level, unemp_l)


###################################################
### code chunk number 13: binary-unemployment.Rnw:85-88 (eval = FALSE)
###################################################
## level <- relevel(level, ref=4)
## bin_l <- glm(unemp_l ~ level, family=binomial)
## summary(bin_l)


###################################################
### code chunk number 14: binary-unemployment.Rnw:93-97 (eval = FALSE)
###################################################
## library(qvcalc)
## qv<-qvcalc(bin_l,"level")
## summary(qv)
## plot(qv)


###################################################
### code chunk number 15: binary-unemployment.Rnw:100-101 (eval = FALSE)
###################################################
## rm(unemployment)


