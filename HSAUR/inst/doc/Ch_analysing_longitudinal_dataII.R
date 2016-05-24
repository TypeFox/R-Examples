### R code from vignette source 'Ch_analysing_longitudinal_dataII.Rnw'
### Encoding: ASCII

###################################################
### code chunk number 1: setup
###################################################
rm(list = ls())
s <- search()[-1]
s <- s[-match(c("package:base", "package:stats", "package:graphics", "package:grDevices",
                "package:utils", "package:datasets", "package:methods", "Autoloads"), s)]
if (length(s) > 0) sapply(s, detach, character.only = TRUE)
if (!file.exists("tables")) dir.create("tables")
if (!file.exists("figures")) dir.create("figures")
set.seed(290875)
options(prompt = "R> ", continue = "+  ",
    width = 63, # digits = 4, 
    SweaveHooks = list(leftpar = function() 
        par(mai = par("mai") * c(1, 1.05, 1, 1))))
HSAURpkg <- require("HSAUR")
if (!HSAURpkg) stop("cannot load package ", sQuote("HSAUR"))
rm(HSAURpkg)
a <- Sys.setlocale("LC_ALL", "C")
book <- TRUE
refs <- cbind(c("AItR", "SI", "CI", "ANOVA", "MLR", "GLM", 
                "DE", "RP", "SA", "ALDI", "ALDII", "MA", "PCA", 
                "MDS", "CA"), 1:15)
ch <- function(x, book = TRUE) {
    ch <- refs[which(refs[,1] == x),]
    if (book) {
        return(paste("Chapter~\\\\ref{", ch[1], "}", sep = ""))
    } else {
        return(paste("Chapter~\\\\ref{", ch[2], "}", sep = ""))
    }
}


###################################################
### code chunk number 2: ALDII-gee
###################################################
library("gee")


###################################################
### code chunk number 3: ALDII-BtheB-data
###################################################
data("BtheB", package = "HSAUR")
BtheB$subject <- factor(rownames(BtheB))
nobs <- nrow(BtheB)
BtheB_long <- reshape(BtheB, idvar = "subject",
    varying = c("bdi.2m", "bdi.4m", "bdi.6m", "bdi.8m"),
    direction = "long")
BtheB_long$time <- rep(c(2, 4, 6, 8), rep(nobs, 4))


###################################################
### code chunk number 4: ALDII-BtheB-geefit-indep
###################################################
osub <- order(as.integer(BtheB_long$subject))
BtheB_long <- BtheB_long[osub,]
btb_gee <- gee(bdi ~ bdi.pre + treatment + length + drug, 
    data = BtheB_long, id = subject, family = gaussian,
    corstr = "independence")


###################################################
### code chunk number 5: ALDII-BtheB-geefit-ex
###################################################
btb_gee1 <- gee(bdi ~ bdi.pre + treatment + length + drug, 
    data = BtheB_long, id = subject, family = gaussian,
    corstr = "exchangeable")


###################################################
### code chunk number 6: ALDII-BtheB-geesummary
###################################################
summary(btb_gee)


###################################################
### code chunk number 7: ALDII-BtheB-gee1summary
###################################################
summary(btb_gee1)


###################################################
### code chunk number 8: ALDII-respiratory-data
###################################################
data("respiratory", package = "HSAUR")
resp <- subset(respiratory, month > "0")
resp$baseline <- rep(subset(respiratory, month == "0")$status, rep(4, 111))
resp$nstat <- as.numeric(resp$status == "good")


###################################################
### code chunk number 9: ALDII-respiratory-fit
###################################################
resp_glm <- glm(status ~ centre + treatment + sex + baseline + 
    age, data = resp, family = "binomial")
resp_gee1 <- gee(nstat ~ centre + treatment + sex + baseline + 
    age, data = resp, family = "binomial", id = subject, 
    corstr = "independence", scale.fix = TRUE, scale.value = 1)
resp_gee2 <- gee(nstat ~ centre + treatment + sex + baseline + 
    age, data = resp, family = "binomial", id = subject, 
    corstr = "exchangeable", scale.fix = TRUE, scale.value = 1)


###################################################
### code chunk number 10: ALDII-resp-glm-summary
###################################################
summary(resp_glm)


###################################################
### code chunk number 11: ALDII-resp-gee1summary
###################################################
summary(resp_gee1)


###################################################
### code chunk number 12: ALDII-resp-gee2-summary
###################################################
summary(resp_gee2)


###################################################
### code chunk number 13: ALDII-resp-confint
###################################################
se <- summary(resp_gee2)$coefficients["treatmenttreatment",
                                      "Robust S.E."]
coef(resp_gee2)["treatmenttreatment"] +  
    c(-1, 1) * se * qnorm(0.975)


###################################################
### code chunk number 14: ALDII-resp-confint-exp
###################################################
exp(coef(resp_gee2)["treatmenttreatment"] + 
    c(-1, 1) * se * qnorm(0.975))


###################################################
### code chunk number 15: ALDII-epilepsy
###################################################
data("epilepsy", package = "HSAUR")
itp <- interaction(epilepsy$treatment, epilepsy$period)
tapply(epilepsy$seizure.rate, itp, mean)
tapply(epilepsy$seizure.rate, itp, var)


###################################################
### code chunk number 16: ALDII-plot1
###################################################
layout(matrix(1:2, nrow = 1))
ylim <- range(epilepsy$seizure.rate)
placebo <- subset(epilepsy, treatment == "placebo")
progabide <- subset(epilepsy, treatment == "Progabide")
boxplot(seizure.rate ~ period, data = placebo,
        ylab = "Number of seizures", 
        xlab = "Period", ylim = ylim, main = "Placebo")
boxplot(seizure.rate ~ period, data = progabide,
        main  = "Progabide", ylab = "Number of seizures", 
        xlab = "Period", ylim = ylim)


###################################################
### code chunk number 17: ALDII-plot2
###################################################
layout(matrix(1:2, nrow = 1))
ylim <- range(log(epilepsy$seizure.rate + 1))
boxplot(log(seizure.rate + 1) ~ period, data = placebo,
        main = "Placebo", ylab = "Log number of seizures", 
        xlab = "Period", ylim = ylim)
boxplot(log(seizure.rate + 1) ~ period, data = progabide,
        main = "Progabide", ylab = "Log number of seizures", 
        xlab = "Period", ylim = ylim)


###################################################
### code chunk number 18: ALDII-epilepsy-gee
###################################################
per <- rep(log(2),nrow(epilepsy))
epilepsy$period <- as.numeric(epilepsy$period)
fm <- seizure.rate ~ base + age + treatment + offset(per)
epilepsy_glm <- glm(fm, data = epilepsy, family = "poisson")
epilepsy_gee1 <- gee(fm, data = epilepsy, family = "poisson", 
    id = subject, corstr = "independence", scale.fix = TRUE, 
    scale.value = 1)
epilepsy_gee2 <- gee(fm, data = epilepsy, family = "poisson", 
    id = subject, corstr = "exchangeable", scale.fix = TRUE, 
    scale.value = 1)
epilepsy_gee3 <- gee(fm, data = epilepsy, family = "poisson", 
    id = subject, corstr = "exchangeable", scale.fix = FALSE, 
    scale.value = 1)


###################################################
### code chunk number 19: ALDII-espilepsy-glm-summary
###################################################
summary(epilepsy_glm)


###################################################
### code chunk number 20: ALDII-espilepsy-gee1-summary
###################################################
summary(epilepsy_gee1)


###################################################
### code chunk number 21: ALDII-espilepsy-gee2-summary
###################################################
summary(epilepsy_gee2)


###################################################
### code chunk number 22: ALDII-espilepsy-gee3-summary
###################################################
summary(epilepsy_gee3)


