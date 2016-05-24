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
    show.signif.stars = FALSE,
    SweaveHooks = list(leftpar = function() 
        par(mai = par("mai") * c(1, 1.05, 1, 1)),
        bigleftpar = function()
        par(mai = par("mai") * c(1, 1.7, 1, 1))))
HSAURpkg <- require("HSAUR3")
if (!HSAURpkg) stop("cannot load package ", sQuote("HSAUR2"))
rm(HSAURpkg)
 ### </FIXME> hm, R-2.4.0 --vanilla seems to need this
a <- Sys.setlocale("LC_ALL", "C")
 ### </FIXME>
book <- TRUE
refs <- cbind(c("AItR", "DAGD", "SI", "CI", "ANOVA", "MLR", "GLM", 
                "DE", "RP", "GAM", "SA", "ALDI", "ALDII", "SIMC", "MA", "PCA", 
                "MDS", "CA"), 1:18)
ch <- function(x) {
    ch <- refs[which(refs[,1] == x),]
    if (book) {
        return(paste("Chapter~\\\\ref{", ch[1], "}", sep = ""))
    } else {
        return(paste("Chapter~", ch[2], sep = ""))
    }
}
if (file.exists("deparse.R"))
    source("deparse.R")

setHook(packageEvent("lattice", "attach"), function(...) {
    lattice.options(default.theme = 
        function()
            standard.theme("pdf", color = FALSE))
    })


###################################################
### code chunk number 2: singlebook
###################################################
book <- FALSE


###################################################
### code chunk number 3: setup
###################################################
options(digits = 3)

if (!interactive()) {
print.summary.gee <- function (x, digits = NULL, quote = FALSE, prefix = "", ...) 
{
    if (is.null(digits)) 
        digits <- options()$digits
    else options(digits = digits)
    cat("...")
    cat("\nModel:\n")
    cat(" Link:                     ", x$model$link, "\n")
    cat(" Variance to Mean Relation:", x$model$varfun, "\n")
    if (!is.null(x$model$M)) 
        cat(" Correlation Structure:    ", x$model$corstr, ", M =", 
            x$model$M, "\n")
    else cat(" Correlation Structure:    ", x$model$corstr, "\n")
    cat("\n...")
    nas <- x$nas
    if (!is.null(nas) && any(nas)) 
        cat("\n\nCoefficients: (", sum(nas), " not defined because of singularities)\n", 
            sep = "")
    else cat("\n\nCoefficients:\n")
    print(x$coefficients, digits = digits)
    cat("\nEstimated Scale Parameter: ", format(round(x$scale, 
        digits)))
    cat("\n...\n")

    invisible(x)
}
}



###################################################
### code chunk number 4: ALDII-gee
###################################################
library("gee")


###################################################
### code chunk number 5: ALDII-BtheB-data
###################################################
data("BtheB", package = "HSAUR3")
BtheB$subject <- factor(rownames(BtheB))
nobs <- nrow(BtheB)
BtheB_long <- reshape(BtheB, idvar = "subject",
    varying = c("bdi.2m", "bdi.3m", "bdi.5m", "bdi.8m"),
    direction = "long")
BtheB_long$time <- rep(c(2, 3, 5, 8), rep(nobs, 4))
names(BtheB_long)[names(BtheB_long) == "treatment"] <- "trt"


###################################################
### code chunk number 6: ALDII-BtheB-geefit-indep
###################################################
osub <- order(as.integer(BtheB_long$subject))
BtheB_long <- BtheB_long[osub,]
btb_gee <- gee(bdi ~ bdi.pre + trt + length + drug, 
    data = BtheB_long, id = subject, family = gaussian,
    corstr = "independence")


###################################################
### code chunk number 7: ALDII-BtheB-geefit-ex
###################################################
btb_gee1 <- gee(bdi ~ bdi.pre + trt + length + drug, 
    data = BtheB_long, id = subject, family = gaussian,
    corstr = "exchangeable")


###################################################
### code chunk number 8: ALDII-BtheB-geesummary
###################################################
summary(btb_gee)


###################################################
### code chunk number 9: ALDII-BtheB-gee1summary
###################################################
summary(btb_gee1)


###################################################
### code chunk number 10: ALDII-respiratory-data
###################################################
data("respiratory", package = "HSAUR3")
resp <- subset(respiratory, month > "0")
resp$baseline <- rep(subset(respiratory, month == "0")$status,
                     rep(4, 111))
resp$nstat <- as.numeric(resp$status == "good")
resp$month <- resp$month[, drop = TRUE]


###################################################
### code chunk number 11: ALDII-respiratory-names
###################################################
names(resp)[names(resp) == "treatment"] <- "trt"
levels(resp$trt)[2] <- "trt"


###################################################
### code chunk number 12: ALDII-respiratory-fit
###################################################
resp_glm <- glm(status ~ centre + trt + gender + baseline 
    + age, data = resp, family = "binomial")
resp_gee1 <- gee(nstat ~ centre + trt + gender + baseline 
    + age, data = resp, family = "binomial", id = subject, 
    corstr = "independence", scale.fix = TRUE, 
    scale.value = 1)
resp_gee2 <- gee(nstat ~ centre + trt + gender + baseline 
    + age, data = resp, family = "binomial", id = subject, 
    corstr = "exchangeable", scale.fix = TRUE, 
    scale.value = 1)


###################################################
### code chunk number 13: ALDII-resp-glm-summary
###################################################
summary(resp_glm)


###################################################
### code chunk number 14: ALDII-resp-gee1summary
###################################################
summary(resp_gee1)


###################################################
### code chunk number 15: ALDII-resp-gee2-summary
###################################################
summary(resp_gee2)


###################################################
### code chunk number 16: ALDII-resp-confint
###################################################
se <- summary(resp_gee2)$coefficients["trttrt",
                                      "Robust S.E."]
coef(resp_gee2)["trttrt"] +  
    c(-1, 1) * se * qnorm(0.975)


###################################################
### code chunk number 17: ALDII-resp-confint-exp
###################################################
exp(coef(resp_gee2)["trttrt"] + 
    c(-1, 1) * se * qnorm(0.975))


###################################################
### code chunk number 18: ALDII-epilepsy
###################################################
data("epilepsy", package = "HSAUR3")
itp <- interaction(epilepsy$treatment, epilepsy$period)
tapply(epilepsy$seizure.rate, itp, mean)
tapply(epilepsy$seizure.rate, itp, var)


###################################################
### code chunk number 19: ALDII-plot1
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
### code chunk number 20: ALDII-plot2
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
### code chunk number 21: ALDII-epilepsy-gee
###################################################
per <- rep(log(2),nrow(epilepsy))
epilepsy$period <- as.numeric(epilepsy$period)
names(epilepsy)[names(epilepsy) == "treatment"] <- "trt"
fm <- seizure.rate ~ base + age + trt + offset(per)
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
### code chunk number 22: ALDII-espilepsy-glm-summary
###################################################
summary(epilepsy_glm)


###################################################
### code chunk number 23: ALDII-espilepsy-gee1-summary
###################################################
summary(epilepsy_gee1)


###################################################
### code chunk number 24: ALDII-espilepsy-gee2-summary
###################################################
summary(epilepsy_gee2)


###################################################
### code chunk number 25: ALDII-espilepsy-gee3-summary
###################################################
summary(epilepsy_gee3)


###################################################
### code chunk number 26: ALDII-respiratory-lmer
###################################################
library("lme4")
resp_lmer <- glmer(status ~ baseline + month + 
    trt + gender + age + centre + (1 | subject), 
    family = binomial(), data = resp)
exp(fixef(resp_lmer))


###################################################
### code chunk number 27: ALDII-resp-lmer-dirty
###################################################
su <- summary(resp_lmer)
if (!interactive()) {
    summary <- function(x) {
        cat("\n...\n")
        cat("Fixed effects:\n")
        lme4V <- packageDescription("lme4")$Version
        if (compareVersion("0.999999-2", lme4V) >= 0) {
            printCoefmat(su@coefs)
        } else {
            printCoefmat(su$coefficients)
        }
        cat("\n...\n")
    }
}


###################################################
### code chunk number 28: ALDII-resp-lmer-summary
###################################################
summary(resp_lmer)


