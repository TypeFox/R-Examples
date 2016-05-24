### R code from vignette source 'Ch_logistic_regression_glm.Rnw'
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
### code chunk number 3: GLM-plasma-plot
###################################################
data("plasma", package = "HSAUR3")
layout(matrix(1:2, ncol = 2))
cdplot(ESR ~ fibrinogen, data = plasma)
cdplot(ESR ~ globulin, data = plasma)


###################################################
### code chunk number 4: GLM-plasma-fit1
###################################################
plasma_glm_1 <- glm(ESR ~ fibrinogen, data = plasma, 
                    family = binomial())


###################################################
### code chunk number 5: GLM-plasma-summary-1
###################################################
summary(plasma_glm_1)


###################################################
### code chunk number 6: GLM-plasma-confint
###################################################
ci <- confint(plasma_glm_1)["fibrinogen",]


###################################################
### code chunk number 7: GLM-plasma-confint
###################################################
confint(plasma_glm_1, parm = "fibrinogen")


###################################################
### code chunk number 8: GLM-plasma-confint
###################################################
print(ci) 


###################################################
### code chunk number 9: GLM-plasma-exp
###################################################
exp(coef(plasma_glm_1)["fibrinogen"])


###################################################
### code chunk number 10: GLM-plasma-exp-ci
###################################################
ci <- exp(confint(plasma_glm_1, parm = "fibrinogen"))


###################################################
### code chunk number 11: GLM-plasma-exp-ci
###################################################
exp(confint(plasma_glm_1, parm = "fibrinogen"))


###################################################
### code chunk number 12: GLM-plasma-exp-ci
###################################################
print(ci)


###################################################
### code chunk number 13: GLM-plasma-fit2
###################################################
plasma_glm_2 <- glm(ESR ~ fibrinogen +  globulin, 
    data = plasma, family = binomial())


###################################################
### code chunk number 14: GLM-plasma-summary-2
###################################################
summary(plasma_glm_2)


###################################################
### code chunk number 15: GLM-plasma-anova-hide
###################################################
plasma_anova <- anova(plasma_glm_1, plasma_glm_2, test = "Chisq")


###################################################
### code chunk number 16: GLM-plasma-anova
###################################################
anova(plasma_glm_1, plasma_glm_2, test = "Chisq")


###################################################
### code chunk number 17: GLM-plasma-predict
###################################################
prob <- predict(plasma_glm_2, type = "response")


###################################################
### code chunk number 18: GLM-plasma-bubble
###################################################
plot(globulin ~ fibrinogen, data = plasma, xlim = c(2, 6), 
     ylim = c(25, 55), pch = ".")
symbols(plasma$fibrinogen, plasma$globulin, circles = prob,
        add = TRUE)


###################################################
### code chunk number 19: GLM-womensrole-fit1
###################################################
data("womensrole", package = "HSAUR3")
fm1 <- cbind(agree, disagree) ~ gender + education
womensrole_glm_1 <- glm(fm1, data = womensrole, 
                        family = binomial())


###################################################
### code chunk number 20: GLM-womensrole-summary-1
###################################################
summary(womensrole_glm_1)


###################################################
### code chunk number 21: GLM-womensrole-probfit
###################################################
role.fitted1 <- predict(womensrole_glm_1, type = "response")


###################################################
### code chunk number 22: GLM-plot-setup
###################################################
myplot <- function(role.fitted)  {
    f <- womensrole$gender == "Female"
    plot(womensrole$education, role.fitted, type = "n", 
         ylab = "Probability of agreeing",
         xlab = "Education", ylim = c(0,1))
    lines(womensrole$education[!f], role.fitted[!f], lty = 1)
    lines(womensrole$education[f], role.fitted[f], lty = 2)  
    lgtxt <- c("Fitted (Males)", "Fitted (Females)")
    legend("topright", lgtxt, lty = 1:2, bty = "n")
    y <-  womensrole$agree / (womensrole$agree + 
                              womensrole$disagree)
    size <- womensrole$agree + womensrole$disagree
    size <- size - min(size)
    size <- (size / max(size)) * 3 + 1
    text(womensrole$education, y, ifelse(f, "\\VE", "\\MA"),
         family = "HersheySerif", cex = size)         
} 


###################################################
### code chunk number 23: GLM-role-fitted1
###################################################
myplot(role.fitted1)


###################################################
### code chunk number 24: GLM-womensrole-fit2
###################################################
fm2 <- cbind(agree,disagree) ~ gender * education
womensrole_glm_2 <- glm(fm2, data = womensrole, 
                        family = binomial())


###################################################
### code chunk number 25: GLM-womensrole-summary-2
###################################################
summary(womensrole_glm_2)


###################################################
### code chunk number 26: GLM-role-fitted2
###################################################
role.fitted2 <- predict(womensrole_glm_2, type = "response")
myplot(role.fitted2)


###################################################
### code chunk number 27: GLM-role-plot2
###################################################
res <- residuals(womensrole_glm_2, type = "deviance")
plot(predict(womensrole_glm_2), res,
     xlab="Fitted values", ylab = "Residuals", 
     ylim = max(abs(res)) * c(-1,1))
abline(h = 0, lty = 2)


###################################################
### code chunk number 28: GLM-polyps-fit1
###################################################
data("polyps", package = "HSAUR3")
polyps_glm_1 <- glm(number ~ treat + age, data = polyps, 
                    family = poisson())


###################################################
### code chunk number 29: GLM-polyps-summary-1
###################################################
summary(polyps_glm_1)


###################################################
### code chunk number 30: GLM-polyp-quasi
###################################################
polyps_glm_2 <- glm(number ~ treat + age, data = polyps, 
                    family = quasipoisson())
summary(polyps_glm_2)


###################################################
### code chunk number 31: GLM-backpain-clogit
###################################################
library("survival")
backpain_glm <- clogit(I(status == "case") ~ 
    driver + suburban + strata(ID), data = backpain)


###################################################
### code chunk number 32: GLM-backpain-print
###################################################
print(backpain_glm)


###################################################
### code chunk number 33: GLM-CHFLS-polr
###################################################
library("MASS")
opts <- options(contrasts = c("contr.treatment", 
                              "contr.helmert"))
CHFLS_polr <- polr(R_happy ~ ., data = CHFLS, Hess = TRUE)
options(opts)


###################################################
### code chunk number 34: GLM-CHFLS-polr
###################################################
summary(CHFLS_polr)


###################################################
### code chunk number 35: GLM-CHFLS-polr-helmert
###################################################
H <- with(CHFLS, contr.helmert(table(R_health)))
rownames(H) <- levels(CHFLS$R_health)
colnames(H) <- paste(levels(CHFLS$R_health)[-1], "- avg")
H


###################################################
### code chunk number 36: GLM-CHFLS-polr-cftest
###################################################
library("multcomp")
op <- options(digits = 2)
cf <- cftest(CHFLS_polr)
cftest <- function(x, digits = max(3, getOption("digits") - 3)) {
    x <- cf
    cat("\n\t", "Simultaneous Tests for General Linear Hypotheses\n\n")
    if (!is.null(x$type)) 
        cat("Multiple Comparisons of Means:", x$type, "Contrasts\n\n\n")
    call <- if (isS4(x$model)) 
        x$model@call
    else x$model$call
    if (!is.null(call)) {
        cat("Fit: ")
        print(call)
        cat("\n")
    }
    pq <- x$test
    mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
    error <- attr(pq$pvalues, "error")
    pname <- switch(x$alternativ, less = paste("Pr(<", ifelse(x$df == 
        0, "z", "t"), ")", sep = ""), greater = paste("Pr(>", 
        ifelse(x$df == 0, "z", "t"), ")", sep = ""), two.sided = paste("Pr(>|", 
        ifelse(x$df == 0, "z", "t"), "|)", sep = ""))
    colnames(mtests) <- c("Estimate", "Std. Error", ifelse(x$df == 
        0, "z value", "t value"), pname)
    type <- pq$type
    if (!is.null(error) && error > .Machine$double.eps) {
        sig <- which.min(abs(1/error - (10^(1:10))))
        sig <- 1/(10^sig)
    }
    else {
        sig <- .Machine$double.eps
    }
    cat("Linear Hypotheses:\n")
    alt <- switch(x$alternative, two.sided = "==", less = ">=", 
        greater = "<=")
    rownames(mtests) <- rownames(mtests)
    printCoefmat(mtests, digits = digits, has.Pvalue = TRUE, 
        P.values = TRUE, eps.Pvalue = sig)
    switch(type, univariate = cat("(Univariate p values reported)"), 
        `single-step` = cat("(Adjusted p values reported -- single-step method)"), 
        Shaffer = cat("(Adjusted p values reported -- Shaffer method)"), 
        Westfall = cat("(Adjusted p values reported -- Westfall method)"), 
        cat("(Adjusted p values reported --", type, "method)"))
    cat("\n\n")
    invisible(x)
}


###################################################
### code chunk number 37: GLM-CHFLS-polr-cftest
###################################################
library("multcomp")
cftest(CHFLS_polr)


###################################################
### code chunk number 38: GLM-CHFLS-polr-cftest
###################################################
options(op)


###################################################
### code chunk number 39: GLM-CHFLS-pred-1
###################################################
CHFLS[1,]


###################################################
### code chunk number 40: GLM-CHFLS-pred-2
###################################################
nd <- CHFLS[rep(1, nlevels(CHFLS$R_health)),]
nd$R_health <- ordered(levels(nd$R_health), 
                       labels = levels(nd$R_health))


###################################################
### code chunk number 41: GLM-CHFLS-pred-3
###################################################
(dens <- predict(CHFLS_polr, newdata = nd, type = "probs"))


###################################################
### code chunk number 42: GLM-CHFLS-pred-plot
###################################################
library("lattice")
D <- expand.grid(R_health = nd$R_health, 
                 R_happy = ordered(LETTERS[1:4]))
D$dens <- as.vector(dens)
barchart(dens ~ R_happy | R_health, data = D, 
         ylab = "Density", xlab = "Happiness",)


###################################################
### code chunk number 43: GLM-findings
###################################################
ci <- round(exp(confint(plasma_glm_1, parm = "fibrinogen")), 2)
ci <- paste("(", paste(ci, collapse = ","), ")", sep = "")


