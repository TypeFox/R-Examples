### R code from vignette source 'Ch_analysing_longitudinal_dataI.Rnw'
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
### code chunk number 3: ALDI-setup
###################################################
library("lme4")
library("multcomp")
residuals <- function(object) {
    y <- getME(object, 'y')
    y - fitted(object)
}


###################################################
### code chunk number 4: ALDI-plot-BtheB
###################################################
data("BtheB", package = "HSAUR3")
layout(matrix(1:2, nrow = 1))
ylim <- range(BtheB[,grep("bdi", names(BtheB))], 
              na.rm = TRUE)
tau <- subset(BtheB, treatment == "TAU")[,
    grep("bdi", names(BtheB))]
boxplot(tau, main = "Treated as Usual", ylab = "BDI", 
        xlab = "Time (in months)", names = c(0, 2, 3, 5, 8), 
        ylim = ylim)
btheb <- subset(BtheB, treatment == "BtheB")[,
    grep("bdi", names(BtheB))]
boxplot(btheb, main = "Beat the Blues", ylab = "BDI", 
        xlab = "Time (in months)", names = c(0, 2, 3, 5, 8), 
        ylim = ylim)


###################################################
### code chunk number 5: ALDI-long-BtheB
###################################################
data("BtheB", package = "HSAUR3")
BtheB$subject <- factor(rownames(BtheB))
nobs <- nrow(BtheB)
BtheB_long <- reshape(BtheB, idvar = "subject", 
    varying = c("bdi.2m", "bdi.3m", "bdi.5m", "bdi.8m"),
    direction = "long") 
BtheB_long$time <- rep(c(2, 3, 5, 8), rep(nobs, 4)) 


###################################################
### code chunk number 6: ALDI-showlong-BtheB
###################################################
subset(BtheB_long, subject %in% c("1", "2", "3"))


###################################################
### code chunk number 7: ALDI-fit-BtheB
###################################################
library("lme4")
BtheB_lmer1 <- lmer(bdi ~ bdi.pre + time + treatment + drug + 
    length + (1 | subject), data = BtheB_long, 
    REML = FALSE, na.action = na.omit)
BtheB_lmer2 <- lmer(bdi ~ bdi.pre + time + treatment + drug + 
    length + (time | subject), data = BtheB_long, 
    REML = FALSE, na.action = na.omit)   
anova(BtheB_lmer1, BtheB_lmer2)


###################################################
### code chunk number 8: ALDI-summary-BtheB
###################################################
summary(BtheB_lmer1)


###################################################
### code chunk number 9: ALDI-summary-BtheB-p
###################################################
cftest(BtheB_lmer1)


###################################################
### code chunk number 10: ALDI-qqnorm-BtheB
###################################################
layout(matrix(1:2, ncol = 2))
qint <- ranef(BtheB_lmer1)$subject[["(Intercept)"]]
qres <- residuals(BtheB_lmer1)
qqnorm(qint, ylab = "Estimated random intercepts",
       xlim = c(-3, 3), ylim = c(-20, 20),
       main = "Random intercepts")
qqline(qint)
qqnorm(qres, xlim = c(-3, 3), ylim = c(-20, 20),
       ylab = "Estimated residuals",
       main = "Residuals")
qqline(qres)


###################################################
### code chunk number 11: ALDI-dropout
###################################################
bdi <- BtheB[, grep("bdi", names(BtheB))]
plot(1:4, rep(-0.5, 4), type = "n", axes = FALSE, 
     ylim = c(0, 50), xlab = "Months", ylab = "BDI")
axis(1, at = 1:4, labels = c(0, 2, 3, 5))
axis(2)
for (i in 1:4) {
    dropout <- is.na(bdi[,i + 1])
    points(rep(i, nrow(bdi)) + ifelse(dropout, 0.05, -0.05), 
           jitter(bdi[,i]), pch = ifelse(dropout, 20, 1))
}


