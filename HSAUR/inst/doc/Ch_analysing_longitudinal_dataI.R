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
### code chunk number 2: ALDI-setup
###################################################
library("Matrix")
library("lme4")


###################################################
### code chunk number 3: ALDI-plot-BtheB
###################################################
data("BtheB", package = "HSAUR")
layout(matrix(1:2, nrow = 1))
ylim <- range(BtheB[,grep("bdi", names(BtheB))], 
              na.rm = TRUE)
tau <- subset(BtheB, treatment == "TAU")[,
    grep("bdi", names(BtheB))]
boxplot(tau, main = "Treated as usual", ylab = "BDI", 
        xlab = "Time (in months)", names = c(0, 2, 4, 6, 8), 
        ylim = ylim)
btheb <- subset(BtheB, treatment == "BtheB")[,
    grep("bdi", names(BtheB))]
boxplot(btheb, main = "Beat the Blues", ylab = "BDI", 
        xlab = "Time (in months)", names = c(0, 2, 4, 6, 8), 
        ylim = ylim)


###################################################
### code chunk number 4: ALDI-long-BtheB
###################################################
data("BtheB", package = "HSAUR")
BtheB$subject <- factor(rownames(BtheB))
nobs <- nrow(BtheB)
BtheB_long <- reshape(BtheB, idvar = "subject", 
    varying = c("bdi.2m", "bdi.4m", "bdi.6m", "bdi.8m"),
    direction = "long") 
BtheB_long$time <- rep(c(2, 4, 6, 8), rep(nobs, 4)) 


###################################################
### code chunk number 5: ALDI-showlong-BtheB
###################################################
subset(BtheB_long, subject %in% c("1", "2", "3"))


###################################################
### code chunk number 6: ALDI-fit-BtheB
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
### code chunk number 7: ALDI-summary-BtheB
###################################################
summary(BtheB_lmer1)


