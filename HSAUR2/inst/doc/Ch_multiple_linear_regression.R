### R code from vignette source 'Ch_multiple_linear_regression.Rnw'
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
HSAURpkg <- require("HSAUR2")
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
### code chunk number 3: MLR-setup
###################################################
library("wordcloud")


###################################################
### code chunk number 4: MLR-hubble-tab
###################################################
data("hubble", package = "gamair")
names(hubble) <- c("galaxy", "velocity", "distance")
toLatex(HSAURtable(hubble, package = "gamair"), pcol = 2,
    caption = paste("Distance and velocity for 24 galaxies."),
    label = "MLR-hubble-tab")


###################################################
### code chunk number 5: MLR-clouds-tab
###################################################
data("clouds", package = "HSAUR2")
toLatex(HSAURtable(clouds), pcol = 1, 
    caption = paste("Cloud seeding experiments in Florida -- see text for",
                    "explanations of the variables."),
    label = "MLR-clouds-tab")


###################################################
### code chunk number 6: MLR-hubble-plot
###################################################
plot(velocity ~ distance, data = hubble)


###################################################
### code chunk number 7: MLR-hubble-beta1
###################################################
sum(hubble$distance * hubble$velocity) / 
    sum(hubble$distance^2)


###################################################
### code chunk number 8: MLR-hubble-lm
###################################################
hmod <- lm(velocity ~ distance - 1, data = hubble)


###################################################
### code chunk number 9: MLR-hubble-lm
###################################################
coef(hmod)


###################################################
### code chunk number 10: MLR-hubble-age
###################################################
Mpc <- 3.09 * 10^19
ysec <- 60^2 * 24 * 365.25
Mpcyear <- Mpc / ysec
1 / (coef(hmod) / Mpcyear)


###################################################
### code chunk number 11: MLR-hubble-lmplot
###################################################
layout(matrix(1:2, ncol = 2))
plot(velocity ~ distance, data = hubble)
abline(hmod)
plot(hmod, which = 1)


###################################################
### code chunk number 12: MLR-clouds-boxplots
###################################################
data("clouds", package = "HSAUR2")
layout(matrix(1:2, nrow = 2))
bxpseeding <- boxplot(rainfall ~ seeding, data = clouds, 
    ylab = "Rainfall", xlab = "Seeding")
bxpecho <- boxplot(rainfall ~ echomotion, data = clouds, 
    ylab = "Rainfall", xlab = "Echo Motion")


###################################################
### code chunk number 13: MLR-clouds-scatterplots
###################################################
layout(matrix(1:4, nrow = 2))
plot(rainfall ~ time, data = clouds)
plot(rainfall ~ cloudcover, data = clouds)
plot(rainfall ~ sne, data = clouds, xlab="S-Ne criterion")
plot(rainfall ~ prewetness, data = clouds)


###################################################
### code chunk number 14: MLR-clouds-outliers
###################################################
rownames(clouds)[clouds$rainfall %in% c(bxpseeding$out, 
                                        bxpecho$out)]


###################################################
### code chunk number 15: MLR-clouds-formula
###################################################
clouds_formula <- rainfall ~ seeding + 
    seeding:(sne + cloudcover + prewetness + echomotion) + 
    time


###################################################
### code chunk number 16: MLR-clouds-modelmatrix
###################################################
Xstar <- model.matrix(clouds_formula, data = clouds)


###################################################
### code chunk number 17: MLR-clouds-contrasts
###################################################
attr(Xstar, "contrasts")


###################################################
### code chunk number 18: MLR-clouds-lm
###################################################
clouds_lm <- lm(clouds_formula, data = clouds)
class(clouds_lm)


###################################################
### code chunk number 19: MLR-clouds-summary
###################################################
summary(clouds_lm)


###################################################
### code chunk number 20: MLR-clouds-coef
###################################################
betastar <- coef(clouds_lm)
betastar


###################################################
### code chunk number 21: MLR-clouds-vcov
###################################################
Vbetastar <- vcov(clouds_lm)


###################################################
### code chunk number 22: MLR-clouds-sd
###################################################
sqrt(diag(Vbetastar))


###################################################
### code chunk number 23: MLR-clouds-lmplot
###################################################
psymb <- as.numeric(clouds$seeding)
plot(rainfall ~ sne, data = clouds, pch = psymb, 
     xlab = "S-Ne criterion")
abline(lm(rainfall ~ sne, data = clouds, 
          subset = seeding == "no"))
abline(lm(rainfall ~ sne, data = clouds, 
          subset = seeding == "yes"), lty = 2)  
legend("topright", legend = c("No seeding", "Seeding"), 
       pch = 1:2, lty = 1:2, bty = "n")


###################################################
### code chunk number 24: MLR-clouds-residfitted
###################################################
clouds_resid <- residuals(clouds_lm)
clouds_fitted <- fitted(clouds_lm)


###################################################
### code chunk number 25: MLR-clouds-residplot
###################################################
plot(clouds_fitted, clouds_resid, xlab = "Fitted values", 
     ylab = "Residuals", type = "n", 
     ylim = max(abs(clouds_resid)) * c(-1, 1))
abline(h = 0, lty = 2)
textplot(clouds_fitted, clouds_resid, words = rownames(clouds), new = FALSE)


###################################################
### code chunk number 26: MLR-clouds-qqplot
###################################################
qqnorm(clouds_resid, ylab = "Residuals")
qqline(clouds_resid)


###################################################
### code chunk number 27: MLR-clouds-cook (eval = FALSE)
###################################################
## plot(clouds_lm)


###################################################
### code chunk number 28: MLR-clouds-cook
###################################################
plot(clouds_lm, which = 4, sub.caption = NULL)


