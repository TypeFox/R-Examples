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
### code chunk number 2: MLR-clouds-boxplots
###################################################
data("clouds", package = "HSAUR")
layout(matrix(1:2, nrow = 2))
bxpseeding <- boxplot(rainfall ~ seeding, data = clouds, 
    ylab = "Rainfall", xlab = "Seeding")
bxpecho <- boxplot(rainfall ~ echomotion, data = clouds, 
    ylab = "Rainfall", xlab = "Echo Motion")


###################################################
### code chunk number 3: MLR-clouds-scatterplots
###################################################
layout(matrix(1:4, nrow = 2))
plot(rainfall ~ time, data = clouds)
plot(rainfall ~ cloudcover, data = clouds)
plot(rainfall ~ sne, data = clouds, xlab="S-Ne criterion")
plot(rainfall ~ prewetness, data = clouds)


###################################################
### code chunk number 4: MLR-clouds-outliers
###################################################
rownames(clouds)[clouds$rainfall %in% c(bxpseeding$out, 
                                        bxpecho$out)]


###################################################
### code chunk number 5: MLR-clouds-formula
###################################################
clouds_formula <- rainfall ~ seeding * (sne + cloudcover +
    prewetness + echomotion) + time


###################################################
### code chunk number 6: MLR-clouds-modelmatrix
###################################################
Xstar <- model.matrix(clouds_formula, data = clouds)


###################################################
### code chunk number 7: MLR-clouds-contrasts
###################################################
attr(Xstar, "contrasts")


###################################################
### code chunk number 8: MLR-clouds-lm
###################################################
clouds_lm <- lm(clouds_formula, data = clouds)
class(clouds_lm)


###################################################
### code chunk number 9: MLR-clouds-summary
###################################################
summary(clouds_lm)


###################################################
### code chunk number 10: MLR-clouds-coef
###################################################
betastar <- coef(clouds_lm)
betastar


###################################################
### code chunk number 11: MLR-clouds-vcov
###################################################
Vbetastar <- vcov(clouds_lm)


###################################################
### code chunk number 12: MLR-clouds-sd
###################################################
sqrt(diag(Vbetastar))


###################################################
### code chunk number 13: MLR-clouds-lmplot
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
### code chunk number 14: MLR-clouds-residfitted
###################################################
clouds_resid <- residuals(clouds_lm)
clouds_fitted <- fitted(clouds_lm)


###################################################
### code chunk number 15: MLR-clouds-residplot
###################################################
plot(clouds_fitted, clouds_resid, xlab = "Fitted values", 
     ylab = "Residuals", type = "n", 
     ylim = max(abs(clouds_resid)) * c(-1, 1))
abline(h = 0, lty = 2)
text(clouds_fitted, clouds_resid, labels = rownames(clouds))


###################################################
### code chunk number 16: MLR-clouds-qqplot
###################################################
qqnorm(clouds_resid, ylab = "Residuals")
qqline(clouds_resid)


###################################################
### code chunk number 17: MLR-clouds-cook (eval = FALSE)
###################################################
## plot(clouds_lm)


###################################################
### code chunk number 18: MLR-clouds-cook
###################################################
plot(clouds_lm, which = 4, sub.caption = NULL)


