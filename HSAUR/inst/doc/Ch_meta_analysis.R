### R code from vignette source 'Ch_meta_analysis.Rnw'
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
### code chunk number 2: MA-smoking-OR
###################################################
library("rmeta")
data("smoking", package = "HSAUR")
smokingOR <- meta.MH(smoking[["tt"]], smoking[["tc"]],
                     smoking[["qt"]], smoking[["qc"]], 
                     names = rownames(smoking))


###################################################
### code chunk number 3: MA-smoking-OR-summary
###################################################
summary(smokingOR)


###################################################
### code chunk number 4: MA-smoking-OR-plot
###################################################
plot(smokingOR, ylab = "")


###################################################
### code chunk number 5: MA-smoking-random
###################################################
smokingDSL <- meta.DSL(smoking[["tt"]], smoking[["tc"]], 
                     smoking[["qt"]], smoking[["qc"]],
                     names = rownames(smoking))
print(smokingDSL)


###################################################
### code chunk number 6: MA-BCG-odds
###################################################
data("BCG", package = "HSAUR")
BCG_OR <- meta.MH(BCG[["BCGVacc"]], BCG[["NoVacc"]],
                  BCG[["BCGTB"]], BCG[["NoVaccTB"]],
                  names = BCG$Study)
BCG_DSL <- meta.DSL(BCG[["BCGVacc"]], BCG[["NoVacc"]],
                  BCG[["BCGTB"]], BCG[["NoVaccTB"]],
                  names = BCG$Study)


###################################################
### code chunk number 7: MA-BCGOR-summary
###################################################
summary(BCG_OR)


###################################################
### code chunk number 8: MA-BCGDSL-summary
###################################################
summary(BCG_DSL)


###################################################
### code chunk number 9: BCG-studyweights
###################################################
studyweights <- 1 / (BCG_DSL$tau2 + BCG_DSL$selogs^2)
y <- BCG_DSL$logs
BCG_mod <- lm(y ~ Latitude + Year, data = BCG, 
              weights = studyweights)


###################################################
### code chunk number 10: MA-mod-summary
###################################################
summary(BCG_mod)


###################################################
### code chunk number 11: BCG-Latitude-plot
###################################################
plot(y ~ Latitude, data = BCG, ylab = "Estimated log-OR")
abline(lm(y ~ Latitude, data = BCG, weights = studyweights))


###################################################
### code chunk number 12: MA-funnel-ex
###################################################
set.seed(290875)
sigma <- seq(from = 1/10, to = 1, length.out = 35)
y <- rnorm(35) * sigma
gr <- (y > -0.5)
layout(matrix(1:2, ncol = 1))
plot(y, 1/sigma, xlab = "Effect size", ylab = "1 / standard error")
plot(y[gr], 1/(sigma[gr]), xlim = range(y),
     xlab = "Effect size", ylab = "1 / standard error")


###################################################
### code chunk number 13: MA-smoking-funnel
###################################################
funnelplot(smokingDSL$logs, smokingDSL$selogs, 
           summ = smokingDSL$logDSL, xlim = c(-1.7, 1.7))
abline(v = 0, lty = 2)


