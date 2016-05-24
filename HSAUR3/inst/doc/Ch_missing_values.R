### R code from vignette source 'Ch_missing_values.Rnw'
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
### code chunk number 3: MV-bp-tab
###################################################
data("bp", package = "HSAUR3")
toLatex(HSAURtable(bp), pcol = 2,
    caption = paste("Blood pressure data."),
    label = "MV-bp-tab")


###################################################
### code chunk number 4: MV-bp-NA
###################################################
sapply(bp, function(x) sum(is.na(x)))


###################################################
### code chunk number 5: MV-bp-msd-cc
###################################################
summary(bp$recovtime, na.rm = TRUE)


###################################################
### code chunk number 6: MV-bp-sd-cc
###################################################
sd(bp$recovtime, na.rm = TRUE)


###################################################
### code chunk number 7: MV-bp-cor-cc
###################################################
with(bp, cor(bloodp, recovtime, use = "complete.obs"))
with(bp, cor(logdose, recovtime, use = "complete.obs"))


###################################################
### code chunk number 8: MV-bp-pairs-cc
###################################################
layout(matrix(1:3, nrow = 1))
plot(bloodp ~ logdose, data = bp)
plot(recovtime ~ bloodp, data = bp)
plot(recovtime ~ logdose, data = bp)


###################################################
### code chunk number 9: MV-bp-lm-cc
###################################################
summary(lm(recovtime ~ bloodp + logdose, data = bp))


###################################################
### code chunk number 10: MV-bp-mice-pkg
###################################################
library("mice")


###################################################
### code chunk number 11: MV-bp-mice
###################################################
imp <- mice(bp, method = "mean", m = 1, maxit = 1)


###################################################
### code chunk number 12: MV-bp-imp-summary
###################################################
with(imp, summary(recovtime))


###################################################
### code chunk number 13: MV-bp-imp-sd
###################################################
with(imp, sd(recovtime))


###################################################
### code chunk number 14: MV-bp-imp-cor
###################################################
with(imp, cor(bloodp, recovtime))
with(imp, cor(logdose, recovtime)) 


###################################################
### code chunk number 15: MV-bp-pairs-imp
###################################################
layout(matrix(1:2, nrow = 1))
plot(recovtime ~ bloodp, data = complete(imp), 
     pch = is.na(bp$recovtime) + 1)
plot(recovtime ~ logdose, data = complete(imp), 
     pch = is.na(bp$recovtime) + 1)
legend("topleft", pch = 1:2, bty = "n",
       legend = c("original", "imputed"))


###################################################
### code chunk number 16: MV-bp-lm-imp
###################################################
with(imp, summary(lm(recovtime ~ bloodp + logdose)))


###################################################
### code chunk number 17: MV-bp-mice
###################################################
imp_ppm <- mice(bp, m = 10, method = "pmm", 
                print = FALSE, seed = 1)


###################################################
### code chunk number 18: MV-bp-pairs-mice
###################################################
layout(matrix(1:2, nrow = 1))
plot(recovtime ~ bloodp, data = complete(imp_ppm),      
     pch = is.na(bp$recovtime) + 1)
plot(recovtime ~ logdose, data = complete(imp_ppm),
     pch = is.na(bp$recovtime) + 1)
legend("topleft", pch = 1:2, bty = "n",
       legend = c("original", "imputed"))


###################################################
### code chunk number 19: MV-bp-mice-out
###################################################
summary(unlist(with(imp_ppm, mean(recovtime))$analyses))
summary(unlist(with(imp_ppm, sd(recovtime))$analyses))


###################################################
### code chunk number 20: MV-bp-mice-cor
###################################################
summary(unlist(with(imp_ppm, 
        cor(bloodp, recovtime))$analyses))
summary(unlist(with(imp_ppm, 
        cor(logdose, recovtime))$analyses))


###################################################
### code chunk number 21: MV-bp-mice-lm
###################################################
fit <- with(imp_ppm, lm(recovtime ~ bloodp + logdose))


###################################################
### code chunk number 22: MV-bp-lm-mice
###################################################
summary(pool(fit))[, c("est", "se", "t", "Pr(>|t|)")]


###################################################
### code chunk number 23: MI-bp-t
###################################################
with(bp, t.test(recovtime, mu = 27))
with(imp, t.test(recovtime, mu = 27))$analyses[[1]]


###################################################
### code chunk number 24: MI-mice-t
###################################################
fit <- with(imp_ppm, lm(I(recovtime - 27) ~ 1))
summary(pool(fit))[, c("est", "se", "t", "Pr(>|t|)")]


###################################################
### code chunk number 25: MI-UStemp-tab
###################################################
data("UStemp", package = "HSAUR3")
toLatex(HSAURtable(UStemp), 
    caption = "Lowest temperatures in Fahrenheit recorded in various months for cities in the US.",
    label = "MI-UStemp-tab", rownames = TRUE)


