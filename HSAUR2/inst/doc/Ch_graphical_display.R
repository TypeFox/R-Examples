### R code from vignette source 'Ch_graphical_display.Rnw'
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
### code chunk number 3: DAGD-USmelanoma-histbox
###################################################
xr <- range(USmelanoma$mortality) * c(0.9, 1.1)
xr


###################################################
### code chunk number 4: DAGD-USmelanoma-histbox
###################################################
layout(matrix(1:2, nrow = 2))
par(mar = par("mar") * c(0.8, 1, 1, 1))
boxplot(USmelanoma$mortality, ylim = xr, horizontal = TRUE, 
        xlab = "Mortality")
hist(USmelanoma$mortality, xlim = xr, xlab = "", main = "", 
     axes = FALSE, ylab = "")
axis(1)


###################################################
### code chunk number 5: DAGD-USmelanoma-boxocean
###################################################
plot(mortality ~ ocean, data = USmelanoma, 
     xlab = "Contiguity to an ocean", ylab = "Mortality")


###################################################
### code chunk number 6: DAGD-USmelanoma-dens
###################################################
dyes <- with(USmelanoma, density(mortality[ocean == "yes"]))
dno <- with(USmelanoma, density(mortality[ocean == "no"]))
plot(dyes, lty = 1, xlim = xr, main = "", ylim = c(0, 0.018))
lines(dno, lty = 2)
legend("topleft", lty = 1:2, legend = c("Coastal State", 
       "Land State"), bty = "n")


###################################################
### code chunk number 7: DAGD-USmelanoma-xy
###################################################
layout(matrix(1:2, ncol = 2))
plot(mortality ~ longitude, data = USmelanoma)
plot(mortality ~ latitude, data = USmelanoma)


###################################################
### code chunk number 8: DAGD-USmelanoma-lat
###################################################
plot(mortality ~ latitude, data = USmelanoma, 
     pch = as.integer(USmelanoma$ocean))
legend("topright", legend = c("Land state", "Coast state"), 
       pch = 1:2, bty = "n")


###################################################
### code chunk number 9: DAGD-USmelanoma-south
###################################################
subset(USmelanoma, latitude < 32)


###################################################
### code chunk number 10: DAGD-CHFLS-happy
###################################################
barplot(xtabs(~ R_happy, data = CHFLS))


###################################################
### code chunk number 11: DAGD-CHFLS-health_happy_xtabs
###################################################
xtabs(~ R_happy + R_health, data = CHFLS)


###################################################
### code chunk number 12: DAGD-CHFLS-health_happy_xtabs2
###################################################
hh <- xtabs(~ R_health + R_happy, data = CHFLS)


###################################################
### code chunk number 13: DAGD-CHFLS-health_happy
###################################################
plot(R_happy ~ R_health, data = CHFLS)


###################################################
### code chunk number 14: DAGD-CHFLS-happy_income
###################################################
layout(matrix(1:2, ncol = 2))
plot(R_happy ~ log(R_income + 1), data = CHFLS)
cdplot(R_happy ~ log(R_income + 1), data = CHFLS)


###################################################
### code chunk number 15: DAGD-CHFLS-RAincome3 (eval = FALSE)
###################################################
## xyplot(jitter(log(A_income + 0.5)) ~ 
##        jitter(log(R_income + 0.5)) | R_edu, data = CHFLS)


###################################################
### code chunk number 16: DAGD-CHFLS-RAincome3
###################################################
library("lattice")
lattice.options(default.theme =
        function()
            standard.theme("pdf", color = FALSE))
print(xyplot(jitter(log(A_income + 0.5)) ~ jitter(log(R_income + 0.5)) | R_edu, data = CHFLS))


###################################################
### code chunk number 17: DAGD-household-tab
###################################################
data("household", package = "HSAUR2")
toLatex(HSAURtable(household), 
  caption = paste("Household expenditure for single men and women."),
  label = "DAGD-household-tab")


###################################################
### code chunk number 18: DAGD-USstates-tab
###################################################
data("USstates", package = "HSAUR2")
toLatex(HSAURtable(USstates), 
  caption = paste("Socio-demographic variables for ten US states."),
  label = "DAGD-USstates-tab")


###################################################
### code chunk number 19: DAGD-suicides2-tab
###################################################
data("suicides2", package = "HSAUR2")
toLatex(HSAURtable(suicides2), 
  caption = paste("Mortality rates per $100,000$ from male suicides."),
  label = "DAGD-suicides2-tab", rownames = TRUE)


###################################################
### code chunk number 20: DAGD-banknote-tab
###################################################
data("banknote", package = "alr3")
banknote$Y <- NULL
banknote <- banknote[c(1:5, 101:200),]
toLatex(HSAURtable(banknote, pkg = "alr3", nrow = 10),
  caption = paste("Swiss bank note data."),
  label = "DAGD-banknote-tab", rownames = FALSE)


