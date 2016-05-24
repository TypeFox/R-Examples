### R code from vignette source 'Ch_simultaneous_inference.Rnw'
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
### code chunk number 3: SIMC-setup
###################################################
library("multcomp")
library("coin")
library("sandwich")
library("lme4")


###################################################
### code chunk number 4: SIMC-alpha-data-figure
###################################################
n <- table(alpha$alength)
levels(alpha$alength) <- abbreviate(levels(alpha$alength), 4)
plot(elevel ~ alength, data = alpha, varwidth = TRUE,
     ylab = "Expression Level", 
     xlab = "NACP-REP1 Allele Length")
axis(3, at = 1:3, labels = paste("n = ", n))


###################################################
### code chunk number 5: SIMC-alpha-aov-tukey
###################################################
library("multcomp")
amod <- aov(elevel ~ alength, data = alpha)
amod_glht <- glht(amod, linfct = mcp(alength = "Tukey"))


###################################################
### code chunk number 6: SIMC-alpha-aov-tukey-K
###################################################
amod_glht$linfct


###################################################
### code chunk number 7: SIMC-alpha-aov-coefvcov
###################################################
coef(amod_glht)
vcov(amod_glht)


###################################################
### code chunk number 8: SIMC-alpha-aov-results
###################################################
confint(amod_glht)
summary(amod_glht)


###################################################
### code chunk number 9: SIMC-aov-tukey-sandwich
###################################################
amod_glht_sw <- glht(amod, linfct = mcp(alength = "Tukey"), 
                     vcov = sandwich)
summary(amod_glht_sw)


###################################################
### code chunk number 10: SIMC-alpha-confint-plot
###################################################
par(mai = par("mai") * c(1, 2.1, 1, 0.5))
layout(matrix(1:2, ncol = 2))
ci1 <- confint(glht(amod, linfct = mcp(alength = "Tukey")))
ci2 <- confint(glht(amod, linfct = mcp(alength = "Tukey"), 
               vcov = sandwich))
ox <- expression(paste("Tukey (ordinary ", bold(S)[n], ")"))
sx <- expression(paste("Tukey (sandwich ", bold(S)[n], ")"))
plot(ci1, xlim = c(-0.6, 2.6), main = ox, 
    xlab = "Difference", ylim = c(0.5, 3.5))
plot(ci2, xlim = c(-0.6, 2.6), main = sx, 
    xlab = "Difference", ylim = c(0.5, 3.5))


###################################################
### code chunk number 11: SIMC-trees-setup
###################################################
trees513 <- subset(trees513, !species %in% c("fir", "ash/maple/elm/lime", "softwood (other)"))
trees513$species <- trees513$species[,drop = TRUE]
levels(trees513$species)[nlevels(trees513$species)] <- "hardwood"


###################################################
### code chunk number 12: SIMC-trees-lmer
###################################################
mmod <- glmer(damage ~ species - 1 + (1 | lattice / plot),
              data = trees513, family = binomial())
K <- diag(length(fixef(mmod)))
K


###################################################
### code chunk number 13: SIMC-trees-K
###################################################
colnames(K) <- rownames(K) <- 
    paste(gsub("species", "", names(fixef(mmod))),
          " (", table(trees513$species), ")", sep = "")
K


###################################################
### code chunk number 14: SIMC-trees-ci
###################################################
ci <- confint(glht(mmod, linfct = K))
ci$confint <- 1 - binomial()$linkinv(ci$confint)
ci$confint[,2:3] <- ci$confint[,3:2]


###################################################
### code chunk number 15: SIMC-trees-plot
###################################################
getOption("SweaveHooks")[["bigleftpar"]]()
plot(ci, xlab = "Probability of Damage Caused by Browsing", 
     xlim = c(0, 0.5), main = "", ylim = c(0.5, 5.5))


###################################################
### code chunk number 16: SIMC-clouds-confband
###################################################
confband <- function(subset, main) {
    mod <- lm(rainfall ~ sne, data = clouds, subset = subset)
    sne_grid <- seq(from = 1.5, to = 4.5, by = 0.25)
    K <- cbind(1, sne_grid)
    sne_ci <- confint(glht(mod, linfct = K))
    plot(rainfall ~ sne, data = clouds, subset = subset,
         xlab = "S-Ne criterion", main = main, 
         xlim = range(clouds$sne), 
         ylim = range(clouds$rainfall))
    abline(mod)
    lines(sne_grid, sne_ci$confint[,2], lty = 2)
    lines(sne_grid, sne_ci$confint[,3], lty = 2)
}


###################################################
### code chunk number 17: SIMC-clouds-lmplot
###################################################
layout(matrix(1:2, ncol = 2))
confband(clouds$seeding == "no", main = "No seeding")
confband(clouds$seeding == "yes", main = "Seeding")


