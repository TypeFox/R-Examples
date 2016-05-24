### R code from vignette source 'Ch_analysis_of_variance.Rnw'
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
### code chunk number 2: ANOVA-weightgain-mean-var
###################################################
data("weightgain", package = "HSAUR")
tapply(weightgain$weightgain, 
       list(weightgain$source, weightgain$type), mean)
tapply(weightgain$weightgain, 
       list(weightgain$source, weightgain$type), sd)


###################################################
### code chunk number 3: ANOVA-weightgain-plot
###################################################
plot.design(weightgain)


###################################################
### code chunk number 4: ANOVA-weightgain-aov
###################################################
wg_aov <- aov(weightgain ~ source * type, data = weightgain)


###################################################
### code chunk number 5: ANOVA-weightgain-aov-summary
###################################################
summary(wg_aov)


###################################################
### code chunk number 6: ANOVA-weightgain-iplot (eval = FALSE)
###################################################
## interaction.plot(weightgain$type, weightgain$source, 
##                  weightgain$weightgain) 


###################################################
### code chunk number 7: ANOVA-weightgain-iplot-nice
###################################################
interaction.plot(weightgain$type, weightgain$source, weightgain$weightgain,
legend = FALSE)
legend(1.5, 95, legend = levels(weightgain$source), title = "weightgain$source", 
       lty = c(2,1), bty = "n")


###################################################
### code chunk number 8: ANOVA-weightgain-coef
###################################################
coef(wg_aov)


###################################################
### code chunk number 9: ANOVA-weightgain-contrasts
###################################################
options("contrasts")


###################################################
### code chunk number 10: ANOVA-weightgain-coef-sum
###################################################
coef(aov(weightgain ~ source + type + source:type, 
    data = weightgain, contrasts = list(source = contr.sum)))


###################################################
### code chunk number 11: ANOVA-foster
###################################################
data("foster", package = "HSAUR")


###################################################
### code chunk number 12: ANOVA-foster-plot
###################################################
plot.design(foster)


###################################################
### code chunk number 13: ANOVA-foster-aov-one (eval = FALSE)
###################################################
## summary(aov(weight ~ litgen * motgen, data = foster))


###################################################
### code chunk number 14: ANOVA-foster-aov-one
###################################################
summary(aov(weight ~ litgen * motgen, data = foster))


###################################################
### code chunk number 15: ANOVA-foster-aov-two (eval = FALSE)
###################################################
## summary(aov(weight ~ motgen * litgen, data = foster))


###################################################
### code chunk number 16: ANOVA-foster-aov-two
###################################################
summary(aov(weight ~ motgen * litgen, data = foster))


###################################################
### code chunk number 17: ANOVA-weightgain-again (eval = FALSE)
###################################################
## summary(aov(weightgain ~ type * source, data = weightgain))


###################################################
### code chunk number 18: ANOVA-foster-aov
###################################################
foster_aov <- aov(weight ~ litgen * motgen, data = foster)


###################################################
### code chunk number 19: ANOVA-foster-tukeyHSD
###################################################
foster_hsd <- TukeyHSD(foster_aov, "motgen")
foster_hsd


###################################################
### code chunk number 20: ANOVA-foster-tukeyHSDplot
###################################################
plot(foster_hsd)


###################################################
### code chunk number 21: ANOVA-water-manova
###################################################
data("water", package = "HSAUR")
summary(manova(cbind(hardness, mortality) ~ location, 
    data = water), test = "Hotelling-Lawley")


###################################################
### code chunk number 22: ANOVA-water-means
###################################################
tapply(water$hardness, water$location, mean)
tapply(water$mortality, water$location, mean)


###################################################
### code chunk number 23: ANOVA-skulls-data
###################################################
data("skulls", package = "HSAUR")
means <- aggregate(skulls[,c("mb", "bh", "bl", "nh")], 
                   list(epoch = skulls$epoch), mean)
means


###################################################
### code chunk number 24: ANOVA-skulls-fig
###################################################
pairs(means[,-1],
    panel = function(x, y) {
        text(x, y, abbreviate(levels(skulls$epoch)))
    })


###################################################
### code chunk number 25: ANOVA-skulls-manova
###################################################
skulls_manova <- manova(cbind(mb, bh, bl, nh) ~ epoch, 
                        data = skulls)
summary(skulls_manova, test = "Pillai")
summary(skulls_manova, test = "Wilks") 
summary(skulls_manova, test = "Hotelling-Lawley")
summary(skulls_manova, test = "Roy")


###################################################
### code chunk number 26: ANOVA-skulls-manova2
###################################################
summary.aov(skulls_manova)


###################################################
### code chunk number 27: ANOVA-skulls-manova3
###################################################
summary(manova(cbind(mb, bh, bl, nh) ~ epoch, data = skulls, 
               subset = epoch %in% c("c4000BC", "c3300BC")))
summary(manova(cbind(mb, bh, bl, nh) ~ epoch, data = skulls, 
               subset = epoch %in% c("c4000BC", "c1850BC")))
summary(manova(cbind(mb, bh, bl, nh) ~ epoch, data = skulls, 
               subset = epoch %in% c("c4000BC", "c200BC")))
summary(manova(cbind(mb, bh, bl, nh) ~ epoch, data = skulls, 
               subset = epoch %in% c("c4000BC", "cAD150")))


