### R code from vignette source 'Ch_simple_inference.Rnw'
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
### code chunk number 2: SI-setup
###################################################
library("vcd")


###################################################
### code chunk number 3: SI-roomwidth-data
###################################################
data("roomwidth", package = "HSAUR")


###################################################
### code chunk number 4: SI-roomwidth-convert
###################################################
convert <- ifelse(roomwidth$unit == "feet", 1, 3.28)


###################################################
### code chunk number 5: SI-roomwidth-summary
###################################################
tapply(roomwidth$width * convert, roomwidth$unit, summary)
tapply(roomwidth$width * convert, roomwidth$unit, sd)


###################################################
### code chunk number 6: SI-roomwidth-boxplot
###################################################
layout(matrix(c(1,2,1,3), nrow = 2, ncol = 2, byrow = FALSE))
boxplot(I(width * convert) ~ unit, data = roomwidth,
        ylab = "Estimated width (feet)",
        varwidth = TRUE, names = c("Estimates in feet",
        "Estimates in metres (converted to feet)"))
feet <- roomwidth$unit == "feet"
qqnorm(roomwidth$width[feet], 
       ylab = "Estimated width (feet)")
qqline(roomwidth$width[feet])
qqnorm(roomwidth$width[!feet], 
       ylab = "Estimated width (metres)")
qqline(roomwidth$width[!feet])


###################################################
### code chunk number 7: SI-roomwidth-formula
###################################################
I(width * convert) ~ unit


###################################################
### code chunk number 8: SI-roomwidth-tt-T-hide
###################################################
tt <- t.test(I(width * convert) ~ unit, data = roomwidth, 
             var.equal = TRUE)


###################################################
### code chunk number 9: SI-roomwidth-tt-T
###################################################
t.test(I(width * convert) ~ unit, data = roomwidth, 
       var.equal = TRUE)


###################################################
### code chunk number 10: SI-roomwidth-tt-F
###################################################
t.test(I(width * convert) ~ unit, data = roomwidth, 
       var.equal = FALSE)


###################################################
### code chunk number 11: SI-roomwidth-wt
###################################################
wilcox.test(I(width * convert) ~ unit, data = roomwidth, 
            conf.int = TRUE)


###################################################
### code chunk number 12: SI-roomwidth-wt-hide
###################################################
pwt <- round(wilcox.test(I(width * convert) ~ unit, data =
             roomwidth)$p.value, 3)


###################################################
### code chunk number 13: SI-waves-data
###################################################
data("waves", package = "HSAUR")


###################################################
### code chunk number 14: SI-wavese-boxplot
###################################################
mooringdiff <- waves$method1 - waves$method2
layout(matrix(1:2, ncol = 2))
boxplot(mooringdiff, ylab = "Differences (Newton metres)", 
        main = "Boxplot")
abline(h = 0, lty = 2)
qqnorm(mooringdiff, ylab = "Differences (Newton metres)")
qqline(mooringdiff)


###################################################
### code chunk number 15: SI-waves-tt
###################################################
t.test(mooringdiff)


###################################################
### code chunk number 16: SI-waves-wt
###################################################
pwt <- round(wilcox.test(mooringdiff)$p.value, 3)


###################################################
### code chunk number 17: SI-waves-wt
###################################################
wilcox.test(mooringdiff)


###################################################
### code chunk number 18: SI-water-data
###################################################
data("water", package = "HSAUR")


###################################################
### code chunk number 19: SI-water-plot
###################################################
nf <- layout(matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE),
             c(2, 1), c(1, 2), TRUE)
psymb <- as.numeric(water$location)
plot(mortality ~ hardness, data = water, pch = psymb)
abline(lm(mortality ~ hardness, data = water))
legend("topright", legend = levels(water$location), 
       pch = c(1,2), bty = "n")
hist(water$hardness)
boxplot(water$mortality)


###################################################
### code chunk number 20: SI-water-cor
###################################################
cor.test(~ mortality + hardness, data = water)


###################################################
### code chunk number 21: SI-water-cor
###################################################
cr <- round(cor.test(~ mortality + hardness, data = water)$estimate, 3)


###################################################
### code chunk number 22: SI-pistonrings-chisq-hide
###################################################
chisqt <- chisq.test(pistonrings)


###################################################
### code chunk number 23: SI-pistonrings-chisq
###################################################
data("pistonrings", package = "HSAUR")
chisq.test(pistonrings)


###################################################
### code chunk number 24: SI-pistonrings-resid
###################################################
chisq.test(pistonrings)$residuals


###################################################
### code chunk number 25: SI-assoc-plot
###################################################
library("vcd")
assoc(pistonrings)


###################################################
### code chunk number 26: SI-rearrests-data
###################################################
data("rearrests", package = "HSAUR")
rearrests


###################################################
### code chunk number 27: SI-rearrests-mcnemar
###################################################
mcs <- round(mcnemar.test(rearrests, correct = FALSE)$statistic, 3)


###################################################
### code chunk number 28: SI-arrests-mcnemar
###################################################
mcnemar.test(rearrests, correct = FALSE)


###################################################
### code chunk number 29: SI-arrests-binom
###################################################
binom.test(rearrests[2], n = sum(rearrests[c(2,3)]))


