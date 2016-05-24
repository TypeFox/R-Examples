### R code from vignette source 'Ch_gam.Rnw'
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
### code chunk number 3: packages
###################################################
library("mgcv")
library("mboost")
library("rpart")
library("wordcloud")


###################################################
### code chunk number 4: GAM-men1500m-plot
###################################################
plot(time ~ year, data = men1500m, xlab = "Year",
     ylab = "Winning time (sec)")


###################################################
### code chunk number 5: GAM-men1500m-lm
###################################################
men1500m1900 <- subset(men1500m, year >= 1900)
men1500m_lm <- lm(time ~ year, data = men1500m1900)
plot(time ~ year, data = men1500m1900, xlab = "Year",
     ylab = "Winning time (sec)")
abline(men1500m_lm)


###################################################
### code chunk number 6: GAM-men1500m-smooth
###################################################
x <- men1500m1900$year
y <- men1500m1900$time
men1500m_lowess <- lowess(x, y)
plot(time ~ year, data = men1500m1900, xlab = "Year",
     ylab = "Winning time (sec)")
lines(men1500m_lowess, lty = 2)
men1500m_cubic <- gam(y ~ s(x, bs = "cr"))
lines(x, predict(men1500m_cubic), lty = 3)


###################################################
### code chunk number 7: GAM-men1500m-quad
###################################################
men1500m_lm2 <- lm(time ~ year + I(year^2), 
                   data = men1500m1900)
plot(time ~ year, data = men1500m1900, xlab = "Year",
     ylab = "Winning time (sec)")
lines(men1500m1900$year, predict(men1500m_lm2))


###################################################
### code chunk number 8: GAM-men1500m-pred
###################################################
predict(men1500m_lm, 
        newdata = data.frame(year = c(2008, 2012)), 
        interval = "confidence")
predict(men1500m_lm2, 
        newdata = data.frame(year = c(2008, 2012)), 
        interval = "confidence")


###################################################
### code chunk number 9: GAM-USairpollution-boost
###################################################
library("mboost")
USair_boost <- gamboost(SO2 ~ ., data = USairpollution)
USair_aic <- AIC(USair_boost)
USair_aic


###################################################
### code chunk number 10: GAM-USairpollution-boostplot
###################################################
USair_gam <- USair_boost[mstop(USair_aic)]
layout(matrix(1:6, ncol = 3))
plot(USair_gam, ask = FALSE)


###################################################
### code chunk number 11: GAM-USairpollution-residplot
###################################################
SO2hat <- predict(USair_gam)
SO2 <- USairpollution$SO2
plot(SO2hat, SO2 - SO2hat, type = "n",         
     xlim = c(-20, max(SO2hat) * 1.1),
     ylim = range(SO2 - SO2hat) * c(2, 1))
textplot(SO2hat, SO2 - SO2hat, rownames(USairpollution), 
         show.lines = FALSE, new = FALSE)
abline(h = 0, lty = 2, col = "grey")


###################################################
### code chunk number 12: GAM-kyphosis-plot
###################################################
layout(matrix(1:3, nrow = 1))
spineplot(Kyphosis ~ Age, data = kyphosis, 
          ylevels = c("present", "absent"))
spineplot(Kyphosis ~ Number, data = kyphosis, 
          ylevels = c("present", "absent"))
spineplot(Kyphosis ~ Start, data = kyphosis, 
         ylevels = c("present", "absent"))


###################################################
### code chunk number 13: GAM-kyphosis-gam
###################################################
(kyphosis_gam <- gam(Kyphosis ~ s(Age, bs = "cr") + 
    s(Number, bs = "cr", k = 3) + s(Start, bs = "cr", k = 3),
    family = binomial, data = kyphosis))


###################################################
### code chunk number 14: GAM-kyphosis-gamplot
###################################################
trans <- function(x) 
    binomial()$linkinv(x)
layout(matrix(1:3, nrow = 1))
plot(kyphosis_gam, select = 1, shade = TRUE, trans = trans)
plot(kyphosis_gam, select = 2, shade = TRUE, trans = trans)
plot(kyphosis_gam, select = 3, shade = TRUE, trans = trans)


###################################################
### code chunk number 15: GAM-womensrole-gam
###################################################
data("womensrole", package = "HSAUR3")
fm1 <- cbind(agree, disagree) ~ s(education, by = gender)
womensrole_gam <- gam(fm1, data = womensrole, 
                      family = binomial())


###################################################
### code chunk number 16: GAM-womensrole-gamplot
###################################################
layout(matrix(1:2, nrow = 1))
plot(womensrole_gam, select = 1, shade = TRUE)
plot(womensrole_gam, select = 1, shade = TRUE)


###################################################
### code chunk number 17: GAM-plot-setup
###################################################
myplot <- function(role.fitted)  {    
    f <- womensrole$gender == "Female"
    plot(womensrole$education, role.fitted, type = "n",
         ylab = "Probability of agreeing",
         xlab = "Education", ylim = c(0,1))
    lines(womensrole$education[!f], role.fitted[!f], lty = 1)
    lines(womensrole$education[f], role.fitted[f], lty = 2)
    lgtxt <- c("Fitted (Males)", "Fitted (Females)")
    legend("topright", lgtxt, lty = 1:2, bty = "n")
    y <-  womensrole$agree / (womensrole$agree +
                              womensrole$disagree)
    size <- womensrole$agree + womensrole$disagree
    size <- size - min(size)
    size <- (size / max(size)) * 3 + 1
    text(womensrole$education, y, ifelse(f, "\\VE", "\\MA"),
         family = "HersheySerif", cex = size)
} 


###################################################
### code chunk number 18: GAM-womensrole-probplot
###################################################
myplot(predict(womensrole_gam, type = "response"))


