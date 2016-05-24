### R code from vignette source 'Ch_cluster_analysis.Rnw'
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
### code chunk number 3: thissetup
###################################################
library("mclust")
library("mvtnorm")
mai <- par("mai")
options(SweaveHooks = list(rmai = function() { par(mai = mai * c(1,1,1,2))}))
data("pottery", package = "HSAUR2")


###################################################
### code chunk number 4: CA-pottery-dist
###################################################
pottery_dist <- dist(pottery[, colnames(pottery) != "kiln"])
library("lattice")
levelplot(as.matrix(pottery_dist), xlab = "Pot Number", 
          ylab = "Pot Number")


###################################################
### code chunk number 5: CA-pottery-distplot
###################################################
trellis.par.set(standard.theme(color = FALSE))
plot(levelplot(as.matrix(pottery_dist), xlab = "Pot Number", ylab = "Pot Number"))


###################################################
### code chunk number 6: CA-pottery-hclust
###################################################
pottery_single <- hclust(pottery_dist, method = "single")
pottery_complete <- hclust(pottery_dist, method = "complete")
pottery_average <- hclust(pottery_dist, method = "average")
layout(matrix(1:3, ncol = 3))
plot(pottery_single, main = "Single Linkage", 
     sub = "", xlab = "")
plot(pottery_complete, main = "Complete Linkage", 
     sub = "", xlab = "")
plot(pottery_average, main = "Average Linkage", 
     sub = "", xlab = "")


###################################################
### code chunk number 7: pottery-cluster
###################################################
pottery_cluster <- cutree(pottery_average, h = 4)
xtabs(~ pottery_cluster + kiln, data = pottery)


###################################################
### code chunk number 8: CA-planets-scatter
###################################################
getOption("SweaveHooks")[["rmai"]]()
data("planets", package = "HSAUR2")
library("scatterplot3d")
scatterplot3d(log(planets$mass), log(planets$period), 
    log(planets$eccen + ifelse(planets$eccen == 0, 0.001, 0)), 
    type = "h", angle = 55, pch = 16, 
    y.ticklabs = seq(0, 10, by = 2), 
    y.margin.add = 0.1, scale.y = 0.7,
    xlab = "log(mass)", ylab = "log(period)", 
    zlab = "log(eccen)")


###################################################
### code chunk number 9: CA-planet-ss
###################################################
rge <- apply(planets, 2, max) - apply(planets, 2, min)
planet.dat <- sweep(planets, 2, rge, FUN = "/")
n <- nrow(planet.dat)
wss <- rep(0, 10)
wss[1] <- (n - 1) * sum(apply(planet.dat, 2, var))
for (i in 2:10)
    wss[i] <- sum(kmeans(planet.dat, 
                         centers = i)$withinss)
plot(1:10, wss, type = "b", xlab = "Number of groups",
     ylab = "Within groups sum of squares")


###################################################
### code chunk number 10: CA-planets-kmeans3
###################################################
planet_kmeans3 <- kmeans(planet.dat, centers = 3)
table(planet_kmeans3$cluster)


###################################################
### code chunk number 11: CA-planets-ccent
###################################################
ccent <- function(cl) {
    f <- function(i) colMeans(planets[cl == i,])
    x <- sapply(sort(unique(cl)), f)
    colnames(x) <- sort(unique(cl))
    return(x)
}


###################################################
### code chunk number 12: CA-planets--kmeans3-ccent
###################################################
ccent(planet_kmeans3$cluster)


###################################################
### code chunk number 13: CA-planets-kmeans5
###################################################
planet_kmeans5 <- kmeans(planet.dat, centers = 5)
table(planet_kmeans5$cluster)
ccent(planet_kmeans5$cluster)


###################################################
### code chunk number 14: CA-planets-mclust
###################################################
library("mclust")
planet_mclust <- Mclust(planet.dat)


###################################################
### code chunk number 15: CA-planets-mclust-plot
###################################################
plot(planet_mclust, planet.dat, what = "BIC", col = "black", 
     ylab = "-BIC", ylim = c(0, 350))  


###################################################
### code chunk number 16: CA-planets-mclust-print
###################################################
print(planet_mclust)


###################################################
### code chunk number 17: CA-planets-mclust-scatter
###################################################
clPairs(planet.dat, 
    classification = planet_mclust$classification, 
    symbols = 1:3, col = "black")


###################################################
### code chunk number 18: CA-planets-mclust-scatterclust
###################################################
getOption("SweaveHooks")[["rmai"]]()
scatterplot3d(log(planets$mass), log(planets$period), 
    log(planets$eccen + ifelse(planets$eccen == 0, 0.001, 0)), 
    type = "h", angle = 55,  scale.y = 0.7, 
    pch = planet_mclust$classification, 
    y.ticklabs = seq(0, 10, by = 2), y.margin.add = 0.1,
    xlab = "log(mass)", ylab = "log(period)", 
    zlab = "log(eccen)")


###################################################
### code chunk number 19: CA-planets-mclust-mu
###################################################
table(planet_mclust$classification)
ccent(planet_mclust$classification)


