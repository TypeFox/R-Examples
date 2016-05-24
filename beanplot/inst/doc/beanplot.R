### R code from vignette source 'beanplot.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("beanplot")
options(prompt = "R> ", continue = "+  ")


###################################################
### code chunk number 2: rug
###################################################
set.seed(1)
par(mfrow = c(1, 2), mai = c(0.8, 0.5, 0.5, 0.5))
plot(density(d <- rnorm(50)), main = "plot(density(d <- rnorm(50))); rug(d)")
rug(d)
beanplot(d, bw = "nrd0", horizontal = TRUE, main = "beanplot(d, bw = \"nrd0\", horizontal = TRUE)")


###################################################
### code chunk number 3: distributions
###################################################
library("beanplot")
par(mfrow = c(1, 2), mai = c(0.5, 0.5, 0.5, 0.1))
mu <- 2
si <- 0.6
c <- 500
bimodal <- c(rnorm(c/2, -mu, si), rnorm(c/2, mu, si))
uniform <- runif(c, -4, 4)
normal <- rnorm(c, 0, 1.5)
ylim <- c(-7, 7)
boxplot(bimodal, uniform, normal, ylim = ylim, main = "boxplot",
  names = 1:3)
beanplot(bimodal, uniform, normal, ylim = ylim, main = "beanplot",
  col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6")


###################################################
### code chunk number 4: singers
###################################################
library("vioplot")
data("singer", package = "lattice")
ylim <- c(55, 80)
par(mfrow = c(2, 1), mai = c(0.8, 0.8, 0.5, 0.5))
data <- split(singer$height, singer$voice.part)
names(data)[1] <- "x"
do.call("vioplot", c(data,
  list(ylim = ylim, names = levels(singer$voice.part), col = "white")))
title(main = "vioplot", ylab = "body height (inch)")
beanplot(height ~ voice.part, data = singer, ll = 0.04, main = "beanplot",
  ylim = ylim, ylab = "body height (inch)")


###################################################
### code chunk number 5: distributionscode (eval = FALSE)
###################################################
## library("beanplot")
## par(mfrow = c(1, 2), mai = c(0.5, 0.5, 0.5, 0.1))
## mu <- 2
## si <- 0.6
## c <- 500
## bimodal <- c(rnorm(c/2, -mu, si), rnorm(c/2, mu, si))
## uniform <- runif(c, -4, 4)
## normal <- rnorm(c, 0, 1.5)
## ylim <- c(-7, 7)
## boxplot(bimodal, uniform, normal, ylim = ylim, main = "boxplot",
##   names = 1:3)
## beanplot(bimodal, uniform, normal, ylim = ylim, main = "beanplot",
##   col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6")


###################################################
### code chunk number 6: singersasym
###################################################
par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
beanplot(height ~ voice.part, data = singer, ll = 0.04,
  main = "beanplot", ylab = "body height (inch)", side = "both",
  border = NA, col = list("black", c("grey", "white")))
legend("bottomleft", fill = c("black", "grey"),
  legend = c("Group 2", "Group 1"))


###################################################
### code chunk number 7: singerscode (eval = FALSE)
###################################################
## library("vioplot")
## data("singer", package = "lattice")
## ylim <- c(55, 80)
## par(mfrow = c(2, 1), mai = c(0.8, 0.8, 0.5, 0.5))
## data <- split(singer$height, singer$voice.part)
## names(data)[1] <- "x"
## do.call("vioplot", c(data,
##   list(ylim = ylim, names = levels(singer$voice.part), col = "white")))
## title(main = "vioplot", ylab = "body height (inch)")
## beanplot(height ~ voice.part, data = singer, ll = 0.04, main = "beanplot",
##   ylim = ylim, ylab = "body height (inch)")


###################################################
### code chunk number 8: singersasymcode (eval = FALSE)
###################################################
## par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
## beanplot(height ~ voice.part, data = singer, ll = 0.04,
##   main = "beanplot", ylab = "body height (inch)", side = "both",
##   border = NA, col = list("black", c("grey", "white")))
## legend("bottomleft", fill = c("black", "grey"),
##   legend = c("Group 2", "Group 1"))


###################################################
### code chunk number 9: orchards
###################################################
beanplot(decrease ~ treatment, data = OrchardSprays, exp(rnorm(20, 3)),
  xlab = "threatment method", ylab = "decrease in potency",
  main = "OrchardSprays")


###################################################
### code chunk number 10: orchardscode (eval = FALSE)
###################################################
## beanplot(decrease ~ treatment, data = OrchardSprays, exp(rnorm(20, 3)),
##   xlab = "threatment method", ylab = "decrease in potency",
##   main = "OrchardSprays")


