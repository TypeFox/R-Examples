### R code from vignette source 'PMCMR.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: PMCMR.Rnw:103-105
###################################################
library("graphics")
boxplot(count ~ spray, data=InsectSprays)


###################################################
### code chunk number 2: PMCMR.Rnw:113-114
###################################################
kruskal.test(count ~ spray, data=InsectSprays)


###################################################
### code chunk number 3: PMCMR.Rnw:119-123
###################################################
require(PMCMR)
data(InsectSprays)
attach(InsectSprays)
posthoc.kruskal.nemenyi.test(x=count, g=spray, dist="Tukey")


###################################################
### code chunk number 4: PMCMR.Rnw:127-128
###################################################
posthoc.kruskal.nemenyi.test(count ~ spray, data=InsectSprays, dist="Tukey")


###################################################
### code chunk number 5: PMCMR.Rnw:133-134
###################################################
(out <- posthoc.kruskal.nemenyi.test(x=count, g=spray, dist="Chisquare"))


###################################################
### code chunk number 6: PMCMR.Rnw:139-140
###################################################
print(out$statistic)


###################################################
### code chunk number 7: PMCMR.Rnw:182-186
###################################################
require(PMCMR)
data(InsectSprays)
attach(InsectSprays)
posthoc.kruskal.dunn.test(x=count, g=spray, p.adjust.method="none")


###################################################
### code chunk number 8: PMCMR.Rnw:191-195
###################################################
require(PMCMR)
data(InsectSprays)
attach(InsectSprays)
posthoc.kruskal.dunn.test(x=count, g=spray, p.adjust.method="bonferroni")


###################################################
### code chunk number 9: PMCMR.Rnw:218-222
###################################################
require(PMCMR)
data(InsectSprays)
attach(InsectSprays)
posthoc.kruskal.conover.test(x=count, g=spray, p.adjust.method="none")


###################################################
### code chunk number 10: PMCMR.Rnw:225-229
###################################################
require(PMCMR)
data(InsectSprays)
attach(InsectSprays)
posthoc.kruskal.conover.test(x=count, g=spray, p.adjust.method="bonferroni")


###################################################
### code chunk number 11: PMCMR.Rnw:247-252
###################################################
require(stats) 
data(PlantGrowth)
attach(PlantGrowth)
kruskal.test(weight, group)
dunn.test.control(x=weight,g=group, p.adjust="bonferroni")


###################################################
### code chunk number 12: PMCMR.Rnw:258-259
###################################################
summary.lm(aov(weight ~ group))


###################################################
### code chunk number 13: PMCMR.Rnw:287-291
###################################################
require(PMCMR)
data(InsectSprays)
attach(InsectSprays)
vanWaerden.test(x=count, g=spray)


###################################################
### code chunk number 14: PMCMR.Rnw:305-309
###################################################
require(PMCMR)
data(InsectSprays)
attach(InsectSprays)
posthoc.vanWaerden.test(x=count, g=spray, p.adjust.method="none")


###################################################
### code chunk number 15: PMCMR.Rnw:362-371
###################################################
## Example from Sachs (1997, p. 402)
require(PMCMR)
x <- c(106, 114, 116, 127, 145, 110, 125,
       143, 148, 151, 136, 139, 149, 160,
       174)
g <- as.factor(c(rep(1,5), rep(2,5), rep(3,5)))
levels(g) <- c("A", "B", "C")
jonckheere.test(x , g, "increasing")
rm(x,g)


###################################################
### code chunk number 16: PMCMR.Rnw:406-415
###################################################
require(PMCMR)
y <- matrix(c(
3.88, 5.64, 5.76, 4.25, 5.91, 4.33, 30.58, 30.14, 16.92,
23.19, 26.74, 10.91, 25.24, 33.52, 25.45, 18.85, 20.45, 
26.67, 4.44, 7.94, 4.04, 4.4, 4.23, 4.36, 29.41, 30.72,
32.92, 28.23, 23.35, 12, 38.87, 33.12, 39.15, 28.06, 38.23,
26.65),nrow=6, ncol=6, 
dimnames=list(1:6,c("A","B","C","D","E","F")))
print(y)


###################################################
### code chunk number 17: PMCMR.Rnw:422-425
###################################################
library("graphics")
groups <- gl(6,6,labels=colnames(y))
boxplot(as.vector(y) ~ groups)


###################################################
### code chunk number 18: PMCMR.Rnw:432-433
###################################################
friedman.test(y)


###################################################
### code chunk number 19: PMCMR.Rnw:438-439
###################################################
posthoc.friedman.nemenyi.test(y)


###################################################
### code chunk number 20: PMCMR.Rnw:456-466
###################################################
require(PMCMR)
y <- matrix(c(
3.88, 5.64, 5.76, 4.25, 5.91, 4.33, 30.58, 30.14, 16.92,
23.19, 26.74, 10.91, 25.24, 33.52, 25.45, 18.85, 20.45, 
26.67, 4.44, 7.94, 4.04, 4.4, 4.23, 4.36, 29.41, 30.72,
32.92, 28.23, 23.35, 12, 38.87, 33.12, 39.15, 28.06, 38.23,
26.65),nrow=6, ncol=6, 
dimnames=list(1:6,c("A","B","C","D","E","F")))
friedman.test(y)
posthoc.friedman.conover.test(y=y, p.adjust="none")


###################################################
### code chunk number 21: PMCMR.Rnw:520-537
###################################################
## Conover (1999, p. 375f):
## Numbers of five brands of a new hand lotion sold in seven stores
## during one week.
y <- matrix(c( 5,  4,  7, 10, 12,
               1,  3,  1,  0,  2,
              16, 12, 22, 22, 35,
               5,  4,  3,  5,  4,
              10,  9,  7, 13, 10,
              19, 18, 28, 37, 58,
              10,  7,  6,  8,  7),
            nrow = 7, byrow = TRUE,
            dimnames =
            list(Store = as.character(1:7),
                 Brand = LETTERS[1:5]))
y
quade.test(y)
posthoc.quade.test(y, dist="TDist", p.adj="none")


###################################################
### code chunk number 22: PMCMR.Rnw:576-587
###################################################
## Example for an incomplete block design:
## Data from Conover (1999, p. 391).
y <- matrix(c(
2,NA,NA,NA,3, NA,  3,  3,  3, NA, NA, NA,  3, NA, NA,
  1,  2, NA, NA, NA,  1,  1, NA,  1,  1,
NA, NA, NA, NA,  2, NA,  2,  1, NA, NA, NA, NA,
  3, NA,  2,  1, NA, NA, NA, NA,  3, NA,  2,  2
), ncol=7, nrow=7, byrow=FALSE,
dimnames=list(1:7, LETTERS[1:7]))
y
durbin.test(y)


###################################################
### code chunk number 23: PMCMR.Rnw:601-602
###################################################
posthoc.durbin.test(y, p.adj="none")


###################################################
### code chunk number 24: PMCMR.Rnw:609-610
###################################################
print(posthoc.durbin.test(y, p.adj="none"))


###################################################
### code chunk number 25: PMCMR.Rnw:612-613
###################################################
summary(posthoc.durbin.test(y, p.adj="none"))


###################################################
### code chunk number 26: PMCMR.Rnw:619-632
###################################################
require(multcompView)
data(InsectSprays)
attach(InsectSprays)
out <- posthoc.kruskal.dunn.test(count ~ spray, p.adjust="bonf")
out.p <- get.pvalues(out)
out.mcV <- multcompLetters(out.p, threshold=0.05)
Rij <- rank(count)
Rj.mean <- tapply(Rij, spray, mean)
xx <- barplot(Rj.mean, ylim=c(0, 1.2* max(Rj.mean)), 
              xlab="Spray", ylab="Mean rank")
yy <- Rj.mean + 3
text(xx, yy, lab=out.mcV$Letters)
detach(InsectSprays)


###################################################
### code chunk number 27: Dunn.tab
###################################################
require(xtable)
require(multcompView)
data(InsectSprays)
attach(InsectSprays)
out <- posthoc.kruskal.dunn.test(count ~ spray, p.adjust="bonf")
out.p <- get.pvalues(out)
out.mcV <- multcompLetters(out.p, threshold=0.05)
Rij <- rank(count)
Rj.mean <- tapply(Rij, spray, mean)
dat <- data.frame(Group = names(Rj.mean),
                  meanRj = Rj.mean,
                  M = out.mcV$Letters)
dat.x <- xtable(dat)
caption(dat.x) <- "Mean ranks ($\\bar{R}_{j}$) of the 
\\texttt{InsectSprays} data set. Different letters (M) 
indicate significant differences ($p < 0.05$) according 
to the Bonferroni-Dunn test (see Chap. \\ref{Dunn})."
colnames(dat.x) <- c("Group", "$\\bar{R}_{j}$", "M")
digits(dat.x) <- 1
label(dat.x) <- "tab:bonf.dunn"
print(dat.x, include.rownames=F, caption.placement="top", 
     	      sanitize.text.function = function(x){x}, 
              table.placement="h")
detach(InsectSprays)


