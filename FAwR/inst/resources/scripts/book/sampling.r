### R code from vignette source 'sampling.rnw'

###################################################
### code chunk number 1: Setup
###################################################
options(repos="http://cran.r-project.org")
if(!require(Hmisc, quietly=TRUE)) install.packages("Hmisc")
if(!require(survey, quietly=TRUE)) install.packages("survey")
if(!require(gmodels, quietly=TRUE)) install.packages("gmodels")
if(!require(lattice, quietly=TRUE)) install.packages("lattice")
source("../../scripts/functions.R")
options(width=60)                    
lattice.options(default.theme = canonical.theme(color = FALSE))


###################################################
### code chunk number 2: ufc
###################################################
Stangle("../../ch2/sweave/ufc.rnw")
source("ufc.R")


###################################################
### code chunk number 3: pref
###################################################
Stangle("../../ch2/sweave/pref.rnw")
source("pref.R")


###################################################
### code chunk number 4: fia
###################################################
Stangle("../../ch2/sweave/fia.rnw")
source("fia.R")


###################################################
### code chunk number 5: SRS.A.hat
###################################################
boot.mean <- function(x, index) {
  mean(x[index])
}


###################################################
### code chunk number 6: SRS.A.hat
###################################################
library(boot)
pf.pt.1 <- subset(pref.point, point==1)
rownames(pf.pt.1) <- 1:nrow(pf.pt.1)
pref.SRS.boot <- boot(pf.pt.1$vol.m3.ha, 
                      boot.mean, 
                      R = 1999)


###################################################
### code chunk number 7: SRS.A.hat
###################################################
pref.SRS.boot


###################################################
### code chunk number 8: bootgraph
###################################################
par(las=1)
plot(pref.SRS.boot, pch=20)


###################################################
### code chunk number 9: SRS.A.hat
###################################################
sd(pref.SRS.boot$t)


###################################################
### code chunk number 10: SRS.A.hat
###################################################
boot.ci(pref.SRS.boot, type="norm")


###################################################
### code chunk number 11: jackbootgraph
###################################################
jack.after.boot(pref.SRS.boot)


###################################################
### code chunk number 12: qq-graph
###################################################
normalized <- qqnorm(pf.pt.1$vol.m3.ha, plot.it = FALSE)
plot(normalized$x, normalized$y, type="n",
     ylab="Sample Quantiles", xlab="Theoretical Quantiles", 
     main="Normal Q-Q Plot")
qqline(pref.point$vol.m3.ha, col="darkgray")
text(x = normalized$x, y = normalized$y, cex = 0.85)


###################################################
### code chunk number 13: sampling.rnw:748-749
###################################################
pf.pt.1[order(pf.pt.1$vol.m3.ha, decreasing = TRUE),][1:5,]


###################################################
### code chunk number 14: qq-graph
###################################################
normalized <- qqnorm(pf.pt.1$vol.m3.ha, plot.it = FALSE)
plot(normalized$x, normalized$y, type="n",
     ylab="Sample Quantiles", xlab="Theoretical Quantiles", 
     main="Normal Q-Q Plot")
qqline(pref.point$vol.m3.ha, col="darkgray")
text(x = normalized$x, y = normalized$y, cex = 0.85)


###################################################
### code chunk number 15: SRS.A.hat
###################################################
bm.2 <- function(x, index)
  c(mean(x[index]), var(x[index])/length(x))


###################################################
### code chunk number 16: SRS.A.hat
###################################################
pref.SRS.b2 <- boot(pf.pt.1$vol.m3.ha, bm.2, 
                        R = 1999)
pref.SRS.b2


###################################################
### code chunk number 17: fig-varstab
###################################################
scatter.smooth(pref.SRS.b2$t[,1], 
               sqrt(pref.SRS.b2$t[,2]), 
               ylab = "Standard Error", 
               xlab = "Bootstrap Mean")


###################################################
### code chunk number 18: fig-sqrtvarstab
###################################################
pref.SRS.b3 <- boot(I(sqrt(pf.pt.1$vol.m3.ha)), 
                        bm.2, 
                        R = 1999)
scatter.smooth(pref.SRS.b3$t[,1], 
               sqrt(pref.SRS.b3$t[,2]), 
               ylab = "Standard Error", 
               xlab = "Bootstrap Mean")


###################################################
### code chunk number 19: varstab
###################################################
scatter.smooth(pref.SRS.b2$t[,1], 
               sqrt(pref.SRS.b2$t[,2]), 
               ylab = "Standard Error", 
               xlab = "Bootstrap Mean")


###################################################
### code chunk number 20: logvarstab
###################################################
pref.SRS.b3 <- boot(I(sqrt(pf.pt.1$vol.m3.ha)), 
                        bm.2, 
                        R = 1999)
scatter.smooth(pref.SRS.b3$t[,1], 
               sqrt(pref.SRS.b3$t[,2]), 
               ylab = "Standard Error", 
               xlab = "Bootstrap Mean")


###################################################
### code chunk number 21: sampling.rnw:939-949
###################################################
ratio.test <- read.csv("../../data/ratio-test.csv")
ratio.test.by.forest <- 
  read.csv("../../data/ratio-test-by-forest.csv")
ratio.test.by.forest$forest <- 
  factor(ratio.test.by.forest$forest)
forests <- levels(fia.plots$forest)
levels(ratio.test$interval) <- 
  c("Basic","BCa","Classical","Jackknife","Linearized","Basic (Log)",
    "Normal (Log)","Studentized (Log)","Normal","Percentile",
    "Basic (Sqrt)","Normal (Sqrt)","Studentized (Sqrt)","Studentized")


###################################################
### code chunk number 22: xyplot-coverage
###################################################
xyplot(contained ~ n | interval,    
       xlab = "Sample size (n)",
       ylab = "Realized Coverage Rate",
       panel = function(x, y) {
         panel.xyplot(x, y, type="l")
         panel.abline(h = 0.95, col = "darkgrey")
       },
       index.cond = list(c(3:5,10,1,9,14,2,6:8)),
       skip = c(rep(FALSE, 11), TRUE, rep(FALSE, 3)),
       data = ratio.test)


###################################################
### code chunk number 23: xyplot-coverage
###################################################
print(
xyplot(contained ~ n | interval,    
       xlab = "Sample size (n)",
       ylab = "Realized Coverage Rate",
       panel = function(x, y) {
         panel.xyplot(x, y, type="l")
         panel.abline(h = 0.95, col = "darkgrey")
       },
       index.cond = list(c(3:5,10,1,9,14,2,6:8)),
       skip = c(rep(FALSE, 11), TRUE, rep(FALSE, 3)),
       data = ratio.test)
      )


###################################################
### code chunk number 24: SRS.A
###################################################
library(survey)
pref.SRS <- svydesign(id = ~1, 
                      data = pref.point, 
                      weight = pref.point$weight)


###################################################
### code chunk number 25: SRS.A
###################################################
svymean(~vol.m3.ha, pref.SRS)


###################################################
### code chunk number 26: SRS.A.hat
###################################################
(SRS.v.hat <- mean(pref.point$vol.m3.ha))
SRS.n <- length(pref.point$vol.m3.ha)
(SRS.v.se <- sd(pref.point$vol.m3.ha) / sqrt(SRS.n))
SRS.v.hat + SRS.v.se * qt(0.975, df = SRS.n - 1) * c(-1, 1)


###################################################
### code chunk number 27: SyRS.A.hat
###################################################
(SyRS.v.hat <- mean(ufc.SyRS.data$vol.m3.ha))


###################################################
### code chunk number 28: SyRS.A.SE
###################################################
SyRS.v.n <- length(ufc.SyRS.data$vol.m3.ha)
vol.first.diffs <- ufc.SyRS.data$vol.m3.ha[1:(SyRS.v.n-1)] - 
  ufc.SyRS.data$vol.m3.ha[2:SyRS.v.n]
SyRS.v.se <- sqrt(1 / (2 * SyRS.v.n * (SyRS.v.n - 1)) * 
                  sum(vol.first.diffs^2))
SyRS.v.se


###################################################
### code chunk number 29: SyRS.A.SE
###################################################
sd(ufc.SyRS.data$vol.m3.ha) / sqrt(SyRS.v.n)


###################################################
### code chunk number 30: ClS.A.hat
###################################################
pref.CS <- svydesign(id = ~cluster, 
                     data = pref.point, 
                     weight = pref.point$weight)


###################################################
### code chunk number 31: sampling.rnw:1258-1259
###################################################
svymean( ~ vol.m3.ha, pref.CS)


###################################################
### code chunk number 32: ClS.A.SE
###################################################
v.clust <- 
  aggregate(x = list(vol.m3.ha = pref.point$vol.m3.ha),
            by = list(cluster = pref.point$cluster),
            FUN = sum)

v.clust$count <- table(pref.point$cluster)
n.clusters <- length(v.clust$count)

(v.bar.c <- sum(v.clust$vol.m3.ha) / sum(v.clust$count))

(se.v.bar.c <- 
  sqrt(sum(v.clust$count^2 / mean(v.clust$count)^2 *
           (v.clust$vol.m3.ha / v.clust$count - 
            v.bar.c)^2) / n.clusters / (n.clusters - 1)))


###################################################
### code chunk number 33: Clus.CI
###################################################
v.bar.c + c(-1, 1) * 2 * se.v.bar.c


###################################################
### code chunk number 34: 2SS.A.hat
###################################################
pref.2SS <- svydesign(id = ~cluster + point, 
                      fpc = ~rep(36/3124, nrow(pref.point)) +
                             rep(0.000001, nrow(pref.point)),
                      weight = pref.point$weight,                      
                      data = pref.point)


###################################################
### code chunk number 35: sampling.rnw:1418-1419
###################################################
svymean(~vol.m3.ha, pref.2SS, na.rm = TRUE)


###################################################
### code chunk number 36: sampling.rnw:1453-1478
###################################################
v.2SS <- 
  aggregate(x = list(vol.m3.ha.bar = pref.point$vol.m3.ha),
            by = list(cluster = pref.point$cluster),
            FUN = mean)
v.2SS$count <- table(pref.point$cluster)

n.2SS <- nrow(v.2SS)
N.2SS <- 3124

m.2SS <- mean(v.2SS$count)
M.2SS <- 10000

f1.2SS <- n.2SS / N.2SS
f2.2SS <- m.2SS / M.2SS

pref.point.temp <- merge(pref.point, v.2SS)

v1 <- var(v.2SS$vol.m3.ha.bar) 
v2 <- 1 / n.2SS / (m.2SS - 1) *
  sum((pref.point.temp$vol.m3.ha - 
       pref.point.temp$vol.m3.ha.bar)^2)

se.v.bar.2SS <- sqrt((1 - f1.2SS)/n.2SS*v1 + 
                     f1.2SS*(1-f2.2SS)/(n.2SS*m.2SS)*v2)
se.v.bar.2SS


###################################################
### code chunk number 37: sampling.rnw:1487-1488
###################################################
sd(v.2SS$vol.m3.ha) / sqrt(n.2SS)


###################################################
### code chunk number 38: sampling.rnw:1529-1533
###################################################
summary(aov.tab <- aov(vol.m3.ha ~ cluster, data=pref.point))
unlist(aov.tab)[5:6]
v1*5
v2


###################################################
### code chunk number 39: StS.A.hat
###################################################
pref.StRS <- svydesign(id = ~1, 
                       strata = ~stratum, 
                       data = pref.point, 
                       weight = pref.point$weight)


###################################################
### code chunk number 40: StS.A.hat
###################################################
svymean(~vol.m3.ha, pref.StRS)


###################################################
### code chunk number 41: Str.A.SE
###################################################
stratum.weights <- rep(1/9, 9)

v.str <- sum(stratum.weights * 
             tapply(pref.point$vol.m3.ha, 
                    pref.point$stratum, 
                    mean))

se.v.str <- sqrt(sum(stratum.weights^2 *
                      tapply(pref.point$vol.m3.ha, 
                             pref.point$stratum, var) / 
                     tapply(pref.point$vol.m3.ha, 
                             pref.point$stratum, length)))


###################################################
### code chunk number 42: StS.A.hat
###################################################
pref.CStRS <- svydesign(id = ~cluster, 
                        strata = ~stratum, 
                        data = pref.point, 
                        weight = pref.point$weight)
svymean(~vol.m3.ha, pref.CStRS)


###################################################
### code chunk number 43: RatE.A.hat
###################################################
pref.SRSc <- svydesign(id = ~1, 
                       data = pref.point.cov, 
                       weight = pref.point$weight)


###################################################
### code chunk number 44: RatE.A.hat
###################################################
volume.over.acndviC <- svyratio(numerator = ~vol.m3.ha, 
                                denominator = ~acndviC, 
                                design = pref.SRSc)
predict(volume.over.acndviC, mean(pref.pixel$ndvic))


###################################################
### code chunk number 45: RatE.A.hat
###################################################
ratio.hat <- sum(pref.point.cov$vol.m3.ha) /  
  sum(pref.point.cov$acndviC)
RatE.v.hat <- mean(pref.pixel$ndvic) * ratio.hat
RatE.v.hat
RatE.v.hat.se <- 
  sqrt((mean(pref.pixel$ndvic))^2 / 
       (mean(pref.point.cov$acndviC))^2 *
       (var(pref.point.cov$vol.m3.ha) + 
        ratio.hat^2 * var(pref.point.cov$acndviC) -
        2 * ratio.hat * cov(pref.point.cov$acndviC, 
                            pref.point.cov$vol.m3.ha)) / 
       nrow(pref.point.cov))
RatE.v.hat.se


###################################################
### code chunk number 46: StS.A.hat
###################################################
pref.CStRSc <- svydesign(id = ~cluster+point, 
                         strata = ~stratum, 
                         data = pref.point.cov, 
                         weight = pref.point$weight)


###################################################
### code chunk number 47: StS.A.hat
###################################################
v.over.a.CStRSc <- svyratio(numerator = ~vol.m3.ha, 
                            denominator = ~acndviC, 
                            design=pref.CStRSc)
predict(v.over.a.CStRSc, mean(pref.pixel$ndvic))


###################################################
### code chunk number 48: RegE.A.hat
###################################################
volume.against.acndviC <- svyglm(vol.m3.ha ~ acndviC, 
                                 design = pref.SRSc)


###################################################
### code chunk number 49: RegE.A.hat
###################################################
library(gmodels)


###################################################
### code chunk number 50: sampling.rnw:1891-1894
###################################################
estimable(volume.against.acndviC,
          conf.int = 0.95,
          cm = rbind(total=c(1, mean(pref.pixel$ndvic))))


###################################################
### code chunk number 51: StS.A.hat
###################################################
pref.CStRSc <- svydesign(id = ~cluster+point, 
                         strata = ~stratum, 
                         weight = pref.point$weight, 
                         data = pref.point.cov)
v.against.a.CStRSc <- svyglm(vol.m3.ha ~ acndviC, 
                             design = pref.CStRSc)
estimable(v.against.a.CStRSc, 
          conf.int = 0.95,
          cm = rbind(total = c(1, mean(pref.pixel$ndvic))))


###################################################
### code chunk number 52: sampling.rnw:1946-1948
###################################################
q.max <- 1000/20
PPP.rand <- q.max * runif(1200)


###################################################
### code chunk number 53: sampling.rnw:1957-1960
###################################################
PPP.guesses <- pref.tree$vol.m3 * 
  rnorm(n = nrow(pref.tree), mean = 1.05, sd = 0.15) + 
  rnorm(n = nrow(pref.tree), mean = -0.05, sd = 0.05)


###################################################
### code chunk number 54: sampling.rnw:1965-1970
###################################################
PPP.sample <- 
  PPP.rand[seq(1,length(PPP.guesses))] < PPP.guesses
PPP.trees <- 
  data.frame(vol.m3 = pref.tree$vol.m3[PPP.sample],
             guess.m3 = PPP.guesses[PPP.sample])


###################################################
### code chunk number 55: sampling.rnw:1985-1987
###################################################
(PPP.hat.1 <- sum(PPP.guesses) / sum(PPP.sample) *
  sum(PPP.trees$vol.m3 / PPP.trees$guess.m3))


###################################################
### code chunk number 56: 3P.A.hat
###################################################
(PPP.hat.2 <- PPP.hat.1 * 
 sum(PPP.guesses) / q.max / sum(PPP.sample))


###################################################
### code chunk number 57: sampling.rnw:2006-2007
###################################################
sum(pref.tree$vol.m3)


###################################################
### code chunk number 58: 3P.A.SE
###################################################
se.PPP.hat <- sd(PPP.trees$vol.m3 / PPP.trees$guess.m3) / 
  sqrt(sum(PPP.sample)) * sum(PPP.guesses)
se.PPP.hat


###################################################
### code chunk number 59: VBAR.A.hat
###################################################
vbar.sample <- sample(1:dim(pref.point)[1], size=30)
vbar.points <- pref.point[vbar.sample,]
vbar.points <- vbar.points[vbar.points$ba.m2.ha > 0,]
vbar.points$ratio <- vbar.points$vol.m3.ha / 
  vbar.points$ba.m2.ha


###################################################
### code chunk number 60: VBAR.A.hat
###################################################
(v.hat.vbar <- 
  mean(pref.point$ba.m2.ha) * mean(vbar.points$ratio))


###################################################
### code chunk number 61: VBAR.A.SE
###################################################
s.bar.g.perc.2 <- (sd(pref.point$ba.m2.ha) / 
                   sqrt(length(pref.point$ba.m2.ha)) / 
                   mean(pref.point$ba.m2.ha))^2
s.bar.r.perc.2 <- (sd(vbar.points$ratio) / 
                   sqrt(length(vbar.points$ratio)) / 
                   mean(vbar.points$ratio))^2
s.bar.v.perc <- sqrt(s.bar.g.perc.2 + s.bar.r.perc.2 -
                     s.bar.g.perc.2 * s.bar.r.perc.2)
s.bar.v <- s.bar.v.perc * v.hat.vbar
s.bar.v


###################################################
### code chunk number 62: sampling.rnw:2143-2145
###################################################
system("rm -fr package-Ch3")
package.skeleton(name = "package-Ch3")


