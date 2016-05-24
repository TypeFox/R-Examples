### R code from vignette source 'qtlhot.Rnw'

###################################################
### code chunk number 1: loadlib
###################################################
library(qtlhot)


###################################################
### code chunk number 2: qtlhot.Rnw:61-69
###################################################
ncross1 <- sim.null.cross(chr.len = rep(100, 4),
                          n.mar = 51,
                          n.ind = 100,
                          type = "bc",
                          n.phe = 1000,
                          latent.eff = 3,
                          res.var = 1,
                          init.seed = 123457)


###################################################
### code chunk number 3: qtlhot.Rnw:74-86
###################################################
cross1 <- include.hotspots(cross = ncross1,
                           hchr = c(2, 3, 4),
                           hpos = c(25, 75, 50),
                           hsize = c(100, 50, 20),
                           Q.eff = 2,
                           latent.eff = 3,
                           lod.range.1 = c(2.5, 2.5),
                           lod.range.2 = c(5, 8),
                           lod.range.3 = c(10, 15),
                           res.var = 1,
                           n.phe = 1000,
                           init.seed = 12345)


###################################################
### code chunk number 4: qtlhot.Rnw:91-94
###################################################
ncor1 <- cor(cross1$pheno)
summary(ncor1[lower.tri(ncor1)])
rm(ncor1)


###################################################
### code chunk number 5: qtlhot.Rnw:96-97
###################################################
if(file.exists("savedperms.RData")) load("savedperms.RData")


###################################################
### code chunk number 6: qtlhot.Rnw:102-106
###################################################
if(!exists("pt")) {
  set.seed(123)
  pt <- scanone(ncross1, method = "hk", n.perm = 1000)
}


###################################################
### code chunk number 7: qtlhot.Rnw:115-119
###################################################
alphas <- seq(0.01, 0.10, by=0.01)
lod.thrs <- summary(pt, alphas)
lod.thrs
lod.thr <- lod.thrs[5]


###################################################
### code chunk number 8: qtlhot.Rnw:124-125
###################################################
scan1 <- scanone(cross1, pheno.col = 1:1000, method = "hk")


###################################################
### code chunk number 9: qtlhot.Rnw:130-132
###################################################
high1 <- highlod(scan1, lod.thr = min(lod.thrs), drop.lod = 1.5)
max(high1, lod.thr = lod.thrs)


###################################################
### code chunk number 10: qtlhot.Rnw:137-139
###################################################
hots1 <- hotsize(high1, lod.thr = lod.thr)
summary(hots1)


###################################################
### code chunk number 11: plotex1counts
###################################################
par(mar=c(4.1,4.1,0.5,0.1))
plot(hots1, cex.lab = 1.5, cex.axis = 1.5)


###################################################
### code chunk number 12: qtlhot.Rnw:163-173
###################################################
if(!exists("hotperm1")) {
  set.seed(12345)
  hotperm1 <- hotperm(cross = cross1,
                  n.quant = 300,
                  n.perm = 100,
                  lod.thrs = lod.thrs,
                  alpha.levels = alphas,
                  drop.lod = 1.5,
                  verbose = FALSE)
}


###################################################
### code chunk number 13: qtlhot.Rnw:187-189
###################################################
names(hotperm1)
summary(hotperm1)


###################################################
### code chunk number 14: figex1slidingbar
###################################################
par(mar=c(4.1,4.1,0.5,4.1))
quant1 <- quantile(hotperm1, 0.05, lod.thr = lod.thr)
plot(high1, quant.level = quant1, sliding = TRUE)


###################################################
### code chunk number 15: plotex1slidingbar
###################################################
par(mar=c(4.1,4.1,0.5,4.1))
quant1 <- quantile(hotperm1, 0.05, lod.thr = lod.thr)
plot(high1, quant.level = quant1, sliding = TRUE)


###################################################
### code chunk number 16: qtlhot.Rnw:222-225
###################################################
hotsq1 <- hotsize(high1, lod = lod.thr, window = 5, quant.level = quant1)
quant.axis <- pmax(1, pretty(c(0, qtl:::summary.scanone(hotsq1)[3,"quant"])))
quant.level <- round(attr(hotsq1, "quant.level")[quant.axis], 2)


###################################################
### code chunk number 17: figex1elias
###################################################
par(mar=c(4.1,4.1,0.5,4.1))
plot(hotsq1)


###################################################
### code chunk number 18: plotex1elias
###################################################
par(mar=c(4.1,4.1,0.5,4.1))
plot(hotsq1)


###################################################
### code chunk number 19: qtlhot.Rnw:253-254
###################################################
summary(hotsq1)


###################################################
### code chunk number 20: qtlhot.Rnw:261-281
###################################################
ncross2 <- sim.null.cross(chr.len = rep(100,4), 
                          n.mar = 51, 
                          n.ind = 100,
                          type = "bc", 
                          n.phe = 1000, 
                          latent.eff = 0, 
                          res.var = 1, 
                          init.seed = 123457)
cross2 <- include.hotspots(cross = ncross2,
                           hchr = c(2, 3, 4),
                           hpos = c(25, 75, 50),
                           hsize = c(100, 50, 20),
                           Q.eff = 2,
                           latent.eff = 0,
                           lod.range.1 = c(2.5, 2.5),
                           lod.range.2 = c(5, 8),
                           lod.range.3 = c(10, 15),
                           res.var = 1,
                           n.phe = 1000,
                           init.seed = 12345)


###################################################
### code chunk number 21: figex2counts
###################################################
scan2 <- scanone(cross2, pheno.col = 1:1000, method = "hk")
high2 <- highlod(scan2, lod.thr = lod.thr, drop.lod = 1.5)
hots2 <- hotsize(high2)
par(mar=c(4.1,4.1,0.1,0.1))
plot(hots2, cex.lab = 1.5, cex.axis = 1.5)


###################################################
### code chunk number 22: plotex2counts
###################################################
scan2 <- scanone(cross2, pheno.col = 1:1000, method = "hk")
high2 <- highlod(scan2, lod.thr = lod.thr, drop.lod = 1.5)
hots2 <- hotsize(high2)
par(mar=c(4.1,4.1,0.1,0.1))
plot(hots2, cex.lab = 1.5, cex.axis = 1.5)


###################################################
### code chunk number 23: qtlhot.Rnw:326-329
###################################################
ncor2 <- cor(cross2$pheno)
summary(ncor2[lower.tri(ncor2)])
rm(ncor2)


###################################################
### code chunk number 24: qtlhot.Rnw:340-350
###################################################
if(!exists("hotperm2")) {
  set.seed(12345)
  hotperm2 <- hotperm(cross = cross2, 
                  n.quant = 300, 
                  n.perm = 100, 
                  lod.thrs = lod.thrs, 
                  alpha.levels = alphas,
                  drop.lod = 1.5, 
                  verbose = FALSE) 
}


###################################################
### code chunk number 25: qtlhot.Rnw:352-353
###################################################
quant2 <- quantile(hotperm2, 0.05, lod.thr = lod.thr)


###################################################
### code chunk number 26: plotex2slidingbar
###################################################
par(mar=c(4.1,4.1,0.5,0.1))
plot(high2, lod.thr = lod.thr, quant.level = quant2, sliding = TRUE)


###################################################
### code chunk number 27: plotex2sigct
###################################################
par(mar=c(4.1,4.1,0.5,0.1))
plot(high2, quant.level = quant2)


###################################################
### code chunk number 28: qtlhot.Rnw:395-397
###################################################
if(!file.exists("savedperms.RData"))
  save(pt, hotperm1, hotperm2, file = "savedperms.RData", compress = TRUE)


