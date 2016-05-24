####################################################################################
## PLOT OF HOTSPOTS
## vignette qtlhot: no parallel need
scan1 <- scanone(cross1, pheno.col = 1:1000, method = "hk")
high1 <- pull.highlods(scan1,, lod.thr = 3.47, drop.lod = 1.5)
hots1 <- hotsize(high1, lod.thr = 3.47, window = 5, quant.level = NL.N.1)
plot(hots1)

## Permutations
NL.N.1 <- NL.N.perm(cross = cross1,
                      Nmax = 300,
                      n.perm = 100,
                      lod.thrs = lod.thrs,
                      drop = 1.5,
                      verbose = FALSE)
hot.scan1 <- qtlhot.scan(cross1, scan1, NL.N.1$max.lod.quant, lod.thrs, probs = seq(.01, .1, by = .01), level = 0.05)
par(mar = c(4.1,4.1,0.1,4.1))
plot(hot.scan1, lodcolumn = 1, add = TRUE, col = "red")
## Add right axis with sliding LOD thresholds.
quant <- attr(hot.scan1, "quant")
tmp <- seq(along = quant)
axis(4, at = tmp, label = round(quant[tmp], 1), las = 1, cex = 0.35)
mtext("sliding LOD thresholds",4, 1, cex=1.5)

## Sliding bar plot
NL.N.1.thrs <- NL.N.summary(NL.N.1[[1]], NL.N.1[[2]], alphas)
NL.1.thr <- NL.N.1.thrs[[1]]
N.1.thr <- NL.N.1.thrs[[2]]
N.1 <- round(N.1.thr[5, 5])
par(mar=c(4.1,4.1,0.1,4.1))
sliding.bar.plot(scan = data.frame(scan1[, 1:2], scandrop1),
                 lod.thr = NL.1.thr[1:N.1, 5],
                 size.thr = 1:N.1,
                 gap = 50,
                 y.axes = seq(1, N.1, by = 10))

#####################################################################################333

## parallel code
## Really want to have plot an summary for object scan.hl, which comes from Phase2 unpermuted data.
per.scan <- scanone(perm.cross, pheno.col=pheno.col, method="hk",
                    addcovar=covars$addcovar, intcovar=covars$intcovar)

scan.hl <- pull.highlods(per.scan, lod = lod.min, drop.lod = drop.lod,
                         restrict.lod = TRUE)

## lod.thrs comes from single trait permutations
## should have probs along with it!
probs <- seq(.8,.99, by = .01)
lod.thrs <- read.csv(file.path(null.dir, "quantiles.csv"), header = TRUE)

## max.lod.quant, max.N, max.N.window come from multiple trait permutations
load("Phase31.RData")

## create quant.level and quant.thr and pull out lod.thr for level
out <- quant.plot(max.lod.quant, max.N, max.N.window, lod.thrs, probs, level = 0.95)

## create object of class scanone with max.N, quant, max.N.window
hotspots <- hotspot.scan(cross, scan.hl, out$lod.thr, out$quant.level, window = 5)
## wrapper for plot.scanone for hotspots
hotspot.plot(hotspots, quant.thr = out$quant.thr, maps = B6BTBR07n.maps, main = main)

