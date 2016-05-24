# $Id: interval_summary.R 169 2006-01-05 05:02:44Z andrewr $

rm(list=ls())

load("../scripts/booted_srs.RData")
require(lattice)

opar <- par(mfrow=c(3,3), mar=c(5,4,0,0))
for (i in 1:length(levels(plots$forest))) {
  plot(density(plots$ba[plots$forest==unique(plots$forest)[i]]),
       main="", ylab="",
       xlab=paste("National forest ID:", unique(plots$forest)[i]))
}
par(opar)

opar <- par(mfrow=c(3,3), mar=c(5,4,0,0))
for (i in 1:length(levels(plots$forest))) {
  qqnorm(plots$ba[plots$forest==unique(plots$forest)[i]],
         main="", ylab="",
         xlab=paste("National forest ID:", unique(plots$forest)[i]))
  qqline(plots$ba[plots$forest==unique(plots$forest)[i]])
}
par(opar)

xyplot(contained ~ n | interval,
       panel = function(x, y) {
         panel.xyplot(x, y, type="l")
         panel.abline(h=0.95, col="darkgrey")
       },
       index.cond=list(c(1:5,7:9,6,10:12)),
       data=overall.results)

xyplot(length ~ n | interval,
       panel = function(x, y) {
         panel.xyplot(x, y, type="l")
       },
       index.cond=list(c(1:5,7:9,6,10:12)),
       data=overall.results)

xyplot(contained ~ length | interval,
       panel = function(x, y) {
         panel.xyplot(x, y, type="l")
         panel.abline(h=0.95, col="darkgrey")
       },
       index.cond=list(c(1:5,7:9,6,10:12)),
       data=overall.results)

xyplot(length ~ n | interval,
       groups=forest,
       auto.key = list(points = TRUE, columns=5,
       title="Nearest National Forest", space="top"),
       xlab="Sample size, n",
       ylab="Mean interval length ()",
       panel = function(x, y, subscripts, groups) {
         for (i in 1:length(levels(results$forest))) {
           panel.xyplot(x[groups[subscripts]==levels(results$forest)[i]],
                        y[groups[subscripts]==levels(results$forest)[i]],
                        col=trellis.par.get()$superpose.line$col[i],
                        type="l")
         }
       },
       index.cond=list(c(1:5,7:9,6,10:12)),
       data=results)

xyplot(contained ~ n | interval,
       groups=forest,
       auto.key = list(points = TRUE, columns=5,
       title="Nearest National Forest", space="top"),
       panel = function(x, y, subscripts, groups) {
         for (i in 1:length(levels(results$forest))) {
           panel.xyplot(x[groups[subscripts]==levels(results$forest)[i]],
                        y[groups[subscripts]==levels(results$forest)[i]],
                        col=trellis.par.get()$superpose.line$col[i],
                        type="l")
         }
         panel.abline(h=0.95, col="darkgrey")
       },
       index.cond=list(c(1:5,7:9,6,10:12)),
       data=results)

xyplot(contained ~ length | interval,
       groups=forest,
       auto.key = list(points = TRUE, columns=5,
       title="Nearest National Forest", space="top"),
       panel = function(x, y, subscripts, groups) {
         for (i in 1:length(levels(results$forest))) {
           panel.xyplot(x[groups[subscripts]==levels(results$forest)[i]],
                        y[groups[subscripts]==levels(results$forest)[i]],
                        col=trellis.par.get()$superpose.line$col[i],
                        type="l")
         }
         panel.abline(h=0.95, col="darkgrey")
       },
       index.cond=list(c(1:5,7:9,6,10:12)),
       data=results)
