# $Id: ratio_summary.R 169 2006-01-05 05:02:44Z andrewr $

rm(list=ls())

load("../scripts/booted_ratio.RData")
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

opar <- par(mfrow=c(3,3), mar=c(5,4,0,0))
for (i in 1:length(levels(plots$forest))) {
  ratio.res <- 
  qqnorm(plots$ba[plots$forest==unique(plots$forest)[i]],
         main="", ylab="",
         xlab=paste("National forest ID:", unique(plots$forest)[i]))
  qqline(plots$ba[plots$forest==unique(plots$forest)[i]])
}
par(opar)

xyplot(ba ~ ht | forest,
       panel = function(x, y) {
         panel.xyplot(x, y)
         panel.loess(x,y, span=1)
       },
       data=plots)

sapply(overall.results, class)
sapply(results, class)

overall.results$interval <-
  factor(overall.results$interval, levels=unique(overall.results$interval))
results$interval <-
  factor(results$interval, levels=unique(results$interval))



means <- tapply(plots$ba, plots$forest, mean)

xyplot(contained ~ n | interval,
       panel = function(x, y) {
         panel.xyplot(x, y, type="l")
         panel.abline(h=0.95, col="darkgrey")
       },
       index.cond=list(c(1:3,7,4:6,8,9:14)),
       skip=c(rep(FALSE, 11), TRUE, rep(FALSE, 3)),
       data=overall.results)

xyplot(length ~ n | interval,
       panel = function(x, y) {
         panel.xyplot(x, y, type="l")
       },
       index.cond=list(c(1:3,7,4:6,8,9:14)),
       skip=c(rep(FALSE, 11), TRUE, rep(FALSE, 3)),
       data=overall.results)

xyplot(contained ~ length | interval,
       panel = function(x, y) {
         panel.xyplot(x, y, type="l")
         panel.abline(h=0.95, col="darkgrey")
       },
       index.cond=list(c(1:3,7,4:6,8,9:14)),
       skip=c(rep(FALSE, 11), TRUE, rep(FALSE, 3)),
       data=overall.results)

#########################################  BY FOREST


xyplot(length ~ n | interval,
       groups=forest,
       auto.key = list(points = TRUE, columns=2, title="Nearest Forest",
         corner=c(1,1), x=0.96, y=0.71),
       panel = function(x, y, subscripts, groups) {
         for (i in 1:length(levels(results$forest))) {
           panel.xyplot(x[groups[subscripts]==levels(results$forest)[i]],
                        y[groups[subscripts]==levels(results$forest)[i]],
                        col=trellis.par.get()$superpose.line$col[i],
                        type="l")
         }
       },
       skip=c(rep(FALSE, 11), TRUE, rep(FALSE, 3)),
       index.cond=list(c(1:3,7,4:6,8,9:14)),
       data=results)


xyplot(contained ~ n | interval,
       groups=forest,
       auto.key = list(points = TRUE, columns=2, title="Nearest Forest",
         corner=c(1,1), x=0.96, y=0.71),
       panel = function(x, y, subscripts, groups) {
         for (i in 1:length(levels(results$forest))) {
           panel.xyplot(x[groups[subscripts]==levels(results$forest)[i]],
                        y[groups[subscripts]==levels(results$forest)[i]],
                        col=trellis.par.get()$superpose.line$col[i],
                        type="l")
         }
         panel.abline(h=0.95, col="darkgrey")
       },
       skip=c(rep(FALSE, 11), TRUE, rep(FALSE, 3)),
       index.cond=list(c(1:3,7,4:6,8,9:14)),
       data=results)



xyplot(contained ~ length | interval,
       groups=forest,
       auto.key = list(points = TRUE, columns=2, title="Nearest Forest",
         corner=c(1,1), x=0.96, y=0.71),
       panel = function(x, y, subscripts, groups) {
         for (i in 1:length(levels(results$forest))) {
           panel.xyplot(x[groups[subscripts]==levels(results$forest)[i]],
                        y[groups[subscripts]==levels(results$forest)[i]],
                        col=trellis.par.get()$superpose.line$col[i],
                        type="l")
         }
         panel.abline(h=0.95, col="darkgrey")
       },
       skip=c(rep(FALSE, 11), TRUE, rep(FALSE, 3)),
       index.cond=list(c(1:3,7,4:6,8,9:14)),
       data=results)


require(e1071)

ks.statistic <- function(x, y, ...,
                 alternative = c("two.sided", "less", "greater"),
                 exact = NULL) 
  ks.test(x, y, ...,
          alternative = c("two.sided", "less", "greater"),
          exact = NULL)$statistic
  
by.forest <- aggregate(x = list(ks.ba=plots$ba),
                       by = list(forest=plots$forest),
                       FUN = ks.statistic,
                       y="pnorm")

results.temp <- results[results$interval %in% c("classical","linearized",
                                                "jackknife","studentized"),]
results.temp$interval <- factor(results.temp$interval)

my.trellis.pars <- trellis.par.get()
my.trellis.pars$superpose.line$lty <- 1:8
my.trellis.pars$superpose.line$col <- rep("black", 8)
trellis.par.set(my.trellis.pars)

xyplot(contained ~ n | forest,
       groups=interval,
       auto.key=list(columns=2, lines=TRUE, points=FALSE,
         title="Interval Style",
         text=c("Classical","Linearized","Jackknife","Studentized Bootstrap")),
       panel = function(x, y, subscripts, groups) {
         for (i in 1:length(levels(results.temp$interval))) {
           panel.xyplot(x[groups[subscripts]==
                          levels(results.temp$interval)[i]],
                        y[groups[subscripts]==
                          levels(results.temp$interval)[i]],
                        lty=trellis.par.get()$superpose.line$lty[i],
                        col=trellis.par.get()$superpose.line$col[i],
                        type="l")
         }
         panel.abline(h=0.95, col="darkgrey")
         panel.text(115, 0.962,
                    sprintf("%.3f",
                            with(by.forest,
                                 ks.ba[forest==
                                       results.temp$forest[subscripts][1]])))
       },
       index.cond=list(order(by.forest$ks.ba)),
       data=results.temp)


ks.normal.ratio.residual <- function(x, y) {  
  ks.test(residuals(lm(x ~ y - 1)), y="pnorm")$statistic
}

by.forest <- aggregate(x = list(ks.ba=plots$ba),
                       by = list(forest=plots$forest),
                       FUN = ks.normal.ratio.residual,
                       y = plots$ht)
