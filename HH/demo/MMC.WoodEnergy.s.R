## This file follows MMC.WoodEnergy-aov.R


if.R(r=
     cat(
"The glht function in the multcomp package in R does not
have an adjust argument.  See file MMC.WoodEnergy.R for
an alternate way to produce all the figures in the article.\n")
     ,s={


## simulation method
energy.multicomp <- multicomp(energy.aov.4, focus="Stove",
                              adjust=list(Wood=c(
                                            "Osage Orange",
                                            "Red Oak",
                                            "Black Walnut",
                                            "White Pine")),
                              plot=FALSE)
oldClass(energy.multicomp) <- c("multicomp.hh", "multicomp")
energy.multicomp$crit.point
## [1] 2.957413

## Tukey method
energy.multicomp <- multicomp(energy.aov.4, focus="Stove",
                              adjust=list(Wood=c(
                                            "Osage Orange",
                                            "Red Oak",
                                            "Black Walnut",
                                            "White Pine")),
                              plot=FALSE,
                              method="tukey", valid.check=FALSE)
oldClass(energy.multicomp) <- c("multicomp.hh", "multicomp")
energy.multicomp$crit.point
##    tukey 
## 2.937812

##
##
energy.multicomp <- multicomp(energy.aov.4, focus="Stove",
                              adjust=list(Wood=c(
                                            "Osage Orange",
                                            "Red Oak",
                                            "Black Walnut",
                                            "White Pine")),
                              crit.point=energy.tpmc,
                              method="tpmc",
                              plot=FALSE)
oldClass(energy.multicomp) <- c("multicomp.hh", "multicomp")
energy.multicomp$crit.point
## [1] 2.925457

old.par <- par(oma=c(0,3,0,0))
plot(energy.multicomp, col.signif='red', lty.signif=1)
par(old.par)
## export.eps(h2("mmc/figure/AMDIII.ex.5.5.multicomp-original.eps"))

energy.multicomp$method <- "tpmc"
energy.multicomp <- multicomp.label.change(energy.multicomp, ".adj1", ".OsgOr")
energy.multicomp <- multicomp.label.change(energy.multicomp, ".adj2", ".RdOak")
energy.multicomp <- multicomp.label.change(energy.multicomp, ".adj3", ".BkWal")
energy.multicomp <- multicomp.label.change(energy.multicomp, ".adj4", ".WPine")
energy.multicomp <- multicomp.reverse(energy.multicomp)
old.par <- par(oma=c(0,3,0,0))
plot(energy.multicomp, col.signif='red', lty.signif=1)
par(old.par)
## export.eps(h2("mmc/figure/AMDIII.ex.5.5.multicomp.eps"))
dimnames(energy.multicomp$lmat)[[1]]

energy.multicomp.order <- energy.multicomp
energy.multicomp.order <-
  multicomp.order(energy.multicomp.order,
                  sort.order=c(2,1,3, 5,6,4, 9,8,7, 10,12,11))
old.par <- par(oma=c(0,3,0,0))
plot(energy.multicomp.order, col.signif='red', lty.signif=1)
par(old.par)
## export.eps(h2("mmc/figure/AMDIII.ex.5.5.multicomp.order.eps"))


energy.new <-
  data.frame(Moist=14,
             Wood=factor(rep(levels(energy$Wood),c(3,3,3,3)),
               levels=levels(energy$Wood)),
             Stove=rep(levels(energy$Stove),4))
energy.new$Energy.adj <-
  predict(energy.aov.4, newdata=energy.new)
energy.new

energy.am.mmc <- list()
for (i in levels(energy.new$Wood))
  energy.am.mmc[[i]] <-
  multicomp.mmc.mean(energy.new$Stove[1:3],
                     energy.nsize[i,],
                     energy.new[energy.new$Wood==i,"Energy.adj"],
                     sqrt(anova(energy.aov.4)[6,"Mean Sq"]),
                     ylabel="adjusted Energy",
                     focus="Stove",
                     plot=FALSE,
                     method="tpmc",
                     crit.point=energy.tpmc)

old.dev <- dev.cur()

trellis.device.hh.color(width=11, height=2.75)
old.par <- par(mfrow=c(1,4), mar=c(5,4,4,4)+.1)
for (i in levels(energy$Wood))
  plot(energy.am.mmc[[i]], ry=c(1.3,7.2), main=i, main2="")
## export.eps(h2("mmc/figure/AMDIII.ex.5.5.am.mmc.eps"))
par(old.par)

old.par <- par(mfrow=c(1,4), mar=c(5,4,4,4)+.1)
## MMC Figure 9a
for (i in levels(energy$Wood))
  plot(energy.am.mmc[[i]], ry=c(1.3,7.2), main=i, main2="",
       xlab="", ylab="", focus="")
## export.eps(h2("mmc/figure/AMDIII.ex.5.5.am.mmc.n.eps"))
par(old.par)
dev.set(old.dev) ## dev.off()

trellis.device.hh.color(width=11, height=2.5)
old.par <- par(mfrow=c(1,4), mar=c(5,4,4,4)+.1)
for (i in levels(energy$Wood))
  plot(energy.am.mmc[[i]],  main=i, main2="")
## export.eps(h2("mmc/figure/AMDIII.ex.5.5.am4.mmc.eps"))
par(old.par)
dev.set(old.dev) ## dev.off()

## back to the standard-size graphics device
old.par <- par(mfrow=c(2,2), mar=c(5,4,4,6)+.1)
for (i in levels(energy$Wood))
  plot(energy.am.mmc[[i]],  main=i, main2="")
## export.eps(h2("mmc/figure/AMDIII.ex.5.5.am4-2x2.mmc.eps"))
par(old.par)

old.par <- par(mfrow=c(2,2), mar=c(5,4,4,6)+.1)
## MMC Figure 9b
for (i in levels(energy$Wood))
  plot(energy.am.mmc[[i]],  main=i, main2="",
       xlab="", ylab="", focus="")
## export.eps(h2("mmc/figure/AMDIII.ex.5.5.am4-2x2.mmc.n.eps"))
par(old.par)

energy.am.mca <- energy.am.mmc[[levels(energy$Wood)[1]]]$mca
for (i in levels(energy$Wood)[-1]) {
  energy.am.mca$table <- rbind(energy.am.mca$table, energy.am.mmc[[i]]$mca$table)
  energy.am.mca$height <- c(energy.am.mca$height, energy.am.mmc[[i]]$mca$height)
}
## MMC Figure 9c
old.par <- par(mar=c(5,10,4,2)+.1)
plot(energy.am.mca, col.signif='red', lty.signif=1)
segments(-1.6, x2=2.8, y1=12.5, y2=12.5, lty=2, xpd=T)
segments(-1.6, x2=2.8, y1=15.5, y2=15.5, lty=2, xpd=T)
segments(-1.6, x2=2.8, y1=18.5, y2=18.5, lty=2, xpd=T)
axis(levels(energy$Wood), side=2, las=1, line=6,
     at=c(21.5, 18.5, 15.5, 12.5)-1.5,
     ticks=F)
## export.eps(h2("mmc/figure/AMDIII.ex.5.5.am.mca.eps"))
par(old.par)

energy.aov.4a <- aov(Energy ~ Moist*(Wood/Stove),
                    data=energy)
anova(energy.aov.4a)

trellis.par.set("superpose.line", tpg.sl.original)
trellis.par.set("superpose.symbol", tpg.ss.original)

})
