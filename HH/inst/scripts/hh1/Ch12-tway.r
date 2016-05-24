## xx rhiz-clov-mmc.s, adjusted mmc for R not yet

## display.s
## display.r

data(display)

interaction2wt(time ~ panel + emergenc, data=display,
               par.strip.text=list(cex=1.1),
               main.cex=1.6,
               scales=list(x=list(cex=.9), y=list(cex=.9, alternating=FALSE)))
## export.eps(hh("tway/figure/display.eps"))



## display$panel has polynomial contrasts as part of "position" for graphing.
## glht in R needs treatment contrasts.
contrasts(display$panel) <- contr.treatment(3)


## emergenc fixed
displayf.aov <- aov(time ~ emergenc * panel, data=display)
anova(displayf.aov)

## emergenc random
displayr.aov <- aov(time ~  panel + Error(emergenc + panel:emergenc),
                    data=display)
summary(displayr.aov)

## multiple comparisons, Tukey method
tmp <- glht(displayf.aov, linfct=mcp(
                            panel="Tukey",
                            `interaction_average` = TRUE,
                            `covariate_average` = TRUE))
plot(tmp)
confint(tmp)
## export.eps(hh("tway/figure/display-mcdiff.eps"))

tmp <- glht(displayf.aov, linfct=mcp(panel="Means",
                            `interaction_average` = TRUE,
                            `covariate_average` = TRUE))
plot(tmp)
confint(tmp)
## export.eps(hh("tway/figure/display-mcmean.eps"))

## multiple comparisons, mmc plot, Tukey method
displayf.mmc <- mmc(displayf.aov, focus="panel")
plot(displayf.mmc, x.offset=1.5)
displayf.mmc
## export.eps(hh("tway/figure/display-mmcdiff.eps"))

## compare with SAS display
tmp <- glht(displayf.aov, linfct=mcp(
                            panel="Tukey",
                            `interaction_average` = TRUE,
                            `covariate_average` = TRUE))
print(qtukey(.95, 3, 12)/sqrt(2))            ##
print(attr(confint(tmp)$confint, "calpha"))  ## R critical value
print(attr(confint(tmp)$confint, "calpha") * sqrt(2))  ## SAS: Critical Value of Studentized Range

ms.res <- summary(displayf.aov)[[1]]["Residuals","Mean Sq"]
sqrt(2*ms.res/8) * attr(confint(tmp)$confint, "calpha")  ## minimum significant difference
## there are 8 observations in each of the three groups



## rhiz-read.s         # rhiz, read data
## Erdman reported the alfalfa data in Table 1.
data(rhiz.alfalfa)
rhiz.alfalfa <- rhiz.alfalfa ## local copy

## Erdman reported the clover data in Table 3.
data(rhiz.clover)
rhiz.clover <- rhiz.clover ## local copy





## rhiz-bwplot.ti.r    # rhiz, Fig 11.2 11.3 11.4 ## preferred
## vertical boxplots conditioned on combination
## our preference

print(position = c(0, .47, 1, .94), more = TRUE,  # top
bwplot(Npg ~ strain | comb, data=rhiz.alfalfa,
         main="alfalfa", layout=c(2,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.75))
)
print(position = c(0, 0, 1, .47), more = TRUE,  # bottom
bwplot(Npg ~ strain | comb, data=rhiz.clover,
         main="clover", layout=c(2,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.75))
)
print(position=c(0, .95, 1, 1), more=FALSE,
      xyplot(0 ~ 0, panel=function(...){},
             main=list("Nitrogen per gram", cex=1.2),
             ylab=expression(""),
             xlab="", scales=list(draw=FALSE),
             par.settings = list(axis.line = list(col = "transparent")))
)
## export.eps(hh("tway/figure/rhnpg.ti.eps"))


print(position = c(0, .47, 1, .94), more = TRUE,  # top
bwplot(nitro ~ strain | comb, data=rhiz.alfalfa,
         main="alfalfa", layout=c(2,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.75))
)
print(position = c(0, 0, 1, .47), more = TRUE,  # bottom
bwplot(nitro ~ strain | comb, data=rhiz.clover,
         main="clover", layout=c(2,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.75))
)
print(position=c(0, .95, 1, 1), more=FALSE,
      xyplot(0 ~ 0, panel=function(...){},
             main=list("Nitrogen", cex=1.2),
             ylab=expression(""),
             xlab="", scales=list(draw=FALSE),
             par.settings = list(axis.line = list(col = "transparent")))
)
## export.eps(hh("tway/figure/rhn.ti.eps"))


print(position = c(0, .47, 1, .94), more = TRUE,  # top
bwplot(weight ~ strain | comb, data=rhiz.alfalfa,
         main="alfalfa", layout=c(2,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.75))
)
print(position = c(0, 0, 1, .47), more = TRUE,  # bottom
bwplot(weight ~ strain | comb, data=rhiz.clover,
         main="clover", layout=c(2,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.75))
)
print(position=c(0, .95, 1, 1), more=FALSE,
      xyplot(0 ~ 0, panel=function(...){},
             main=list("Weight", cex=1.2),
             ylab=expression(""),
             xlab="", scales=list(draw=FALSE),
             par.settings = list(axis.line = list(col = "transparent")))
)
## export.eps(hh("tway/figure/rhw.ti.eps"))

## The HH text (page 357) suggests three alternate arrangements for
## the boxplots, and references three files for the details.
##
## horizontal boxplots conditioned on strain
## rhiz-bwplot.r       # rhiz, (Fig 11.2 11.3 11.4) alternate
## rhiz-bwplot.s       # rhiz, (Fig 11.2 11.3 11.4) alternate
##
## vertical boxplots conditioned on strain
## rhiz-bwplot.t.r     # rhiz, (Fig 11.2 11.3 11.4) alternate
## rhiz-bwplot.t.s     # rhiz, (Fig 11.2 11.3 11.4) alternate
##
## horizontal boxplots conditioned on combination
## rhiz-bwplot.i.r     # rhiz, (Fig 11.2 11.3 11.4) alternate
## rhiz-bwplot.i.s     # rhiz, (Fig 11.2 11.3 11.4) alternate
##
## All three are collected together in this directory in file
## hh("scripts/Ch12-rhiz-bwplot-alternate.r")
cat("See file Ch12-rhiz-bwplot-alternate.r\n")



## rhiz-alf-aov.s      # rhiz, Fig 11.5 11.6 Table 11.5

rhiz.alfalfa.aov <- aov(Npg ~ strain * comb, data=rhiz.alfalfa)
summary(rhiz.alfalfa.aov)

alf.means <- model.tables(rhiz.alfalfa.aov, type="means", se=TRUE)
alf.means$tables$strain
alf.means$n["strain"]
alf.means$se$strain

alfalfa.mca <- glht(rhiz.alfalfa.aov, linfct=mcp(strain="Tukey",
                                        `interaction_average` = TRUE,
                                        `covariate_average` = TRUE))
old.omd <- par(omd=c(.15,1, 0,1))
plot(alfalfa.mca)
par(old.omd)
print(confint(alfalfa.mca))
## export.eps(hh("tway/figure/alfmc.ps"))


old.fin <- par()$fin
par(fin=c(old.fin[1], 3))

plot(y=c(1.01,1,.99,1,1,1), x=alf.means$tables$strain,
     pch=16,
     yaxt="n", ylim=c(.85,1.10), ylab="",
     xaxt="n", xlab="")
rug(alf.means$tables$strain, ticksize=-.1)
axis(3)
mtext("Npg", 3, line=3)
lines(x=c(29.7,31.5), y=c(.95,.95))
lines(x=c(30.5,32.5), y=c(.90,.90))
strain.labels <- paste(dimnames(alf.means$tables$strain)[[1]],
                       format(alf.means$tables$strain, digits=6),
                       sep=": ")

axis(1, at=alf.means$tables$strain[c(2,4, 6)], labels=strain.labels[c(2,4,  6)], line=1, tick=FALSE, adj=.45)
axis(1, at=alf.means$tables$strain[c(    5  )], labels=strain.labels[c(    5  )], line=2, tick=FALSE, adj=.45)
axis(1, at=alf.means$tables$strain[3         ], labels=strain.labels[3         ], line=2, tick=FALSE, adj=.45)
axis(1, at=alf.means$tables$strain[1         ], labels=strain.labels[1         ], line=3, tick=FALSE, adj=.4)

par(fin=old.fin)
## export.eps(hh("tway/figure/alfmeans.eps"))




## rhiz-clov-aov.s     # rhiz, Fig 11.7 Table 11.6 11.7 11.8 11.9
rhiz.clover.aov <- aov(Npg ~ strain * comb, data=rhiz.clover)
summary(rhiz.clover.aov)

clov.means <- model.tables(rhiz.clover.aov, type="means", se=TRUE)

## tway/figure/clovint.eps
old.mar <- par(mar=c(0,2,0,0)+par("mar"))
interaction.plot(rhiz.clover$strain, rhiz.clover$comb, rhiz.clover$Npg,
                 lty=c(1,6), cex=1.3)
par(old.mar)
## export.eps(hh("tway/figure/clovint.eps"))


## tway/figure/clovint2wt.eps
position(rhiz.clover$comb) <- c(2.5,4.5)
interaction2wt(Npg ~ strain + comb, data=rhiz.clover,
               par.strip.text=list(cex=1.4),
               main.cex=1.6,
               scales=list(x=list(cex=.6), y=list(cex=.9, alternating=FALSE)),
               key.in=list(cex=.8, size=1, cex.title=1.2, between=.5),
               xlim=c(.5, 6.5))
## export.eps(hh("tway/figure/clovint2wt.eps"))
## export.eps(hh("tway/figure/clovint2wt.color.eps"))



rhiz.clover.nest.aov <- aov(Npg ~ comb/strain, data=rhiz.clover)
summary(rhiz.clover.nest.aov)

names(coef(rhiz.clover.nest.aov))
summary(rhiz.clover.nest.aov,
        split=list("comb:strain"=
          list(clover=c(1,3,5,7,9),
               "clover+alfalfa"=c(2,4,6,8,10))))

## Look at the contrasts, their generated dummy variables,
## and their regression coefficients.
## Abbreviate their names for presentation.
tmp <- abbreviate(names(coef(rhiz.clover.nest.aov)))
tmp

contrasts(rhiz.clover$comb)
contrasts(rhiz.clover$strain)

cnx <- aov(Npg ~ comb/strain, data=rhiz.clover, x=TRUE)$x
dimnames(cnx)[[2]] <- tmp
cnx

cnxb <- round(coef(summary.lm(rhiz.clover.nest.aov)), 3)
dimnames(cnxb)[[1]] <- tmp
cnxb




## rhiz-alf-mmc.s
old.omd <- par(omd=c(0,.85 ,0,1))
## the above value of omd sets fig to the invalid value
## > par()$fig
## [1] -0.08823529  1.08823529  0.00000000  1.00000000
par(fig=c(0,1,0,1)) ## this fig restores it.
alf.mmc <- mmc(rhiz.alfalfa.aov, focus="strain")
plot(alf.mmc, x.offset=.3, ry=c(29, 33))
## export.eps(hh("tway/figure/alfalfa.mmc.eps"))
## this plot is hard to read, we need to break the ties and retain the order
plotMatchMMC(alf.mmc$mca)
print(alf.mmc)
## export.eps(hh("tway/figure/alfalfa.mmc-tiebreak.eps"))

## construct an orthogonal set of interpretable contrasts
zapsmall(alf.mmc$none$lmat[2:6,]) ## strain rows
##                              c 15 12  1 10  7
alf.comp <- cbind("1,7,10-c"=c(-3, 0, 0, 1, 1, 1),
                  "1,10-7"  =c( 0, 0, 0, 1, 1,-2),
                  "1-10"    =c( 0, 0, 0, 1,-1, 0),
                  "15-12"   =c( 0, 1,-1, 0, 0, 0),
            "1,7,10,c-12,15"=c( 1,-2,-2, 1, 1, 1))
dimnames(alf.comp)[[1]] <- dimnames(alf.mmc$none$lmat)[[2]]
alf.comp
alf.mmc <- mmc(rhiz.alfalfa.aov, focus="strain",
                         focus.lmat=alf.comp)
plot(alf.mmc, x.offset=.3, ry=c(29, 33))  ### these are at the wrong heights
## export.eps(hh("tway/figure/alfalfa.lmat.mmc.eps"))
## this plot is hard to read, we need to break the ties and retain the order
plotMatchMMC(alf.mmc$lmat, col.signif="blue")
print(alf.mmc)
## export.eps(hh("tway/figure/alfalfa.lmat.mmc-tiebreak.eps"))
par(old.omd)


## rhiz-clov-mmc.s

## rhiz.clover.nest.aov <- aov(Npg ~ comb/strain, data=rhiz.clover)
summary(rhiz.clover.nest.aov)
summary(rhiz.clover.nest.aov,
        split=list("comb:strain"=
          list(clover=c(1,3,5,7,9),
               "clover+alfalfa"=c(2,4,6,8,10))))
model.tables(rhiz.clover.nest.aov, type="means", se=TRUE)
## This gives the anova table we are interested,
## R glht won't work with a nested model.

#####  see HH-R.package/HH/demo/MMC.WoodEnergy.R   #####
#####  see also 8004.s09/0212/ironpot.s            #####
#####  hh/dsgn/code/filmcoat3.s                    #####
#####  hh/dsgntwo/answer/bean4.s                   #####

## use the crossed model from rhiz-clov-aov.s with the adjust argument
rhiz.clover.aov <- aov(Npg ~ strain * comb, data=rhiz.clover)
summary(rhiz.clover.aov)
model.tables(rhiz.clover.aov, type="means", se=TRUE)

## compare simple effects with clover+alfalfa to simple effects with clover
## force these graphs to have the same xlim and ylim as the previous graphs
## simple effects with clover
## rhiz.clover.aov <- aov(Npg ~ strain * comb, data=rhiz.clover)
rhiz.clover.comb.aov <- list()
rhiz.clover.comb.aov[["clover"]] <-
  aov(Npg ~ strain, data=rhiz.clover, subset=(comb=="clover"))
rhiz.clover.comb.aov[["clover+alfalfa"]] <-
  aov(Npg ~ strain, data=rhiz.clover, subset=(comb=="clover+alfalfa"))

## sapply(rhiz.clover.comb.aov, summary)

ResidMS <- function(x)
  summary(x)[[1]]["Residuals","Mean Sq"]
ResidMSAvg <- mean(sapply(rhiz.clover.comb.aov, ResidMS))

## individually scaled
clover.str.clov.mmc <- lapply(rhiz.clover.comb.aov, mmc,
                              linfct=mcp(strain="Tukey"),
                              calpha=qtukey(1-.05/2, 6, 48)/sqrt(2))
old.omd <- par(omd=c(0,.87,0,1))
for (i in levels(rhiz.clover$comb)) {
  plot(clover.str.clov.mmc[[i]], ry=c(16,40), main=i, main2="")
  plotMatchMMC(clover.str.clov.mmc[[i]]$mca, col.signif=8, lty.signif=1, main=i)
}
par(old.omd)

## scaled to common value
calpha <- qtukey(1-.05/2, 6, 48)/sqrt(2)
clover.commonstrMS.clov.mmc <- list()
clover.commonstrMS.clov.mmc[["clover"]] <-
  mmc(rhiz.clover.comb.aov[["clover"]], linfct=mcp(strain="Tukey"),
      calpha=sqrt(ResidMSAvg/ResidMS(rhiz.clover.comb.aov[["clover"]])) *
      calpha)
clover.commonstrMS.clov.mmc[["clover+alfalfa"]] <-
  mmc(rhiz.clover.comb.aov[["clover+alfalfa"]], linfct=mcp(strain="Tukey"),
      calpha=sqrt(ResidMSAvg/ResidMS(rhiz.clover.comb.aov[["clover+alfalfa"]])) *
      calpha)
cat("\nclover:\n",
    paste("calpha = ", round(calpha, 2),
          ", ResidMSAvg = ", round(ResidMSAvg, 2),
          ", ResidMS ", "clover", " = ",
          round(ResidMS(rhiz.clover.comb.aov[["clover"]]), 2),
          "\n", sep=""),
    "Quantile below is sqrt(ResidMSAvg/ResidMS) * calpha\n",
    "stderr below is based on ResidMS, lower and upper bounds are based on ResidMSAvg\n")
print(clover.commonstrMS.clov.mmc[["clover"]])
cat("\nclover+alfalfa:\n",
    paste("calpha = ", round(calpha, 2),
          ", ResidMSAvg = ", round(ResidMSAvg, 2),
          ", ResidMS ", "clover+alfalfa", " = ",
          round(ResidMS(rhiz.clover.comb.aov[["clover+alfalfa"]]), 2),
          "\n", sep=""),
    "Quantile below is sqrt(ResidMSAvg/ResidMS) * calpha\n",
    "stderr below is based on ResidMS, lower and upper bounds are based on ResidMSAvg\n")
print(clover.commonstrMS.clov.mmc[["clover+alfalfa"]])

old.omd <- par(omd=c(0,.87,0,1))
for (i in levels(rhiz.clover$comb)) {
  plot(clover.commonstrMS.clov.mmc[[i]], ry=c(16,40), main=i,
       main2=paste("Average Residual MS =", round(ResidMSAvg, 2)))
  plotMatchMMC(clover.commonstrMS.clov.mmc[[i]]$mca, col.signif=8, lty.signif=1, main=i, adjusted=TRUE)
}
par(old.omd)

## construct an orthogonal basis set of contrasts

zapsmall(clover.str.clov.mmc$clover$none$lmat) ## strain rows
##                               5  1  c  7 13  4
clover.comp <- cbind("5-rest"=c( 5,-1,-1,-1,-1,-1),
                   "1-4.7.13"=c( 0, 3, 0,-1,-1,-1),
                 "c-1.4.7.13"=c( 0,-1, 4,-1,-1,-1),
                       "7-13"=c( 0, 0, 0, 1,-1, 0),
                     "4-7.13"=c( 0, 0, 0, 1, 1,-2))
dimnames(clover.comp)[[1]] <- dimnames(clover.str.clov.mmc$clover$none$lmat)[[2]]

clover.str.clov.mmc$clover <- mmc(rhiz.clover.comb.aov$clover,
                                  linfct=mcp(strain="Tukey"),
                                  calpha=qtukey(1-.05/2, 6, 48)/sqrt(2),
                                  focus.lmat=clover.comp)
plot(clover.str.clov.mmc$clover, ry=c(16,40), main=i, main2="")
plotMatchMMC(clover.str.clov.mmc$clover$lmat, col.signif=8, lty.signif=1, main=i)





## plasma.s
data(plasma)
plasma <- plasma ## local copy

old.par <- par(mfrow=c(2,1))
interaction.plot(plasma$id, plasma$time, plasma$plasma, fixed=TRUE, cex=1.3)
interaction.plot(plasma$time, plasma$id, plasma$plasma, cex=1.3)
par(old.par)
##trellis.device(postscript,file="tway/figure/plasmaint.eps",horizontal=TRUE)
## I see anomolies for id=3 at 8pm and for id=6 at 11am.

position(plasma$time) <- c(1,3,5,7,9)+.5
print(position=c(0,0,.97,1),
      interaction2wt(plasma ~ time + id, data=plasma,
                     par.strip.text=list(cex=1.4),
                     main.cex=1.6,
                     scales=list(x=list(cex=.9), y=list(cex=.9, alternating=FALSE)))
      )
## export.eps(hh("tway/figure/plasma2wt.eps"))

plasma.aov <- aov(plasma ~ id + time, data=plasma)
summary(plasma.aov)


## workstation.s
data(workstation)

workstation.aov <- aov(devices ~ method / station,
                       data=workstation, qr=TRUE)
summary(workstation.aov)
model.tables(workstation.aov, se=TRUE)

bwplot(devices ~ station | method, data=workstation,
       ylab=list(cex=1.4),
       xlab=list("station %in% method", cex=1.4),
       strip=function(..., strip.names)
       strip.default(..., strip.names=c(TRUE, TRUE)),
       par.strip.text=list(cex=1.4),
       scales=list(x=list(cex=1),y=list(cex=1.2)),
       layout=c(3,1))
## export.eps(hh("tway/figure/workstation.eps"))


## feed.s
data(feed)

print(position=c(0,0,.96,1),
      interaction2wt(retained ~ supp + temp, data=feed,
                     ##par.strip.text=list(cex=1.4),
                     main.cex=1.6,
                     scales=list(x=list(cex=.9), y=list(cex=.9, alternating=FALSE)))
      )
## export.eps(hh("tway/figure/feed-i2wt.eps"))
## export.eps(hh("tway/figure/feed-i2wt-color.eps"))




feed.int.aov <- aov(retained ~ temp * supp, data=feed)
anova(feed.int.aov)
model.tables(feed.int.aov, type="means", se=TRUE)

old.par <- par(mfrow=c(2,1))
interaction.plot(feed$temp, feed$supp, feed$retained, fixed=TRUE, cex=1.3)
interaction.plot(feed$supp, feed$temp, feed$retained, fixed=TRUE, cex=1.3)
par(old.par)
## export.eps(hh("tway/figure/feed.eps"))

feed.aov <- aov(retained ~ temp + supp, data=feed)
anova(feed.aov)
model.tables(feed.aov, type="means", se=TRUE)

tapply(feed$retained, feed[,c("supp", "temp")], mean)


summary(feed.aov, split=
        list(temp=list(linear=1,quadratic=2),
             supp=list(linear=1,quadratic=2,cubic=3)))


