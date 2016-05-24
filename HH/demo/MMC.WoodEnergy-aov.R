## The energy data is from
##
## Milliken, G. A. and Johnson, D. E. (2002).
## Analysis of Messy Data Volume III: Analysis of Covariance,
## volume III, Chapman & Hall/CRC.

## There are three files for the MMC.WoodEnergy example.

## MMC.WoodEnergy-aov.R Data
##
## Description and data plots (Figure 8) and analysis of variance

## MMC.WoodEnergy.R
##
## MMC plots (Figure 9) in both S-Plus and R.

## MMC.WoodEnergy.s.R
##
## MMC plots as used in the MMC article.  These are only in S-Plus as
## they use the S-Plus multicomp adjust argument.  The glht function
## in the R multcomp package doesn't have an equivalent argument.

data(energy)

tpg.sl.original <- trellis.par.get("superpose.line")
tpg.sl <- as.data.frame(tpg.sl.original)
if.R(r=tpg.sl$col <- tpg.sl.original$col, s={})
tpg.sl <- rbind(tpg.sl[1:6,], tpg.sl[1:6,])
tpg.sl$lty <- rep(c(1,5,8), 4)
trellis.par.set("superpose.line", 
                if.R(r={tpg.sl.uc <- unclass(tpg.sl)
                        attr(tpg.sl.uc,"row.names") <- NULL
                        tpg.sl.uc},
                     s=tpg.sl))

tpg.ss.original <- trellis.par.get("superpose.symbol")
tpg.ss <- as.data.frame(tpg.ss.original)
if.R(r={tpg.ss$col <- tpg.ss.original$col
        tpg.ss$fill <- tpg.ss.original$fill},
     s={})
tpg.ss <- rbind(tpg.ss[1:6,], tpg.ss[1:6,])
trellis.par.set("superpose.symbol",
                if.R(r={tpg.ss.uc <- unclass(tpg.ss)
                        attr(tpg.ss.uc,"row.names") <- NULL
                        tpg.ss.uc},
                     s=tpg.ss))

xyplot(Energy ~ Moist | Wood, groups=Stove, data=energy,
       panel=panel.superpose,
       type="b", pch=16,
       layout=c(4,1),
       par.strip.text=list(cex=1.2),
       key=list(border=T, title="Stove",
         text=list(levels(energy$Stove), col=tpg.sl$col[1:3]),
         lines=tpg.sl[1:3,]))

energy.abline <-
xyplot(Energy ~ Moist | Wood, groups=Stove, data=energy,
       panel=function(x, y, subscripts, groups, ...) {
         panel.superpose(x=x, y=y, subscripts=subscripts, groups=groups, ...)
         for (i in 1:length(levels(groups)))
           panel.abline(lm(y[groups[subscripts]==levels(groups)[i]] ~
                           x[groups[subscripts]==levels(groups)[i]]),
                        col=tpg.sl$col[i],
                        lty=tpg.sl$lty[i])
         },
       type="p", pch=c("A","B","C"),
       scales=list(cex=1, alternating=1),
       between=list(x=c(1,1,1)),
       layout=c(4,1),
       par.strip.text=list(cex=1.2),
       key=list(border=T, title="Stove",
         text=list(levels(energy$Stove), col=tpg.sl$col[1:3]),
         lines=tpg.sl[1:3,],
         space="right"))
print(energy.abline, position=c(-.05,.4, 1,1))
## export.eps(h2("mmc/figure/AMDIII.ex.5.5.abline.eps"))



## MMC Figure 8
energy.am.abline <-
  xyplot(Energy ~ Moist | Wood, groups=Stove, data=energy,
         panel=function(x, y, subscripts, groups, ...) {
           panel.superpose(x=x, y=y, subscripts=subscripts, groups=groups, ...)
           for (i in 1:length(levels(groups))) {
             panel.abline(lm(y[groups[subscripts]==levels(groups)[i]] ~
                             x[groups[subscripts]==levels(groups)[i]]),
                          col=tpg.sl$col[i],
                          lty=tpg.sl$lty[i])
             panel.abline(v=14, lty=2)
           }
         },
         type="p", pch=c("A","B","C"),
         scales=list(cex=1, alternating=1),
         between=list(x=c(1,1,1)),
         layout=c(4,1),
         par.strip.text=list(cex=1.2),
         key=list(border=T, title="Stove",
           text=list(levels(energy$Stove),col=tpg.sl$col[1:3]),
           lines=tpg.sl[1:3,],
           space="right"))
print(energy.am.abline, position=c(-.05,.3, 1,1))
## export.eps(h2("mmc/figure/AMDIII.ex.5.5.am.abline.eps"))

energy.am.abline.2 <- energy.am.abline
if.R(s={
  energy.am.abline.2$key <- NULL
  energy.am.abline.2$between$x[] <- 4
},r={
  energy.am.abline.2 <- update(energy.am.abline.2, legend=NULL)
})
print(energy.am.abline.2, position=c(0,.6, 1,1))
## export.eps(h2("mmc/figure/AMDIII.ex.5.5.am.2.abline.eps"))




tpg.sl <- as.data.frame(tpg.sl.original)
if.R(r=tpg.sl$col <- tpg.sl.original$col, s={})
tpg.sl <- rbind(tpg.sl[1:3,], tpg.sl[1:3,], tpg.sl[1:3,], tpg.sl[1:3,])
trellis.par.set("superpose.line", 
                if.R(r={tpg.sl.uc <- unclass(tpg.sl)
                        attr(tpg.sl.uc,"row.names") <- NULL
                        tpg.sl.uc},
                     s=tpg.sl))

energy.aov.4 <- aov(Energy ~ Moist + Stove*Wood + Moist:Stove:Wood,
                    data=energy)
anova(energy.aov.4)


## Cheung and Chan method
energy.nsize <- tapply(energy$Energy, energy[,c("Wood","Stove")], length)
energy.nsize
date() ## tpmc takes several minutes
energy.tpmc <-
  try(
  tpmc(NGROUP=4,
       NK=3,
       DF=energy.aov.4$df.residual,
       NSIZE=energy.nsize)
      )
if (class(energy.tpmc)==if.R(s="Error", r="try-error"))
  energy.tpmc <- 2.925457  ## value calculated on machine where tpmc is compiled
date()
