## This file follows MMC.WoodEnergy-aov.R
if (!exists("energy.aov.4"))
  stop('Please run demo("MMC.WoodEnergy-aov", package="HH") before running this demo.',
       call.=FALSE)

## multicomp main effect for Stove
if.R(r={
  energy.glht <- glht(energy.aov.4, focus="Stove",
                      linfct=mcp(Stove="Tukey",
                      `interaction_average`=TRUE, `covariate_average`=TRUE))
  ## warning from glht.matrix about not-estimable coefficient may be ignored,
  ## tmp <- model.matrix(energy.aov.4)
  ## cbind(rowSums(tmp[,14:25]), tmp[,2])  ## three-way interaction and covariate are equal
  Stove.means <- model.tables(energy.aov.4, type="means",
                              cterms="Stove")$tables$Stove
  height.mca <- Stove.means %*% abs(t(contrMat(Stove.means, "Tukey")))
  energy.multicomp <- as.multicomp(energy.glht, focus="Stove",
                                   lmat.rows=3:4, height=height.mca)
},s={
  energy.multicompS <- multicomp(energy.aov.4, focus="Stove",
                                plot=FALSE,
                                method="tukey", valid.check=FALSE)
  oldClass(energy.multicompS) <- c("multicomp.hh", "multicomp")
  zapsmall(energy.multicompS$lmat)
  energy.multicompS
})


## separate ANOVA for each Wood
energy.aov.4W <- list()
for (i in levels(energy$Wood)) {
  energy.aov.4W[[i]] <- aov(Energy ~ Moist + Stove + Moist:Stove,
                      data=energy, subset=(Wood==i))
  print(anova(energy.aov.4W[[i]]))
}


## Separate multicomp for each Wood, each with its own
## critical values.
## calculate, print, plot
if.R(r={
  energy.mca.4W <- list()
  par(mfrow=c(2,2))
  for (i in levels(energy$Wood)) {
    energy.mca.4W[[i]] <-
      glht(energy.aov.4W[[i]],
           linfct=mcalinfct(energy.aov.4W[[i]], "Stove"))
    print(confint(energy.mca.4W[[i]]))
    plot(energy.mca.4W[[i]], xlim=c(-3,3), main=i, ylim=c(.5,3.5))
  }
  par(mfrow=c(1,1))
},s={
  energy.mca.4W <- list()
  par(mfrow=c(2,2))
  for (i in levels(energy$Wood)) {
    energy.mca.4W[[i]] <-
      multicomp(energy.aov.4W[[i]], method="tukey", focus="Stove")
    oldClass(energy.mca.4W[[i]]) <- c("multicomp.hh", "multicomp")
    print(energy.mca.4W[[i]])
    frame()
    par(new=TRUE)
    plot(energy.mca.4W[[i]], xlim=c(-3,3), main=i, comparisons.per.page=6)
  }
  par(mfrow=c(1,1))
})


## Separate calculated multicomp for each Wood, with forced adjustment
## for Moist=14 and common critical value.

## Common contrast matrix for all Wood types.
if.R(r={
  mca.14 <- mcalinfct(energy.aov.4W[[1]], "Stove")
  non.zero <- mca.14[,5:6] != 0
  mca.14[,5:6][non.zero] <- 14 * sign(mca.14[,5:6][non.zero])
},s={
  mca.14 <- zapsmall(energy.mca.4W[[1]]$lmat)
  non.zero <- mca.14[6:8,] != 0
  mca.14[6:8,][non.zero] <- 14 * sign(mca.14[6:8,][non.zero])
}
     )


## calculate
energy.mca.4W <- list()
for (i in levels(energy$Wood)) {
  energy.mca.4W[[i]] <-
    if.R(r={
      glht(energy.aov.4W[[i]],
           linfct=mca.14)
    },s={
      tmp <- multicomp(energy.aov.4W[[i]],
                       focus="Stove",
                       lmat=mca.14, comparisons="none",
                       crit.point=energy.tpmc)
      oldClass(tmp) <- c("multicomp.hh", "multicomp")
      tmp
    })
}

## display full multicomp or confint printed output
if.R(s=
     energy.mca.4W
     ,r=
     lapply(energy.mca.4W, confint, calpha=energy.tpmc)
     )

## display just the tables
if.R(s=
     lapply(energy.mca.4W, `[[`, "table")
     ,r=
     lapply(energy.mca.4W, function(x) confint(x, calpha=energy.tpmc)$confint[,])
     )


## plot multicomp or confint plots for all four Woods on one page
par(mfrow=c(2,2))
for (i in levels(energy$Wood))
  if.R(r=
       plot(confint(energy.mca.4W[[i]], calpha=energy.tpmc),
            xlim=c(-3,3), main=i, xlab="", ylim=c(.5,3.5))
       ,
       s={
         frame()
         par(new=TRUE)
         plot(energy.mca.4W[[i]], xlim=c(-3,3), main=i,
              comparisons.per.page=6, xlabel.print=FALSE)
       }
       )
par(mfrow=c(1,1))


## Do all four as MMC plots
## calculate
energy.mmc.4W <- list()
for (i in levels(energy$Wood)) {
  energy.mmc.4W[[i]] <-
    if.R(r={
      mmc(energy.aov.4W[[i]], calpha=energy.tpmc,
          linfct=mca.14, lmat.rows=3:4, focus="Stove")
    },s={
      multicomp.mmc(energy.aov.4W[[i]],
                    focus="Stove",
                    lmat=mca.14, comparisons="none",
                    lmat.rows=3:5,
                    crit.point=energy.tpmc,
                    plot=FALSE)
    })
}

## MMC Figure 9a, MMC plots in common scale
par(mfrow=c(2,2))
for (i in levels(energy$Wood))
  plot(energy.mmc.4W[[i]], ry=c(1.3,7.2), main=i, main2="",
       xlab="", ylab="", focus="",
       col.lmat.signif='red')
par(mfrow=c(1,1))

## MMC Figure 9b, MMC plots in individual scales
par(mfrow=c(2,2))
for (i in levels(energy$Wood))
  plot(energy.mmc.4W[[i]], main=i, main2="",
       xlab="", ylab="", focus="",
       col.lmat.signif='red')
par(mfrow=c(1,1))


## MMC Figure 9c, multicomp or confint plots in common scale
par(mfrow=c(4,2))
for (i in levels(energy$Wood))
  if.R(r={
    old.mar <- par(mar=c(4.1, 4.1, 1.1, 2.1))
    plot(x=0:1, y=0:1, type="n", bty="n",
         xaxt="n", yaxt="n", xlab="", ylab="")
    text(.6, .5, i, cex=1.5, adj=.5)
    plot(confint(energy.mmc.4W[[i]]$mca$glht, calpha=energy.tpmc),
         xlim=c(-.4, 2.8), ylim=c(0,4), main="", xlab="")
    par(old.mar)
  },s={
    old.cex <- par(cex=.9)
    plot(x=0:1, y=0:1, type="n", bty="n",
         xaxt="n", yaxt="n", xlab="", ylab="")
    text(.7, .7, i, cex=1.2, adj=0)
    frame()
    par(new=TRUE)
    plot(energy.mmc.4W[[i]]$lmat, xlim=c(-.4, 2.8), main="",
         comparisons.per.page=4, xlabel.print=FALSE, cex=.8,
         col.signif='red')
    par(old.cex)
  }
       )
par(mfrow=c(1,1))
