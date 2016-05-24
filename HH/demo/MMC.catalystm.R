## The catalystm data comes from Montgomery (1997)
## Montgomery, D. C. (1997).
## Design and Analysis of Experiments, Wiley, 4th edition.



data(catalystm)
catalystm1.aov <- aov(concent ~ catalyst, data=catalystm)

catalystm.mmc <-
if.R(r=
     mmc(catalystm1.aov, linfct = mcp(catalyst = "Tukey"))
    ,s=
     multicomp.mmc(catalystm1.aov, plot=FALSE)
)

old.mar <- if.R(s=par(mar=c(5,8,4,4)+.1),
                r=par(mar=c(12,4,4,3)+.1))
plot(catalystm.mmc, x.offset=1.6, ry.mmc=c(50.5,57),
     print.lmat=FALSE)

catalystm.lmat <- cbind("AB-D" =c( 1, 1, 0,-2),
                        "A-B"  =c( 1,-1, 0, 0),
                        "ABD-C"=c( 1, 1,-3, 1))
dimnames(catalystm.lmat)[[1]] <- levels(catalystm$catalyst)

catalystm.mmc <-
if.R(r=
     mmc(catalystm1.aov, linfct = mcp(catalyst = "Tukey"),
              focus.lmat=catalystm.lmat)
     ,s=
     multicomp.mmc(catalystm1.aov, focus.lmat=catalystm.lmat,
                   plot=FALSE)
)

plot(catalystm.mmc, x.offset=1.6, ry.mmc=c(50.5,57))



lty.contr0 <- if.R(r=3, s=2)
lty.iso    <- if.R(r=2, s=8)

## MMC Figure 1, pairwise contrasts
plot(catalystm.mmc, x.offset=1.6, ry.mmc=c(50.5,57),
     print.lmat=FALSE,
     col.mca.signif='red', lty.mca.not.signif=4,
     lty.contr0=lty.contr0, col.contr0='darkgray',
     lty.iso=lty.iso, col.iso='darkgray')

## MMC Figure 5, user-specified contrasts
plot(catalystm.mmc, x.offset=1.6, ry.mmc=c(50.5,57),
     col.mca.signif='red', lty.mca.not.signif=4,
     lty.contr0=lty.contr0, col.contr0='darkgray',
     lty.iso=lty.iso, col.iso='darkgray',
     col.lmat.signif='blue')

## both pairwise contrasts and user-specified contrasts
plot(catalystm.mmc, x.offset=1.6, ry.mmc=c(50.5,57),
     print.mca=TRUE,
     col.mca.signif='red', lty.mca.not.signif=4,
     lty.contr0=lty.contr0, col.contr0='darkgray',
     lty.iso=lty.iso, col.iso='darkgray',
     col.lmat.signif='blue')

par(old.mar)

## See the help file for plot.mmc.multicomp for more examples



## illustrate the construction of the isomeans grid and the contrasts
##
## This illustration is no longer in the demo.
## It is now in the file "HH/scripts/Ch07-mccomp.r".
