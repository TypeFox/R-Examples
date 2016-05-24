## The cc176 data comes from Cochran and Cox, page 176
## Cochran, W. G. and Cox, G. M. (1957).
## Experimental Designs, Wiley.

if.R(r=
old.contr <- options(contrasts=c("contr.treatment", "contr.treatment"))
,s={})

data(cc176)

current.levels <- c("galvanic","faradic","60.cycle","25.cycle")
cc176 <-
  data.frame(wt.d=cc176[,1],
             wt.n=cc176[,2],
             n.treats=ordered(rep(c(1,3,6), length=96)),
             current=ordered(rep(rep(current.levels,
               rep(3,4)),8), levels=current.levels),
             minutes=ordered(rep(rep(c(1,2,3,5),rep(12,4)), 2)),
             rep=rep(c("I","II"), c(48,48)))
contrasts(cc176$n.treats)
contrasts(cc176$n.treats) <-
  if.R(s=
       contr.poly(c(1,3,6))
       ,r=
       contr.poly(3, scores=c(1,3,6))
       )
contrasts(cc176$n.treats)
if.R(r={
     contrasts(cc176$current) <- "contr.treatment"
     contrasts(cc176$n.treats) <- "contr.treatment"
   }
     ,s={})


if.R(r=
options(old.contr)
,s={})

cc176.aov <- aov(wt.d ~ rep + wt.n + n.treats + wt.n*current, data=cc176)
summary(cc176.aov)

summary(cc176.aov,
        split=list(n.treats=list(n.treats.lin=1, n.treats.quad=2)),
        expand.split=FALSE)



## almost MMC Figure 3, mmc plot of current, default settings for display
if.R(s={
  old.par <- par(mar=c(5,4,4,8)+.1)
  cc176.mmc <-
    multicomp.mmc(cc176.aov, focus="current",
                  ry=c(56,72), x.offset=1.5,
                  valid.check=FALSE, ## Tukey method (slightly narrower bounds)
                  plot=FALSE)
  print(cc176.mmc)
  plot(cc176.mmc)
  par(old.par)
}
,r={
  cc176.mmc <-
    mmc(cc176.aov, linfct=mcp(current="Tukey",
        `covariate_average`=TRUE, `interaction_average`=TRUE))
  print(cc176.mmc)
  plot(cc176.mmc)
})


## MMC Figure 3, mmc plot of current, control of display parameters
if.R(
     s={
       old.par <- par(mar=c(5,4,4,8)+.1)
       plot(cc176.mmc, iso.name=FALSE, col.iso=16,
            print.mca=TRUE, print.lmat=FALSE)
       par(old.par)
     },r={
       old.par <- par(mar=c(15,4,4,8)+.1)
       plot(cc176.mmc,  x.offset=1.5,
            iso.name=FALSE, col.iso=16,
            col.mca.signif="red", lty.mca.not.signif=2,
            print.mca=TRUE, print.lmat=FALSE)
       par(old.par)
     })



## MMC Figure 4, current with orthogonal contrasts and control of display
cc176.orth <- cbind(  "g-f"=c( 1,-1, 0, 0),
                    "60-25"=c( 0, 0, 1,-1),
                    "gf-AC"=c( 1, 1,-1,-1))
dimnames(cc176.orth)[[1]] <- levels(cc176$current)

if.R(s={

  cc176.mmc <-
    multicomp.mmc(cc176.aov, focus="current", method="tukey",
                  focus.lmat=cc176.orth, valid.check=FALSE,
                  plot=FALSE)
##                  ry=c(56,72), x.offset=1.5,

  print(cc176.mmc)
  old.par <- par(mar=c(5,4,4,8)+.1)
       plot(cc176.mmc, iso.name=FALSE, col.iso=16,
            col.mca.signif="red", lty.mca.not.signif=2,
            print.mca=FALSE, print.lmat=TRUE)
  par(old.par)
}
,r={
  cc176.mmc <-
    mmc(cc176.aov,
        linfct=mcp(current="Tukey",
          `covariate_average`=TRUE, `interaction_average`=TRUE),
        focus.lmat=cc176.orth,
        )
  print(cc176.mmc)


  old.par <- par(mar=c(15,4,4,8)+.1)
  plot(cc176.mmc, x.offset=1.5,
       iso.name=FALSE, col.iso=16,
       col.mca.signif="red", lty.mca.not.signif=2,
       col.lmat.signif="blue", lty.lmat.not.signif=2,
       print.mca=FALSE, print.lmat=TRUE)
  par(old.par)
  })



## MMC Figure 10, current, pairwise and orthogonal contrasts
if.R(
     s={
       old.par <- par(mar=c(5,4,4,8)+.1)
       plot(cc176.mmc, iso.name=FALSE, col.iso=16,
            print.mca=TRUE, print.lmat=TRUE,
            col.lmat.signif=6)
       par(old.par)
     },r={
       old.par <- par(mar=c(15,4,4,8)+.1)
       plot(cc176.mmc,  x.offset=1.5,
            iso.name=FALSE, col.iso=16,
            col.mca.signif="red", lty.mca.not.signif=2,
            col.lmat.signif="blue", lty.lmat.not.signif=2,
            print.mca=TRUE, print.lmat=TRUE)
       par(old.par)
     })
