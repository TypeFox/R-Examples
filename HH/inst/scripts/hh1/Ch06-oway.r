library(HH)

#### oway/code/catalystm.s
#### oway/code/catalystm.r
data(catalystm)
     
if.R(s={
tpg <- trellis.par.get("box.rectangle")
tpg.old <- tpg
tpg$lwd <- 1
trellis.par.set("box.rectangle", tpg)
},
     r={})

if.R(s=
t(bwplot(catalyst ~ concent, data=catalystm,
         scales=list(cex=1.5),
         xlab=list("concentration", cex=1.5),
         ylab=list("catalyst",cex=1.5)))
,r=
bwplot(concent ~ catalyst, data=catalystm,
         scales=list(cex=1.5),
         ylab=list("concentration", cex=1.5),
         xlab=list("catalyst",cex=1.5))
)
## export.eps(hh("oway/figure/catalystm1.eps"))

catalystm1.aov <- aov(concent ~ catalyst, data=catalystm)
anova(catalystm1.aov)

model.tables(catalystm1.aov, "means")

hov(concent ~ catalyst, data=catalystm)
print(position=c(0,.2,1,1),
      hovPlot(concent ~ catalyst, data=catalystm)
      )
## export.eps(hh("oway/figure/catalystm-hov.eps"))

if.R(s=trellis.par.set("box.rectangle", tpg.old),
     r={})


#### mcomp/code/catalystm-mmc3.s
#### mcomp/code/catalystm-mmc3.r
## Use glht.mmc with R.
## Use multicomp.mmc with S-Plus.

## data and ANOVA

catalystm.mca <-
if.R(r=glht(catalystm1.aov, linfct = mcp(catalyst = "Tukey")),
     s=multicomp(catalystm1.aov, plot=FALSE))
plot(catalystm.mca) ## S-Plus always and R for small number of levels
## export.eps(hh("oway/figure/catalystm2.eps"))
if.R(r=confint(catalystm.mca),
     s=catalystm.mca)

## pairwise comparisons
old.omd <- par(omd=c(0,.95,0,1))
catalystm.mmc <-
if.R(r=mmc(catalystm1.aov, linfct = mcp(catalyst = "Tukey")),
     s=multicomp.mmc(catalystm1.aov, plot=FALSE))
catalystm.mmc
if.R(s=plot(catalystm.mmc, x.offset=1),
     r=plot(catalystm.mmc, ry=c(50,58), x.offset=1.8))
## export.eps(hh("mcomp/figure/catalystm-mmc-mca.eps"))
## tiebreaker
plotMatchMMC(catalystm.mmc$mca)
plot(catalystm.mmc$none, xlab="concent")



## user-specified contrasts       A  B  C  D
catalystm.lmat <- cbind("AB-D" =c(1, 1, 0,-2),
                        "A-B"  =c(1,-1, 0, 0),
                        "ABD-C"=c(1, 1,-3, 1))
dimnames(catalystm.lmat)[[1]] <- levels(catalystm$catalyst)

catalystm.lmat

catalystm.mmc <-
  if.R(r=mmc(
         catalystm1.aov,
         linfct = mcp(catalyst = "Tukey"),
         focus.lmat=catalystm.lmat)
       ,s=multicomp.mmc(catalystm1.aov, lmat=rbind(0,catalystm.lmat),
                     plot=FALSE)
)

catalystm.mmc
plot(catalystm.mmc, ry=c(50,58), x.offset=1.2)
## tiebreaker
plotMatchMMC(catalystm.mmc$lmat, col.signif='blue')
par(old.omd)


#### oway/code/batch.s
data(batch)

if.R(s=
t(bwplot(Batch ~ Calcium, data=batch,
         scales=list(cex=1.5),
         xlab=list(cex=1.5),
         ylab=list("Batch", cex=1.5)))
,r=
bwplot(Calcium ~ Batch, data=batch,
         scales=list(cex=1.5),
         xlab=list(cex=1.5),
         ylab=list(cex=1.5))
)
## export.eps(hh("oway/figure/batch-data.eps"))

## anova
batch1.aov <- aov(Calcium ~ Batch, data=batch)
anova(batch1.aov)

## means
model.tables(batch1.aov, "means")

## homogeneity of variance
hov(Calcium ~ Batch, data=batch)
print(position=c(0,.2,1,1),
      hovPlot(Calcium ~ Batch, data=batch)
      )
## export.eps(hh("oway/figure/batch-hov.eps"))


#### oway/code/turkey-oway.s
data(turkey)
turkey <- turkey ## local copy

if.R(s=
     t(bwplot(diet ~ wt.gain, data=turkey,
              scales=list(cex=1.5),
              xlab=list("Weight Gain", cex=1.5),
              ylab=list("diet",cex=1.5)))
     ,r=
     bwplot(wt.gain ~ diet, data=turkey,
            scales=list(cex=1.5),
            ylab=list("Weight Gain", cex=1.5),
            xlab=list("diet",cex=1.5))
     )
## export.eps(hh("oway/figure/turkey.f1.eps"))

turkey.aov <- aov(wt.gain ~ diet, data=turkey)
summary(turkey.aov)
model.tables(turkey.aov, type="means", se=TRUE)
tapply(turkey$wt.gain, turkey$diet, mean)

contrasts(turkey$diet)

contrasts(turkey$diet) <-
  cbind(control.vs.treatment=c(1,-.25,-.25,-.25,-.25),
        A.vs.B              =c(0, .5,  .5, -.5, -.5 ),
        amount              =c(0, .5, -.5,  .5, -.5 ),
        A.vs.B.by.amount    =c(0, .5, -.5, -.5,  .5 ))
contrasts(turkey$diet)

tapply(turkey$wt.gain, turkey$diet, mean) %*% contrasts(turkey$diet)

turkey2.aov <- aov(wt.gain ~ diet, data=turkey)
summary(turkey2.aov)
summary(turkey2.aov,
        split=list(diet=list(
                     control.vs.treatment=1,
                     A.vs.B=2,
                     amount=3,
                     A.vs.B.by.amount=4)))
if.R(s=
multicomp(turkey2.aov,
          focus="diet",
          lmat=rbind(0,contrasts(turkey$diet)),
          method="scheffe",
          bounds="both",
          plot=TRUE,
          comparisons="none")
## export.eps(hh("oway/figure/turkey.scheffe.eps"))
     ,r={

       ## Means
       turkey.Means.multcomp <- glht(turkey2.aov,
                               linfct=mcp(diet="Means"))

       ## construct linfct matrix to correspond to constrasts
       linfct.contrasts <- t(cbind(control=0,
                                   turkey.Means.multcomp$linfct[,2:5]) %*%
                             contrasts(turkey$diet))
       linfct.contrasts[1,] <- linfct.contrasts[1,]*5
       turkey.multcomp <- glht(turkey2.aov,
                               linfct=mcp(diet=linfct.contrasts))
       turkey.multcomp$linfct

       scheffe.quantile <- sqrt(4*qf(.95, 4, 25))
         
       ctm <- confint(turkey.multcomp, calpha=scheffe.quantile)
       print(ctm)
       old.oma <- par(oma=c(0,4.5,0,0))
       plot(ctm)  ## plot the result of the confint in order to get the calpha
       par(old.oma)
     }
     )



#### splus.library/hov-bf.s
#### splus.library/hov.plot.s
#### The functions here are now included in the HH package.


#### oway/code/patient.r
data(patient)
