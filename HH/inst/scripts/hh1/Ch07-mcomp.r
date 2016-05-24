library(HH)

#### mcomp/code/weightloss.r
#### mcomp/code/weightloss.s
data(weightloss)

if.R(r=
     bwplot(loss ~ group, data=weightloss,
            scales=list(cex=1.5),
            ylab=list("Weight Loss", cex=1.5),
            xlab=list("group",cex=1.5))
     ,s=
     t(bwplot(group ~ loss, data=weightloss,
              scales=list(cex=1.5),
              xlab=list("Weight Loss", cex=1.5),
              ylab=list("group",cex=1.5)))
     )
## export.eps(hh("mcomp/figure/weightloss-data.eps"))

weightloss.aov <- aov(loss ~ group, data=weightloss)
summary(weightloss.aov)

if.R(
     s={
weightloss.dunnett <-
multicomp(weightloss.aov,
          method="dunnett", comparisons="mcc",
          bounds="lower", control=4, plot=TRUE)
     },
     r={
weightloss.dunnett <- 
glht(weightloss.aov,
     linfct = mcp(group=contrMat(rep(10,5), base=4)),
     alternative = 'greater')
weightloss.dunnett
plot(weightloss.dunnett)
     }
     )
## export.eps(hh("mcomp/figure/weightloss.dunnet.eps"))

if.R(r={
  weightloss.mmc <-
    mmc(weightloss.aov,
        linfct = mcp(group=contrMat(rep(10,5), base=4)),
        alternative = 'greater')
  weightloss.mmc
  plot(weightloss.mmc)
}, s={
  weightloss.mmc <-
    multicomp.mmc(weightloss.aov,
                  method="dunnett", comparisons="mcc",
                  bounds="lower", control=4, plot=TRUE)
}
     )
## export.eps(hh("mcomp/figure/weightloss.dunnet.mmc.eps"))
plotMatchMMC(weightloss.mmc$mca)


#### oway/code/turkey-oway.s
## Ch06-oway.r includes this code
#### oway/code/turkey-oway.s
data(turkey)
turkey.aov <- aov(wt.gain ~ diet, data=turkey)
scheffe.quantile <- sqrt(4*qf(.95, 4, 25))
contrasts(turkey$diet) <-
  cbind(control.vs.treatment=c(1,-.25,-.25,-.25,-.25),
        A.vs.B              =c(0, .5,  .5, -.5, -.5 ),
        amount              =c(0, .5, -.5,  .5, -.5 ),
        A.vs.B.by.amount    =c(0, .5, -.5, -.5,  .5 ))

#### mcomp/code/turkey-mmc.s
## follows oway/code/turkey-oway.s
## follows the turkey section of Ch06-oway.r

old.par <-      ## resize the graphics device before plotting.
  if.R(r=par(omd=c(0,.82,0,1)),
       s=par(oma=c(0,1,0,3.5)))

if.R(r={
  turkey.mmc <- mmc(turkey.aov, calpha=scheffe.quantile)  ## original treatment contrasts
  plot(turkey.mmc, ry=c(3.5,10))
  turkey.mmc
}, s={
  turkey.mmc <- multicomp.mmc(turkey2.aov, ## constructed contrasts
                              method="scheffe", x.offset=.9)
})
## export.eps(hh("mcomp/figure/turkey-mca.eps"))

## We need to break the visual ties from the overprinting of this
## graph.  Standard R multiple comparisons plot with contrasts
## ordered to match the mmc plot ordering.
plotMatchMMC(turkey.mmc$mca)
## export.eps(hh("mcomp/figure/turkey-mca2.eps"))

## Look at the generated lmat.  We reordered the columns to match the
## order of the mmc plot and changed the signs to be larger minus
## smaller.
zapsmall(turkey.mmc$mca$lmat)


## construct the lmat for the contrasts
if.R(r={
  contrasts(turkey$diet) ## these are the constructed contrasts
  dimnames(contrasts(turkey$diet))[[2]][1] <- "cvt"
  contrasts(turkey$diet) ## these are the constructed contrasts
  turkey.mmc <- mmc(turkey.aov, calpha=scheffe.quantile, focus="diet",
                    focus.lmat=contrasts(turkey$diet))
  plot(turkey.mmc, ry=c(3.5,10))
  turkey.mmc
}, s={
contrasts(turkey$diet)
turkey.lmat <- rbind("(Intercept)"=0, contrasts(turkey$diet))
turkey.lmat

turkey.mmc <- multicomp.mmc(turkey2.aov, method="scheffe", x.offset=.9,
                     lmat=turkey.lmat)
})
## export.eps(hh("mcomp/figure/turkey-lmat.eps"))

## We need to break the visual ties from the overprinting of this lmat
## graph.
plotMatchMMC(turkey.mmc$lmat, col.signif='blue', ylabel.inside=TRUE)
## export.eps(hh("mcomp/figure/turkey-mcp2.eps"))

par(old.par)


#### mcomp/code/inconsistent.s
## Demonstration of inconsistent confidence intervals when sample
## sizes are unequal.

group <- LETTERS[1:4]
n <- c(5,100,100,5)
ybar <- c(2, 2.1, 2.8, 3)
s=.8

if.R(r={
  group <- factor(group)
  inconsistent.aov <- aovSufficient(ybar ~ group, weights=n, sd=s)
  anova(inconsistent.aov)
  inconsistent.glht <- glht(inconsistent.aov,
                            linfct=mcp(group="Tukey"),
                            vcov.=vcovSufficient,
                            df=inconsistent.aov$df.residual)
  crit.point <- qtukey(.95, 4, 206)/sqrt(2)
  confint(inconsistent.glht, calpha=crit.point)
  plot(inconsistent.glht,
       xlab="simultaneous confidence limits, Tukey method",
       sub=paste("critical point =", round(crit.point, 4)))
}, s={
  mca <- multicomp.mean(group, n, ybar, s,
                        ylabel="response", focus="group", plot=FALSE)
  mca <- multicomp.reverse(mca)
  plot(mca)
  title(main=paste("critical point =", round(mca$crit.point,4)))
  mca
  ## export.eps(hh("mcomp/figure/unequal1.eps"))
}
)
     
if.R(r={
  inconsistent.mmc <- mmc(inconsistent.aov,
                          linfct=mcp(group="Tukey"),
                          vcov.=vcovSufficient,
                          df=inconsistent.aov$df.residual,
                          calpha=qtukey(.95, 4, 206)/sqrt(2))
  print(inconsistent.mmc)
  plot(inconsistent.mmc$none)
  plot(inconsistent.mmc, ry=c(1.65, 3.25), x.offset=.5)
  plotMatchMMC(inconsistent.mmc$mca)
}, s={
  none <-
    multicomp.mean(group, n, ybar, s,
                   ylabel="response", focus="group", plot=TRUE,
                   comparisons="none",
                   method="tukey", valid.check=FALSE)
  title(main=paste("critical point =", round(mca$crit.point,4)))
  ## export.eps(hh("mcomp/figure/unequal2.eps"))
  
  mca$height <-
    matrix(outer(none$table[,"estimate"], none$table[,"estimate"], "+"
                 )[outer(1:4, 1:4, ">")],
           1, 6, dimnames=list(NULL, dimnames(mca$table)[[1]]))
  

  tmp <- list(mca=mca, none=none)
  oldClass(tmp) <- c("mmc.multicomp", "list" )
  
  plot.mmc.multicomp(tmp, x.offset=.4, ry=c(1.9,3))
  ## export.eps(hh("mcomp/figure/unequal3.eps"))
  
  tmp$mca <- multicomp.order(tmp$mca)
  tmp$none <- multicomp.order(tmp$none)
  plot(tmp$mca)
  title(main="Same information as plot unequal1,\n with rows ordered by height on mmc plot ")
  ## export.eps(hh("mcomp/figure/unequal4.eps"))
})


#### splus.library/mmc.multicomp.s
## These functions are included in the HH package

#### mcomp/code/catalystm-mmc3.r
#### mcomp/code/catalystm-mmc3.s
#### mcomp/code/catalystm.s
## The entire set of catalystm examples are now in file Ch06-oway.r



#### mcomp/code/mmc.explain.r
#### mcomp/code/mmc.explain.s
tpg.col <- function(col.in=1) {
  if.R(s=
       tpg <- attr(get(".Device", where=0), "trellis.settings")
       ,r=
       tpg <- trellis.par.get()
       )
  old.tpg <- tpg
  
  tpg$axis.line$col <- col.in
  trellis.par.set("axis.line", tpg$axis.line)
  
  tpg$add.text$col <- col.in
  trellis.par.set("add.text", tpg$add.text)
  
  tpg$plot.line$col <- col.in
  trellis.par.set("plot.line", tpg$plot.line)
  
  tpg$plot.symbol$col <- col.in
  trellis.par.set("plot.symbol", tpg$plot.symbol)
}




mmc.explain <- function(group, n, ybar, ms.5, crit.point=1.96, ylabel="ylabel",
                        method="user-defined",
                        xlim,
                        ylim=xlim,
                        exit=FALSE,
                        col.in=c(1,1,1),
                        browser.in=FALSE) {
  tpg.col(col.in[1])
  par(col=col.in[1])
  ybar.mat <- matrix(ybar, length(ybar), length(ybar))
  ry <- range(ybar)
  if (missing(xlim)) xlim <- ry + c(-1,1)*.7*diff(ry)
  xyplot(as.vector(ybar.mat) ~ as.vector(t(ybar.mat)),
         aspect=1, xlim=xlim, ylim=ylim,
         xlab=list("h=ybar", col=col.in[1]),
         ylab=list("v=ybar", col=col.in[1]),
         main=list(paste("Multiple comparisons of response variable:", ylabel),
           col=col.in[1]),
         par.settings=list( ## R lattice needs this, S-Plus ignores it
           clip=list(panel="off"),
           layout.widths=list(left.padding=10),
           layout.heights=list(bottom.padding=10)),
         crit.point=crit.point,
         uy=ybar,
         uy.labels=group,
         ms.5=ms.5,
         exit=exit,
         col.in=col.in,
         browser.in=browser.in,
         scales=list(col=col.in[1]),
         panel=function(x, y, ...,
           crit.point, uy, uy.labels, ms.5, browser.in, exit, col.in) {

           segments.me <- if.R(r=lsegments,
                               s=segments)
           text.me <- if.R(r=ltext,
                           s=text)
           mtext.me <- if.R(r=ltext,
                            s=mtext)


           ## (ybar, ybar) points
           tpg.col(col.in[1])
           panel.xyplot(x, y, ...,
                        col=col.in[1])
           if (browser.in) {
             mtext.me(get("xlab", frame=sys.parent())$label,
                   side=1, line=2.9, col=col.in[1])
             mtext.me(get("ylab", frame=sys.parent())$label,
                   side=2, line=2.9, col=col.in[1])
           }
           ## square in (h,v) and (d,m) coordinates
            for (i in uy) {
              segments.me(min(uy), i, max(uy), i, lty=2, col=col.in[1])
              segments.me(i, min(uy), i, max(uy), lty=2,
                       col=col.in[1])
            }
           ## labels for constant v and h lines
           text.me(x=uy, y=rep(min(y)-.5, length(uy)), uy.labels,
                col=col.in[1])
           text.me(y=uy, x=rep(min(y)-.5, length(uy)), uy.labels,
                col=col.in[1])

           ## means on m axis
           panel.abline(a=0, b=1, col=col.in[1])
           text.me(x=max(y)+c(1.1,.6), y=max(y)+c(.5,1.5),
                c("d = h-v = 0","(m-axis)"),
                srt=45, adj=0,
                col=col.in[1])
           ## ticks on m axis
           segments.me(x1=uy-.025*diff(range(y)),
                    y1=uy+.025*diff(range(y)),
                    x2=uy+.025*diff(range(y)),
                    y2=uy-.025*diff(range(y)),
                    col=col.in[1])           

           ## (0,0) point and axes
           segments.me(x1=44.5,
                    y1=44.25,
                    x2=44.5,
                    y2=47,
                    xpd=TRUE,
                    col=col.in[1])
           segments.me(x1=44.25,
                    y1=44.5,
                    x2=47,
                    y2=44.5,
                    xpd=TRUE,
                    col=col.in[1])
           text.me(x=c(44,44.5),y=c(44.5,44), c("0","0"),
                    col=col.in[1])
           segments.me(x1=44.5,
                    y1=44.5,
                    x2=45.5,
                    y2=45.5,
                    xpd=TRUE,
                    col=col.in[1])
           segments.me(x1=45.5,
                    y1=45.5,
                    x2=par("usr")[1],
                    y2=par("usr")[1],
                    xpd=TRUE,
                    lty=2,
                    col=col.in[1])
           segments.me(x1=40.5,
                    y1=48.5,
                    x2=46.5,
                    y2=42.5,
                    xpd=TRUE,
                    col=col.in[1])
           text.me(x=c(42.1,41.9), y=c(47.5,46.5),
                c("(d-axis)","m = (h+v)/2 = 0"), srt=-45,
                col=col.in[1])
           if (browser.in) browser()
           if (exit==1) return()
           ## export.eps(hh("mcomp/figure/mmc1-a.eps"))


           tpg.col(col.in[2])
           ## index for means on m axis
           text.me(x=uy-.3,
                y=uy+.3,
                uy.labels, col=col.in[2])

           ## ticks in m coordinates
           panel.abline(a=8, b=1,
                        col=col.in[2])
           segments.me(x1=uy-8/2,
                    y1=uy+8/2,
                    x2=uy-8/2+.05*diff(range(y)),
                    y2=uy+8/2-.05*diff(range(y)),
                    col=col.in[2])           
           ## m labels, ybar
           text.me(x=uy-8/2-.4,
                y=uy+8/2+.4,
                round(uy,1),
                col=col.in[2])
           text.me(x=uy[1]-8/2-.4+.8,
                y=uy[1]+8/2+.4+.8,
                "m",
                col=col.in[2])
           text.me(x=max(y)-9.3, y=max(y)+1,
                "d = h-v = constant",
                srt=45, adj=0,
                col=col.in[2])
           ## perpendiculars from A to m-axis
           segments.me(uy[1],uy[1], uy[1]-8/2, uy[1]+8/2, lty=3,
                col=col.in[2]) ## to m-axis

           if (browser.in) browser()
           if (exit==4) return()
           ## export.eps(hh("mcomp/figure/mmc1-b0.eps"))

           
           ## differences on d axis
           panel.abline(a=2*min(y)-.2*diff(range(y))-.8, b=-1,
                        col=col.in[2])
           text.me(x=min(y)-2.4, y=min(y)-1.5, srt=-45,
                "m = (h+v)/2 = constant",
                col=col.in[2])
           segments.me(x1=uy-(uy-min(uy))/2-.1*diff(range(y))-.4,
                    y1=min(uy)-(uy-min(uy))/2-.1*diff(range(y))-.4,
                    x2=uy-(uy-min(uy))/2-.05*diff(range(y))-.4,
                    y2=min(uy)-(uy-min(uy))/2-.05*diff(range(y))-.4,
                    col=col.in[2])
           text.me(x=uy-(uy-min(uy))/2+.0*diff(range(y))-.6,
                y=min(uy)-(uy-min(uy))/2+.0*diff(range(y))-.4,
                c(paste(uy.labels[1], uy[1]-min(uy), sep=": "),
                  uy.labels[2:4]),
                adj=0,
                col=col.in[2])
           segments.me(x1= (0:6)/2+min(uy)-.10*diff(range(y))-.4,
                    y1=-(0:6)/2+min(uy)-.10*diff(range(y))-.4,
                    x2= (0:6)/2+min(uy)-.15*diff(range(y))-.4,
                    y2=-(0:6)/2+min(uy)-.15*diff(range(y))-.4,
                    col=col.in[2])
           text.me(x= (0:7)/2+min(uy)-.18*diff(range(y))-.4,
                y=-(0:7)/2+min(uy)-.18*diff(range(y))-.4,
                c(0:6, "d"),
                col=col.in[2])
           ## perpendicular from A to d-axis
           segments.me(uy[1],uy[4], uy[1]-3.8, uy[4]-3.8, lty=3,
                col=col.in[2]) ## to d-axis
           if (exit==2) return()
           ## export.eps(hh("mcomp/figure/mmc1-b.eps"))

           tpg.col(col.in[3])
           ## CI for ybar2-ybar4
           segments.me(x1=uy[2]-ms.5*sqrt(1/n[2]+1/n[4])*crit.point/2,
                    y1=uy[4]+ms.5*sqrt(1/n[2]+1/n[4])*crit.point/2,
                    x2=uy[2]+ms.5*sqrt(1/n[2]+1/n[4])*crit.point/2,
                    y2=uy[4]-ms.5*sqrt(1/n[2]+1/n[4])*crit.point/2,
                col=col.in[3])
         }
)
}



## example program using catalystm data for constructing MMC plot
## based on the HH book, Section 7.2.2

## Repeated from Ch06-oway.r
data(catalystm)
catalystm1.aov <- aov(concent ~ catalyst, data=catalystm)

## pairwise comparisons
catalystm.mmc <-
if.R(r=mmc(catalystm1.aov, linfct = mcp(catalyst = "Tukey")),
     s=multicomp.mmc(catalystm1.aov, plot=FALSE))
catalystm.mmc
old.omd <- par(omd=c(0,.95,0,1))  ## look ahead to determine these values
if.R(r=plot(catalystm.mmc, ry=c(50,58), x.offset=1.8),
     s=plot(catalystm.mmc, x.offset=1))
## export.eps(hh("mcomp/figure/catalystm-mmc-mca.eps"))

## new for Ch07-mcomp.r
## just the B-D contrast
`lmat.B-D` <- catalystm.mmc$mca$lmat[,"B-D", drop=FALSE]
if.R(r=dimnames(`lmat.B-D`)[[1]] <- levels(catalystm$catalyst),
     s=dimnames(`lmat.B-D`)[[1]] <- c("(I)", levels(catalystm$catalyst)))

if.R(r={catalystm.mmc <-
          mmc(catalystm1.aov, focus.lmat=`lmat.B-D`)
        print(catalystm.mmc)
        plot(catalystm.mmc, ry=c(50,58), x.offset=1.8)
      }, s={
        catalystm.mmc <- multicomp.mmc(catalystm1.aov,
                                       lmat=`lmat.B-D`, plot=FALSE)
        plot(catalystm.mmc, x.offset=1)
      }
     )
par(old.omd)
## export.eps(hh("mcomp/figure/mmc2.eps"))

group <- levels(catalystm$catalyst)
n <- c(4,4,4,4)
ybar <- tapply(catalystm$concent, catalystm$catalyst, mean)

if.R(s={
  ms.5 <- summary(catalystm1.aov)$"Mean Sq"[2]^.5
  crit.point <- catalystm.mmc$none$crit.point
},r={
  ms.5 <- summary.aov(catalystm1.aov)[[1]][2,"Mean Sq"]^.5
  crit.point <- catalystm.mmc$mca$crit.point
})

##mmc1 <-
xlim.explain <- if.R(r=c(45.5,62.5),
                     s=c(46.5,61.0))
mmc.explain(group, n, ybar, ms.5, crit.point,
            ylabel="concent", method="tukey",
            xlim=xlim.explain)

##if.R(r={mmc1 <- update(mmc1, par.settings=list(clip=list(panel="off")))
##       print(mmc1, position=c(.1,.1,1,1))},
##     s=print(mmc1))
## export.eps(hh("mcomp/figure/mmc1.eps"))


## #with this version of the call, the routine will stop
## #twice, at intermediate stages.  We manually snapshot these stages.
## mmc.explain(group, n, ybar, ms.5, crit.point,
##             ylabel="concent", method="tukey",
##             browser.in=TRUE)



## With these four calls to mmc.explain and one call to par,
## each of the stopping points is colored differently.

##mmc1a <-
mmc.explain(group, n, ybar, ms.5, crit.point,
            ylabel="concent", method="tukey",
            xlim=xlim.explain,
            exit=1, col.in=c(1,1,1))
##if.R(r={mmc1a <- update(mmc1a, par.settings=list(clip=list(panel="off")))
##        print(mmc1a, position=c(.1,.1,1,1))},
##     s=print(mmc1a))
## export.eps(hh("mcomp/figure/mmc1-a.eps"))

##mmc1b0 <-
mmc.explain(group, n, ybar, ms.5, crit.point,
            ylabel="concent", method="tukey",
            xlim=xlim.explain,
            exit=4, col.in=c(65,1,1))
##if.R(r={mmc1b0 <- update(mmc1b0, par.settings=list(clip=list(panel="off")))
##        print(mmc1b0, position=c(.1,.1,1,1))},
##     s=print(mmc1b0))
## export.eps(hh("mcomp/figure/mmc1-b0.eps"))

##mmc1b <-
mmc.explain(group, n, ybar, ms.5, crit.point,
            ylabel="concent", method="tukey",
            xlim=xlim.explain,
            exit=2, col.in=c(65,1,1))
##if.R(r={mmc1b <- update(mmc1b, par.settings=list(clip=list(panel="off")))
##        print(mmc1b, position=c(.1,.1,1,1))},
##     s=print(mmc1b))
## export.eps(hh("mcomp/figure/mmc1-b.eps"))

##mmc1 <-
mmc.explain(group, n, ybar, ms.5, crit.point,
            ylabel="concent", method="tukey",
            xlim=xlim.explain,
            exit=3, col.in=c(65,65,1))
##if.R(r={mmc1 <- update(mmc1, par.settings=list(clip=list(panel="off")))
##        print(mmc1, position=c(.1,.1,1,1))},
##     s=print(mmc1))
## export.eps(hh("mcomp/figure/mmc1.eps"))

par(col=1) ## restore col value set inside mmc.explain
## if (col.in[3]==1) then all will be normal when you return from this call


#### mcomp/code/pulmonary.s
#### mcomp/code/pulmonary-book.s
#### mcomp/code/pulmonary2.s

## please see file HH-R.package/HH/demo/MMC.pulmonary.R

#### HH-R.package/HH/man/aov.sufficient.Rd
     
