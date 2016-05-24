## c:/HOME/rmh/hh/dsgn/:
## cc176.s
## cc176-ancova-plot.s
## tires.r
## tires.s
## filmcoat.s
## filmcoat2.s
## gunload.s
## turkey.factors.s
## turkey.aov2.s
## turkey.means.s
## turkey.f2.r
## turkey.f2.s
## turkey.aov3.s
## turkey.aov3-match.s
## contrasts.s
## vulcan.s
## tires.latin.s
## tires.latin.2008.s
## -rwx------+   1 rmh None    80272 2004-06-08 16:18 dsgn.tex

library(HH)

## cc176.s
data(cc176)
cc176 <- cc176

## y=wt.d analysis ignoring x=wt.n
## (get the same y^2 ANOVA table as Cochran and Cox)
cc176.aov <- aov(wt.d ~ rep +        n.treats*minutes*current, data=cc176)
summary(cc176.aov)

## y=wt.d with x=wt.n as covariate
## (get essentially the same ANOVA as the approximate (y-bx)^2 ANOVA table
## in Cochran and Cox)
cc176.aov <- aov(wt.d ~ rep + wt.n + n.treats*minutes*current, data=cc176)
summary(cc176.aov)

summary(cc176.aov,
        split=list(n.treats=list(n.treats.lin=1, n.treats.quad=2)),
        expand.split=FALSE)
## dsgn/transcript/cc176-1.st

## repeat with just main effects of n.treats and current
cc176.aov <- aov(wt.d ~ rep + wt.n + n.treats + current, data=cc176)
summary(cc176.aov)

summary(cc176.aov,
        split=list(n.treats=list(n.treats.lin=1, n.treats.quad=2)),
        expand.split=FALSE)

## Full model including interaction of covariate with treatments.
## Four-way interaction used as the error term
## anova table is in file dsgn/transcript/cc176-full.st
cc176.aov <- aov(wt.d ~ rep + wt.n * n.treats*minutes*current
                            - wt.n : n.treats:minutes:current,
                 data=cc176)
summary(cc176.aov)
summary(cc176.aov,
        split=list(n.treats=list(n.treats.lin=1, n.treats.quad=2)),
        expand.split=FALSE)

## repeat with main effects of n.treats and current and wt.n interaction
cc176.aov <- aov(wt.d ~ rep + wt.n * (n.treats+current), data=cc176)
summary(cc176.aov)

## repeat with main effects of n.treats and current and wt.n interaction
## anova table is in file dsgn/transcript/cc176-model.st
cc176.aov <- aov(wt.d ~ rep + wt.n + n.treats + wt.n*current, data=cc176)
summary(cc176.aov)

summary(cc176.aov,
        split=list(n.treats=list(n.treats.lin=1, n.treats.quad=2)),
        expand.split=FALSE)


## mmc plot of current
if.R(s={
old.par <- par(mar=c(5,4,4,8)+.1)
multicomp.mmc(cc176.aov, focus="current") ## Sidak method
## export.eps(hh("dsgn/figure/cc176-6.eps"))

cc176.mmc <- multicomp.mmc(cc176.aov, focus="current",
                           valid.check=FALSE, ## Tukey method (slightly narrower bounds)
                           x.offset=1)
par(old.par)
print(cc176.mmc)
}
,r={
  cc176t <- cc176
  for (i in names(cc176t))
    if (is.factor(cc176t[[i]]))
      contrasts(cc176t[[i]]) <- contr.treatment(length(levels(cc176t[[i]])))
  sapply(cc176t, class)
  cc176t.aov <- aov(wt.d ~ rep + wt.n + n.treats + wt.n*current, data=cc176t)
  print(summary(cc176t.aov))
  cc176.mmc <- mmc(cc176t.aov, focus="current")
  old.omd <- par(omd=c(0,.85,0,1))
  plot(cc176.mmc, ry=c(55,72), x.offset=2.5)
  par(old.omd)
  print(cc176.mmc)
}
## export.eps(hh("dsgn/figure/cc176-7.eps"))
)

##                               g  f 60 25
current.lmat <- cbind("cc-gf"=c(-1,-1, 1, 1),
                      "25-60"=c( 0, 0,-1, 1),
                        "g-f"=c( 1,-1, 0, 0))
dimnames(current.lmat)[[1]] <- levels(cc176$current)
current.lmat

if.R(s={
old.par <- par(mar=c(5,4,4,8)+.1)
cc176.mmc <- multicomp.mmc(cc176.aov, focus="current",
                           lmat=cc176.mmc$none$lmat[,levels(cc176$current)] %*%
                           current.lmat,
                           valid.check=FALSE, ## Tukey method
                           x.offset=1)
par(old.par)
print(cc176.mmc)
}
,r={
  cc176.mmc <- mmc(cc176t.aov, focus="current", focus.lmat=current.lmat)
  old.omd <- par(omd=c(0,.85,0,1))
  plot(cc176.mmc, ry=c(55,72), x.offset=2.5)
  par(old.omd)
  print(cc176.mmc)
})



tpgss.old <- trellis.par.get("superpose.symbol")
tpgss <- Rows(trellis.par.get("superpose.symbol"),1:4)
tpgss$pch <- levels(cc176$minutes)
trellis.par.set("superpose.symbol", tpgss)

tmp <-
xyplot(wt.d ~ wt.n | n.treats*current, data=cc176,
       group=minutes,
       panel=function(x,y,...) {
         panel.superpose(x,y,...)
         panel.abline(lm(cc176$wt.d ~ cc176$wt.n), lty=3)
       },
       cex=1, par.strip.text=list(cex=1.4),
       strip=function(...) strip.default(strip.names=c(TRUE,TRUE), ...),
       scales=list(cex=1, alternating=FALSE),
       xlab=list(cex=1.2), ylab=list(cex=1.2),
       between=list(x=1, y=1)
       ##,
       ##key=list(
       ##  space="right", text=list(levels(cc176$minutes), adj=1), points=tpgss,
       ## border=1, title="minutes", cex.title=1.2,
       ##  cex=1)
       )
tmp
if (FALSE)
  print(tmp, position=c(0,-.3,1,1.1))  ## oversize on S-Plus screen, ok in file.
## export.eps(hh("dsgn/figure/cc176-1.eps"))


xyplot(wt.d ~ wt.n | n.treats * rep * current, data=cc176,
       group=minutes, panel=panel.superpose,
       layout=c(6,4),
       cex=1, par.strip.text=list(cex=.8),
       strip=function(...) strip.default(strip.names=c(TRUE,TRUE), ...),
       scales=list(cex=1, alternating=FALSE),
       xlab=list(cex=1.2), ylab=list(cex=1.2),
       between=list(x=c(1,1,3,1,1), y=1),
       key=list(
         space="right", text=list(levels(cc176$minutes), adj=1), points=tpgss,
         border=1, title="minutes", cex.title=1.2,
         cex=1))
## export.eps(hh("dsgn/figure/cc176-1-rep.eps"))


## adjust y for x
cc176$y.adj <- cc176$wt.d  -
  (cc176$wt.n - mean(cc176$wt.n))*coef(cc176.aov)["wt.n"]
## duplicate CC Table 5.17
cc176.means <- tapply(cc176$y.adj, cc176[,c("current","n.treats")], mean)
cc176.means
apply(cc176.means, 1, mean)
## dsgn/transcript/cc176-1.st


position(cc176$current) <- c(1.1,2.7,4.3,5.9)
position(cc176$n.treats) <- c(1,3,6)
tmp2 <-
interaction2wt(y.adj ~ current + n.treats, data=cc176,
               main.cex=1.6,
               scales=list(x=list(cex=.7), y=list(cex=.9, alternating=FALSE)))
if.R(r=tmp2,
     s=print(tmp2, position=c(0,.05,1,1)))
## export.eps(hh("dsgn/figure/cc176-2.eps"))

tmp3 <-
interaction2wt(y.adj ~ current + minutes + n.treats, data=cc176,
               main.cex=1.6,
               scales=list(x=list(cex=.6), y=list(cex=.9, alternating=FALSE)))
if.R(r=tmp3,
     s=print(tmp3, position=c(0,.05,1,1)))
## export.eps(hh("dsgn/figure/cc176-3.eps"))



## What does a three-way interaction look like?
## Here is one of the (3! = 6) possible plots
## If there were a significant 3-way interaction, the patterns
## in boxplots in the same row or column would not be the same.
## For example, we note that (y.adj ~ minutes) has a downslope in
## the galvanic by 3 panel and an uphill slope in the faradic by 3 panel.
## The ANOVA table tells us this difference in slope is not significant.
##
if.R(s=
t(bwplot(minutes ~ y.adj | n.treats + current, data=cc176,
         strip=function(...) strip.default(strip.names=c(TRUE,TRUE), ...),
         ylab=list(labels="minutes"),
         par.strip.text=list(cex=1.4),
         scales=list(
           x=list(cex=1, alternating=FALSE),
           y=list(cex=.8, alternating=FALSE)),
         between=list(x=1, y=1)))
,r=
bwplot(y.adj ~ minutes | n.treats + current, data=cc176,
         strip=function(...) strip.default(strip.names=c(TRUE,TRUE), ...),
         xlab=list(labels="minutes"),
         par.strip.text=list(cex=1.4),
         scales=list(
           y=list(cex=1, alternating=FALSE),
           x=list(cex=.8, alternating=FALSE)),
         between=list(x=1, y=1))
)
## export.eps(hh("dsgn/figure/cc176-4.eps"))


## restore trellis par
trellis.par.set("superpose.symbol", tpgss.old)





## cc176-ancova-plot.s
## follows dsgn/code/cc176.s

## ANCOVA plot

## trellis has 7 symbols, we need 12
tpgss.old <- trellis.par.get("superpose.symbol")
tpgsl.old <- trellis.par.get("superpose.line")

tpgss <- Rows(tpgss.old, rep(1:3, 4))
tpgss$pch <- rep(levels(cc176$n.treats), 4)

tpgsl <- Rows(tpgsl.old, rep(1:3, each=4))

trellis.par.set("superpose.symbol", tpgss)
trellis.par.set("superpose.line", tpgsl)

cc176$n.c <- interaction(cc176$n.treat, cc176$current)
if.R(r={},
     s=contrasts(cc176$n.c) <- contr.treatment(levels(cc176$n.c)))

cc176.a <-
ancova(wt.d ~ n.c, x=wt.n, data=cc176,
       par.strip.text=list(cex=1.2),
       layout=c(3,5), skip=c(rep(FALSE,12),TRUE,FALSE,TRUE),
       scales=list(cex=1, alternating=FALSE),
       between=list(x=.5, y=c(.5,.5,.5,1.5)))
cc176.b <-
ancova(wt.d ~ wt.n, groups=n.c, data=cc176,
       par.strip.text=list(cex=1.2),
       layout=c(3,5), skip=c(rep(FALSE,12),TRUE,FALSE,TRUE),
       scales=list(cex=1, alternating=FALSE),
       between=list(x=.5, y=c(.5,.5,.5,1.5)))
cc176.c <-
ancova(wt.d ~ n.c + wt.n, data=cc176,
       par.strip.text=list(cex=1.2),
       layout=c(3,5), skip=c(rep(FALSE,12),TRUE,FALSE,TRUE),
       scales=list(cex=1, alternating=FALSE),
       between=list(x=.5, y=c(.5,.5,.5,1.5)))
cc176.d <-
ancova(wt.d ~ wt.n + n.c, data=cc176,
       par.strip.text=list(cex=1.2),
       layout=c(3,5), skip=c(rep(FALSE,12),TRUE,FALSE,TRUE),
       scales=list(cex=1, alternating=FALSE),
       between=list(x=.5, y=c(.5,.5,.5,1.5)))
cc176.e <-
ancova(wt.d ~ wt.n * n.c, data=cc176,
       par.strip.text=list(cex=1.2),
       layout=c(3,5), skip=c(rep(FALSE,12),TRUE,FALSE,TRUE),
       scales=list(cex=1, alternating=FALSE),
       between=list(x=.5, y=c(.5,.5,.5,1.5)))

summary(cc176.a)
summary(cc176.b)
summary(cc176.c)
summary(cc176.d)
summary(cc176.e)


## This is a very busy set of graphs.
## We will tinker with the trellis objects to adjust the appearance.
## Suppress the key.

if.R(s={
cc176.key <- attr(cc176.a,"trellis")$key
attr(cc176.a,"trellis")$key <- NULL
attr(cc176.b,"trellis")$key <- NULL
attr(cc176.c,"trellis")$key <- NULL
attr(cc176.d,"trellis")$key <- NULL
attr(cc176.e,"trellis")$key <- NULL
}
,r={
cc176.key <- attr(cc176.a,"trellis")$legend$right$args$key
attr(cc176.a,"trellis")$legend <- NULL
attr(cc176.b,"trellis")$legend <- NULL
attr(cc176.c,"trellis")$legend <- NULL
attr(cc176.d,"trellis")$legend <- NULL
attr(cc176.e,"trellis")$legend <- NULL
}
)
## Suppress the superpose panel by putting it on its
## own page and then don't display that page.
attr(cc176.a,"trellis")$layout <- c(3,4,1)
attr(cc176.b,"trellis")$layout <- c(3,4,1)
attr(cc176.c,"trellis")$layout <- c(3,4,1)
attr(cc176.d,"trellis")$layout <- c(3,4,1)
attr(cc176.e,"trellis")$layout <- c(3,4,1)
if.R(s={},
     r={
attr(cc176.a,"trellis")$y.between <- attr(cc176.e,"trellis")$y.between[1:3]
attr(cc176.b,"trellis")$y.between <- attr(cc176.b,"trellis")$y.between[1:3]
attr(cc176.c,"trellis")$y.between <- attr(cc176.c,"trellis")$y.between[1:3]
attr(cc176.d,"trellis")$y.between <- attr(cc176.d,"trellis")$y.between[1:3]
attr(cc176.e,"trellis")$y.between <- attr(cc176.e,"trellis")$y.between[1:3]
}
)

## smaller strip labels
attr(cc176.a,"trellis")$par.strip.text$cex <- 1.2
attr(cc176.b,"trellis")$par.strip.text$cex <- 1.2
attr(cc176.c,"trellis")$par.strip.text$cex <- 1.2
attr(cc176.d,"trellis")$par.strip.text$cex <- 1.2
attr(cc176.e,"trellis")$par.strip.text$cex <- 1.2

## larger xlab and ylab
attr(cc176.a,"trellis")$xlab <- list(attr(cc176.a,"trellis")$xlab, cex=1.2)
attr(cc176.b,"trellis")$xlab <- list(attr(cc176.b,"trellis")$xlab, cex=1.2)
attr(cc176.c,"trellis")$xlab <- list(attr(cc176.c,"trellis")$xlab, cex=1.2)
attr(cc176.d,"trellis")$xlab <- list(attr(cc176.d,"trellis")$xlab, cex=1.2)
attr(cc176.e,"trellis")$xlab <- list(attr(cc176.e,"trellis")$xlab, cex=1.2)
attr(cc176.a,"trellis")$ylab <- list(attr(cc176.a,"trellis")$ylab, cex=1.2)
attr(cc176.b,"trellis")$ylab <- list(attr(cc176.b,"trellis")$ylab, cex=1.2)
attr(cc176.c,"trellis")$ylab <- list(attr(cc176.c,"trellis")$ylab, cex=1.2)
attr(cc176.d,"trellis")$ylab <- list(attr(cc176.d,"trellis")$ylab, cex=1.2)
attr(cc176.e,"trellis")$ylab <- list(attr(cc176.e,"trellis")$ylab, cex=1.2)

## modify the main titles
attr(cc176.a,"trellis")$main <- list(paste("a. Separate horizontal lines:", attr(cc176.a,"trellis")$main$main), cex=1.4)
attr(cc176.b,"trellis")$main <- list(paste("b. Identical lines:", attr(cc176.b,"trellis")$main$main), cex=1.4)
attr(cc176.c,"trellis")$main <- list(paste("c,d. Parallel lines:", attr(cc176.c,"trellis")$main$main), cex=1.4)
attr(cc176.d,"trellis")$main <- list(paste("c,d. Parallel lines:", attr(cc176.d,"trellis")$main$main), cex=1.4)
attr(cc176.e,"trellis")$main <- list(paste("e. Separate lines:", attr(cc176.e,"trellis")$main$main), cex=1.4)

## y ~ a                     ## different horizontal line in each group
print(position=c(0,0,  .75, .95), more=FALSE, attr(cc176.a,"trellis"))
## export.eps(hh("dsgn/figure/cc176-5a.eps"))  ## first page only

## y ~ x                     ## constant line across all groups
print(position=c(0,0,  .75, .95), more=FALSE, attr(cc176.b,"trellis"))
## export.eps(hh("dsgn/figure/cc176-5b.eps"))  ## first page only

## y ~ x + a  or  y ~ a + x  ## constant slope, different intercepts
## cc176.c and cc176.d have the same graph
print(position=c(0,0,  .75, .95), more=FALSE, attr(cc176.d,"trellis"))
## export.eps(hh("dsgn/figure/cc176-5d.eps"))  ## first page only

## y ~ x * a  or  y ~ a * x  ## different slopes, and different intercepts
print(position=c(0,0,  .75, .95), more=FALSE, attr(cc176.e,"trellis"))
## export.eps(hh("dsgn/figure/cc176-5e.eps"))  ## first page only

## print a single key
if.R(r={
  grid.newpage()
  cc176.key$x <- .5
  cc176.key$y <- .5
  cc176.key$corner <- c(.5,0)
  invisible(draw.key(cc176.key, draw=TRUE))
}, s={
  frame()
  cc176.key$x <- .5
  cc176.key$y <- .5
  cc176.key$corner <- c(.5,0)
  do.call("key", cc176.key)
  par(mfrow=c(1,1))
})
# export.eps(hh("dsgn/figure/cc176-5key.eps"))

## restore trellis par
trellis.par.set("superpose.symbol", tpgss.old)
trellis.par.set("superpose.line", tpgsl.old)



## tires.r
## tires.s
data(tires)
tires <- tires

tires.aov <- aov(wear ~ car + position + brand, data=tires)

summary(tires.aov)

tapply(tires$wear, tires$car, "mean")
tapply(tires$wear, tires$position, "mean")
tapply(tires$wear, tires$brand, "mean")

if.R(s={
  ## multicomp(tires.aov, method="tukey", focus="car")
  ## multicomp(tires.aov, method="tukey", focus="position")
  multicomp(tires.aov, method="tukey", focus="brand")
  tires.brand.mmc <- multicomp.mmc(tires.aov, method="tukey", focus="brand")
  print(tires.brand.mmc)
},r={
  tires.brand.glht <- glht(tires.aov, linfct=mcp(brand="Tukey"))
  confint(tires.brand.glht)
  tires.mmc.brand <- mmc(tires.aov, linfct=mcp(brand="Tukey"))
  print(tires.mmc.brand)
  plot(tires.mmc.brand, ry=c(10.2, 14.5), x.offset=.8)
})




## filmcoat.s
data(filmcoat)

## interaction plot
filmcoat.2wt <-
  interaction2wt(data=filmcoat, coat ~ temprt+pressure)
if.R(r=print(filmcoat.2wt),
     s=print(filmcoat.2wt, position=c(0.05,0,1,1)))
## export.eps(hh("dsgn/figure/filmcoat.int.eps"))

## means
tapply(filmcoat$coat, filmcoat[,"temprt"], mean)
tapply(filmcoat$coat, filmcoat[,"pressure"], mean)
tapply(filmcoat$coat, filmcoat[,c("temprt","pressure")], mean)

## anova
film.aov1 <- aov(coat ~ temprt*pressure, data=filmcoat)
summary(film.aov1)

contrasts(filmcoat$temprt)
contrasts(filmcoat$pressure)


## split for single degree of freedom contrasts requires the nesting notation
## The order is not obvious, so we try both.  aov2 does temprt/pressure
film.aov2 <- aov(coat ~ temprt/pressure, data=filmcoat)
summary(film.aov2)
summary.lm(film.aov2)$coef      ## needed to get the order used in split
if.R(r=   ## SIMILAR S-Plus/R CHANGES needed for split= elsewhere in this file
     summary(film.aov2,
             split=list("temprt:pressure"=
               list(t.low=c(1,4), t.med=c(2,5), t.high=c(3,6))))
     ,s=
     summary(film.aov2,
             split=list("pressure %in% temprt"=
               list(t.low=c(1,4), t.med=c(2,5), t.high=c(3,6))))

     )

### Sometimes we must look at the generated x matrix
film.aov2x <- aov(coat ~ temprt/pressure, data=filmcoat, x=TRUE)
if.R(r=film.aov2x$x[, "temprtt.low:pressurep.med"],
     s=film.aov2x$x[, "temprtt.lowpressurep.med"])


## split for single degree of freedom contrasts requires the nesting notation
## The order is not obvious, so we try both.  aov3 does pressure/temprt
film.aov3 <- aov(coat ~ pressure/temprt, data=filmcoat)
summary(film.aov3)
if.R(r=
     summary(film.aov3,
             split=list("pressure:temprt"=
               list(p.low=c(1,4), p.med=c(2,5), p.high=c(3,6))))
     , s=
     summary(film.aov3,
             split=list("temprt %in% pressure"=
               list(p.low=c(1,4), p.med=c(2,5), p.high=c(3,6))))
)

## multiple comparisons
## differences from crossed models

## We are saving the results of multicomp() and then plotting and
## printing them separately.  The plot.multicomp function in the
## hh/splus.library gives us control over the placement of the graph
## on the device.  The default placement crowds the very long labels
## in this example.
if.R(r={},
     s={
mcout5.filmcoat <-
  multicomp(film.aov1, focus = "temprt",
            adjust = list(pressure = c("p.low","p.med","p.high") ),
            method = "sim")
mcout5.filmcoat <-
  multicomp.label.change(mcout5.filmcoat, old=".adj1", new="@p.low")
mcout5.filmcoat <-
  multicomp.label.change(mcout5.filmcoat, old=".adj2", new="@p.med")
mcout5.filmcoat <-
  multicomp.label.change(mcout5.filmcoat, old=".adj3", new="@p.high")
mcout5.filmcoat <-
  multicomp.label.change(mcout5.filmcoat, old="-", new=" - ", how.many=1)
old.par <- par(fig=c(0,1,-.4,1))
par(new=FALSE) ## setting fig, incorrectly sets new=TRUE
plot(mcout5.filmcoat,
     plt=c(.3,.9,.1,.9), x.label.adj=1.2)
par(old.par)
par(new=FALSE) ## setting fig, incorrectly sets new=TRUE
## export.eps(hh("dsgn/figure/mcout5.eps"))
mcout5.filmcoat

mcout6.filmcoat <-
  multicomp(film.aov1, focus = "pressure",
            adjust = list(temprt = c("t.low","t.med","t.high") ),
            method = "sim")
mcout6.filmcoat <-
  multicomp.label.change(mcout6.filmcoat, old=".adj1", new="@t.low")
mcout6.filmcoat <-
  multicomp.label.change(mcout6.filmcoat, old=".adj2", new="@t.med")
mcout6.filmcoat <-
  multicomp.label.change(mcout6.filmcoat, old=".adj3", new="@t.high")
mcout6.filmcoat <-
  multicomp.label.change(mcout6.filmcoat, old="-", new=" - ", how.many=1)
old.par <- par(fig=c(0,1,-.4,1)) ## setting fig incorrectly sets new=TRUE
par(new=FALSE) ## setting fig incorrectly sets new=TRUE
plot(mcout6.filmcoat,
     plt=c(.3,.9,.1,.9), x.label.adj=1.1, xrange.include=-12)
par(old.par)
## export.eps(hh("dsgn/figure/mcout6.eps"))
mcout6.filmcoat
})


## filmcoat2.s  ## see also
## follows filmcoat.s

if.R(r={
## separate ANOVA for each pressure
filmcoat.aov.3p <- sapply(levels(filmcoat$pressure),
                          function(i) aov(coat ~ temprt,
                                          data=filmcoat,
                                          subset=(pressure==i)),
                          simplify=FALSE,
                          USE.NAMES=TRUE)
print(lapply(filmcoat.aov.3p, anova))

## Separate multicomp for each pressure
filmcoat.mca.3p <- lapply(filmcoat.aov.3p, glht, linfct=mcp(temprt="Tukey"))

## print and plot multicomps for each pressure, each with its own critical value.
print(lapply(filmcoat.mca.3p, confint))
old.par <- par(mfrow=c(3,1), omd=c(.1,1,0,1))
for (i in levels(filmcoat$pressure))
  plot(filmcoat.mca.3p[[i]], xlim=c(-13,9), main=i, ylim=c(.5,3.5))
par(old.par)


## print and plot multicomps for each pressure, with common critical value.
crit.val <- qtukey(.95, 3, 18, 3)/sqrt(2)
filmcoat.mca.3p.common <- lapply(filmcoat.mca.3p, confint, calpha=crit.val)
print(filmcoat.mca.3p.common)
old.par <- par(mfrow=c(3,1), omd=c(.1,1,0,1))
for (i in levels(filmcoat$pressure))
  plot(filmcoat.mca.3p.common[[i]], xlim=c(-13,9), main=i, ylim=c(.5,3.5))
par(old.par)


# separate MMC plots, with common critical value
filmcoat.mmc.3p <- lapply(filmcoat.aov.3p, mmc, calpha=crit.val)
print(filmcoat.mmc.3p)
old.par <- par(mfrow=c(3,1), omd=c(0,.87,0,1))
for (i in levels(filmcoat$pressure))
  plot(filmcoat.mmc.3p[[i]], ry=c(34,46), main=i, main2="")
par(old.par)


# separate MMC plots for each pressure, with common stderr and confidence-interval width
ResidMS <- function(x)
  summary(x)[[1]]["Residuals","Mean Sq"]
ResidMSAvg <- mean(sapply(filmcoat.aov.3p, ResidMS))

filmcoat.mmc.3pc <- list()
filmcoat.mmc.3pc[["p.low"]] <- mmc(filmcoat.aov.3p[["p.low"]], calpha=crit.val * sqrt(ResidMSAvg/ResidMS(filmcoat.aov.3p[["p.low"]])))
filmcoat.mmc.3pc[["p.med"]] <- mmc(filmcoat.aov.3p[["p.med"]], calpha=crit.val * sqrt(ResidMSAvg/ResidMS(filmcoat.aov.3p[["p.med"]])))
filmcoat.mmc.3pc[["p.high"]] <- mmc(filmcoat.aov.3p[["p.high"]], calpha=crit.val * sqrt(ResidMSAvg/ResidMS(filmcoat.aov.3p[["p.high"]])))

old.par <- par(mfrow=c(3,1), omd=c(0,.87,0,1))
for (i in levels(filmcoat$pressure)) {
  cat("\n", i, "\n",
      paste("calpha = ", round(crit.val, 2),
            ", ResidMSAvg = ", round(ResidMSAvg, 2),
            ", ResidMS ", i, " = ",
            round(ResidMS(filmcoat.aov.3p[[i]]), 2),
            "\n", sep=""),
      "Quantile below is sqrt(ResidMSAvg/ResidMS) * calpha\n",
      "stderr below is based on ResidMS, lower and upper bounds are based on ResidMSAvg\n")
  print(filmcoat.mmc.3pc[[i]])
  plot(filmcoat.mmc.3pc[[i]], ry=c(34,46), main=i, main2=paste("Average Residual MS =", round(ResidMSAvg, 2)))
}
par(old.par)



## separate ANOVA for each temprt
filmcoat.aov.3t <- sapply(levels(filmcoat$temprt),
                          function(i) aov(coat ~ pressure,
                                          data=filmcoat,
                                          subset=(temprt==i)),
                          simplify=FALSE,
                          USE.NAMES=TRUE)
print(lapply(filmcoat.aov.3t, anova))

## Separate multicomp for each temprt
filmcoat.mca.3t <- lapply(filmcoat.aov.3t, glht, linfct=mcp(pressure="Tukey"))

## print and plot multicomps for each temprt, each with its own critical value.
print(lapply(filmcoat.mca.3t, confint))
old.par <- par(mfrow=c(3,1), omd=c(.1,1,0,1))
for (i in levels(filmcoat$temprt))
  plot(filmcoat.mca.3t[[i]], xlim=c(-13,9), main=i, ylim=c(.5,3.5))
par(old.par)


## print and plot multicomps for each temprt, with common critical value.
crit.val <- qtukey(.95, 3, 18, 3)/sqrt(2)
filmcoat.mca.3t.common <- lapply(filmcoat.mca.3t, confint, calpha=crit.val)
print(filmcoat.mca.3t.common)
old.par <- par(mfrow=c(3,1), omd=c(.1,1,0,1))
for (i in levels(filmcoat$temprt))
  plot(filmcoat.mca.3t.common[[i]], xlim=c(-13,9), main=i, ylim=c(.5,3.5))
par(old.par)


# separate MMC plots, with common critical value
filmcoat.mmc.3t <- lapply(filmcoat.aov.3t, mmc, calpha=crit.val)
print(filmcoat.mmc.3t)
old.par <- par(mfrow=c(3,1), omd=c(0,.87,0,1))
for (i in levels(filmcoat$temprt))
  plot(filmcoat.mmc.3t[[i]], ry=c(34,46), main=i, main2="")
par(old.par)


# separate MMC plots for each temprt, with common stderr and confidence-interval width
ResidMS <- function(x)
  summary(x)[[1]]["Residuals","Mean Sq"]
ResidMSAvg <- mean(sapply(filmcoat.aov.3p, ResidMS))

filmcoat.mmc.3tc <- list()
filmcoat.mmc.3tc[["t.low"]] <- mmc(filmcoat.aov.3t[["t.low"]], calpha=crit.val * sqrt(ResidMSAvg/ResidMS(filmcoat.aov.3t[["t.low"]])))
filmcoat.mmc.3tc[["t.med"]] <- mmc(filmcoat.aov.3t[["t.med"]], calpha=crit.val * sqrt(ResidMSAvg/ResidMS(filmcoat.aov.3t[["t.med"]])))
filmcoat.mmc.3tc[["t.high"]] <- mmc(filmcoat.aov.3t[["t.high"]], calpha=crit.val * sqrt(ResidMSAvg/ResidMS(filmcoat.aov.3t[["t.high"]])))

old.par <- par(mfrow=c(3,1), omd=c(0,.87,0,1))
for (i in levels(filmcoat$temprt)) {
  cat("\n", i, "\n",
      paste("calpha = ", round(crit.val, 2),
            ", ResidMSAvg = ", round(ResidMSAvg, 2),
            ", ResidMS ", i, " = ",
            round(ResidMS(filmcoat.aov.3t[[i]]), 2),
            "\n", sep=""),
      "Quantile below is sqrt(ResidMSAvg/ResidMS) * calpha\n",
      "stderr below is based on ResidMS, lower and upper bounds are based on ResidMSAvg\n")
  print(filmcoat.mmc.3tc[[i]])
  plot(filmcoat.mmc.3tc[[i]], ry=c(34,46), main=i, main2=paste("Average Residual MS =", round(ResidMSAvg, 2)))
}
par(old.par)


},
     s={
par(mfrow=c(1,1))

mcout5a.filmcoat <-
  multicomp(film.aov1, focus = "temprt",
            adjust = list(pressure = "p.low"),
            crit.point=mcout5.filmcoat$crit.point)
mcout5a.filmcoat <-
  multicomp.label.change(mcout5a.filmcoat, old=".adj1", new="@p.low")
old.par <- par(fig=c(0,1,-.4,1))
par(new=FALSE) ## setting fig, incorrectly sets new=TRUE
plot(mcout5a.filmcoat,
     plt=c(.3,.9,.1,.9), x.label.adj=1.1, xlim=c(-6,14))
par(old.par)
par(new=FALSE) ## setting fig, incorrectly sets new=TRUE
## export.eps(hh("dsgn/figure/filmcoata.t.p.low.eps"))

mcout5b.filmcoat <-
  multicomp(film.aov1, focus = "temprt",
            adjust = list(pressure = "p.med"),
            crit.point=mcout5.filmcoat$crit.point)
mcout5b.filmcoat <-
  multicomp.label.change(mcout5b.filmcoat, old=".adj1", new="@p.med")
old.par <- par(fig=c(0,1,-.4,1))
par(new=FALSE) ## setting fig, incorrectly sets new=TRUE
plot(mcout5b.filmcoat,
     plt=c(.3,.9,.1,.9), x.label.adj=1.1, xlim=c(-6,14))
par(old.par)
par(new=FALSE) ## setting fig, incorrectly sets new=TRUE
## export.eps(hh("dsgn/figure/filmcoata.t.p.med.eps"))


mcout5c.filmcoat <-
  multicomp(film.aov1, focus = "temprt",
            adjust = list(pressure = "p.high"),
            crit.point=mcout5.filmcoat$crit.point)
mcout5c.filmcoat <-
  multicomp.label.change(mcout5c.filmcoat, old=".adj1", new="@p.high")
old.par <- par(fig=c(0,1,-.4,1))
par(new=FALSE) ## setting fig, incorrectly sets new=TRUE
plot(mcout5c.filmcoat,
     plt=c(.3,.9,.1,.9), x.label.adj=1.1, xlim=c(-6,14))
par(old.par)
par(new=FALSE) ## setting fig, incorrectly sets new=TRUE
## export.eps(hh("dsgn/figure/filmcoata.t.p.high.eps"))



mcout6a.filmcoat <-
  multicomp(film.aov1, focus = "pressure",
            adjust = list(temprt = "t.low"),
            crit.point=mcout6.filmcoat$crit.point)
mcout6a.filmcoat <-
  multicomp.label.change(mcout6a.filmcoat, old=".adj1", new="@t.low")
old.par <- par(fig=c(0,1,-.4,1))
par(new=FALSE) ## setting fig, incorrectly sets new=TRUE
plot(mcout6a.filmcoat,
     plt=c(.3,.9,.1,.9), x.label.adj=1.1, xlim=c(-14,6))
par(old.par)
par(new=FALSE) ## setting fig, incorrectly sets new=TRUE
## export.eps(hh("dsgn/figure/filmcoata.p.t.low.eps"))

mcout6b.filmcoat <-
  multicomp(film.aov1, focus = "pressure",
            adjust = list(temprt = "t.med"),
            crit.point=mcout6.filmcoat$crit.point)
mcout6b.filmcoat <-
  multicomp.label.change(mcout6b.filmcoat, old=".adj1", new="@t.med")
old.par <- par(fig=c(0,1,-.4,1))
par(new=FALSE) ## setting fig, incorrectly sets new=TRUE
plot(mcout6b.filmcoat,
     plt=c(.3,.9,.1,.9), x.label.adj=1.1, xlim=c(-14,6))
par(old.par)
par(new=FALSE) ## setting fig, incorrectly sets new=TRUE
## export.eps(hh("dsgn/figure/filmcoata.p.t.med.eps"))


mcout6c.filmcoat <-
  multicomp(film.aov1, focus = "pressure",
            adjust = list(temprt = "t.high"),
            crit.point=mcout6.filmcoat$crit.point)
mcout6c.filmcoat <-
  multicomp.label.change(mcout6c.filmcoat, old=".adj1", new="@t.high")
old.par <- par(fig=c(0,1,-.4,1))
par(new=FALSE) ## setting fig, incorrectly sets new=TRUE
plot(mcout6c.filmcoat,
     plt=c(.3,.9,.1,.9), x.label.adj=1.1, xlim=c(-14,8))
par(old.par)
par(new=FALSE) ## setting fig, incorrectly sets new=TRUE
## export.eps(hh("dsgn/figure/filmcoata.p.t.high.eps"))







old.par <- par(mar=c(5,4,4,13.5)+.1)
mcout5a.filmcoat.mmc <-
  multicomp.mmc(film.aov1, focus = "temprt",
                adjust = list(pressure = "p.low"),
                crit.point=mcout5.filmcoat$crit.point,
                lmat.rows=2:4,
                ry=c(34,46),
                plot=FALSE)
mcout5a.filmcoat.mmc <-
  multicomp.label.change(mcout5a.filmcoat.mmc, old=".adj1", new="@p.low")
plot(mcout5a.filmcoat.mmc,
     main="temperature at p.low",
     method="simulation")
## export.eps(hh("dsgn/figure/filmcoatb.t.p.low.eps"))

mcout5b.filmcoat.mmc <-
  multicomp.mmc(film.aov1, focus = "temprt",
                adjust = list(pressure = "p.med"),
                crit.point=mcout5.filmcoat$crit.point,
                lmat.rows=2:4,
                ry=c(34,46),
                plot=FALSE)
mcout5b.filmcoat.mmc <-
  multicomp.label.change(mcout5b.filmcoat.mmc, old=".adj1", new="@p.med")
plot(mcout5b.filmcoat.mmc,
     main="temperature at p.med",
     method="simulation")
## export.eps(hh("dsgn/figure/filmcoatb.t.p.med.eps"))

mcout5c.filmcoat.mmc <-
  multicomp.mmc(film.aov1, focus = "temprt",
                adjust = list(pressure = "p.high"),
                crit.point=mcout5.filmcoat$crit.point,
                lmat.rows=2:4,
                ry=c(34,46),
                plot=FALSE)
mcout5c.filmcoat.mmc <-
  multicomp.label.change(mcout5c.filmcoat.mmc, old=".adj1", new="@p.high")
plot(mcout5c.filmcoat.mmc,
     main="temperature at p.high",
     method="simulation")
## export.eps(hh("dsgn/figure/filmcoatb.t.p.high.eps"))


mcout6a.filmcoat.mmc <-
  multicomp.mmc(film.aov1, focus = "pressure",
                adjust = list(temprt = "t.low"),
                crit.point=mcout6.filmcoat$crit.point,
                lmat.rows=5:7,
                ry=c(34,46),
                plot=FALSE)
mcout6a.filmcoat.mmc <-
  multicomp.label.change(mcout6a.filmcoat.mmc, old=".adj1", new="@t.low")
plot(mcout6a.filmcoat.mmc,
     main="pressure at t.low",
     method="simulation")
## export.eps(hh("dsgn/figure/filmcoatb.p.t.low.eps"))

mcout6b.filmcoat.mmc <-
  multicomp.mmc(film.aov1, focus = "pressure",
                adjust = list(temprt = "t.med"),
                crit.point=mcout6.filmcoat$crit.point,
                lmat.rows=5:7,
                ry=c(34,46),
                plot=FALSE)
mcout6b.filmcoat.mmc <-
  multicomp.label.change(mcout6b.filmcoat.mmc, old=".adj1", new="@t.med")
plot(mcout6b.filmcoat.mmc,
     main="pressure at t.med",
     method="simulation")
## export.eps(hh("dsgn/figure/filmcoatb.p.t.med.eps"))

mcout6c.filmcoat.mmc <-
  multicomp.mmc(film.aov1, focus = "pressure",
                adjust = list(temprt = "t.high"),
                crit.point=mcout6.filmcoat$crit.point,
                lmat.rows=5:7,
                ry=c(34,46),
                plot=FALSE)
mcout6c.filmcoat.mmc <-
  multicomp.label.change(mcout6c.filmcoat.mmc, old=".adj1", new="@t.high")
plot(mcout6c.filmcoat.mmc,
     main="pressure at t.high",
     method="simulation")
## export.eps(hh("dsgn/figure/filmcoatb.p.t.high.eps"))

par(old.par)


## plot 3 on the same page
mcout5a.filmcoat.mmc <-
  multicomp.label.change(mcout5a.filmcoat.mmc, new="", old="@p.low")
mcout5b.filmcoat.mmc <-
  multicomp.label.change(mcout5b.filmcoat.mmc, new="", old="@p.med")
mcout5c.filmcoat.mmc <-
  multicomp.label.change(mcout5c.filmcoat.mmc, new="", old="@p.high")
##
mcout6a.filmcoat.mmc <-
  multicomp.label.change(mcout6a.filmcoat.mmc, new="", old="@t.low")
mcout6b.filmcoat.mmc <-
  multicomp.label.change(mcout6b.filmcoat.mmc, new="", old="@t.med")
mcout6c.filmcoat.mmc <-
  multicomp.label.change(mcout6c.filmcoat.mmc, new="", old="@t.high")
##
par(mfrow=c(1,3), mar=c(30,3,10,6)+.1)
plot(mcout6a.filmcoat.mmc,  main="pressure at t.low", main2="")
plot(mcout6b.filmcoat.mmc,  main="pressure at t.med",  main2="")
plot(mcout6c.filmcoat.mmc, main="pressure at t.high", main2="")
mtext("filmcoat: simple effects of pressure", side=3, outer=TRUE, cex=1.8, line=-3)
##
plot(mcout5a.filmcoat.mmc,  main="temperature at p.low",  main2="")
plot(mcout5b.filmcoat.mmc,  main="temperature at p.med",  main2="")
plot(mcout5c.filmcoat.mmc, main="temperature at p.high", main2="")
mtext("filmcoat: simple effects of temperature", side=3, outer=TRUE, cex=1.8, line=-3)
par(mar=c(5,4,4,2)+.1) ## two separate lines for par() required.
par(mfrow=c(1,1))      ## two separate lines for par() required.
})








## gunload.s
## Analysis of Gunload Data
data(gunload)
gunload.aov <- aov(rounds ~ method*group + Error((team %in% group)/method),
                   data=gunload, qr=TRUE)
summary(gunload.aov)

## model.tables for models that have an Error() term require a bug fix
## in S-Plus 6.2 and earlier.  The library("HH") fixes that bug.
## Another bug requires you to set se=FALSE for models with an Error() term.

model.tables(gunload.aov, type="means", se=FALSE)

if.R(s=
     print(position=c(0,.1,1,1),
           t(bwplot(team  ~ rounds | group*method, data=gunload,
                    xlab=list(cex=1.4),
                    ylab=list("team %in% group", cex=1.4),
                    par.strip.text=list(cex=1.4),
                    scales=list(x=list(cex=1),y=list(cex=1.2))
                    ))
           )
     ,r=
     bwplot(rounds  ~ team | group*method, data=gunload,
            ylab=list(cex=1.4),
            xlab=list("team %in% group", cex=1.4),
            par.strip.text=list(cex=1.4),
            scales=list(x=list(cex=1),y=list(cex=1.2))
            )
     )






## turkey.factors.s
## follows Ch06-oway.r  section labeled oway/code/turkey-oway.s
data(turkey)
turkey <- turkey
turkey[c(1,7,13,19,25),]

turkey$trt.vs.control <- factor(rep(c("control","treatment"), c(6,24)))
contrasts(turkey$trt.vs.control) <- c(4,-1)

turkey$additive <- factor(rep(c("control","A","B"), c(6,12,12)),
                          levels=c("control","A","B"))
contrasts(turkey$additive) <- c(0,1,-1)

turkey$amount <- factor(rep(c(0,1,2,1,2), c(6,6,6,6,6)))
contrasts(turkey$amount) <- c(0,1,-1)

turkey[c(1,7,13,19,25),]



## turkey.aov2.s
## follows dsgn/code/turkey.factors.s

turkey3.aov <- aov(wt.gain ~ trt.vs.control / (additive*amount),
                   data=turkey, x=TRUE)
summary(turkey3.aov)


## turkey.means.s
## follows dsgn/code/turkey.factors.s

tapply(turkey$wt.gain,
       turkey[,c("additive","amount")],
       mean)


## turkey.f2.s
## turkey.f2.r
## follows dsgn/code/turkey.factors.s

## bwplot(additive ~ wt.gain | amount,
##       data=turkey, layout=c(3,1))

additive.rev <- ordered(turkey$additive, rev(levels(turkey$additive)))

print(position=c(0,.3, 1,1),
bwplot(additive.rev ~ wt.gain | amount,
       data=turkey, layout=c(3,1),
       scales=list(cex=1.4, alternating=FALSE),
       xlab=list(cex=1.4),
       ylab=list("additive",cex=1.4),
       par.strip.text=list(cex=1.4),
       strip=function(...)
       strip.default(..., style = 1,
                     strip.names = c(TRUE, TRUE)))
)
## export.eps(hh("dsgn/figure/turkey.f2h.eps"))


if.R(r={
print(position=c(.25,0, .75,1),
  bwplot(wt.gain ~ additive | amount,
         data=turkey, layout=c(1,3),
         scales=list(cex=1.4, alternating=FALSE),
         between=list(y=1),
         ylab=list(cex=1.4),
         xlab=list("additive",cex=1.4),
         par.strip.text=list(cex=1.4),
         strip=function(...)
         strip.default(..., style = 1,
                       strip.names = c(TRUE, TRUE)))
      )
}, s={
  print(position=c(.25,0, .75,1),
        t(bwplot(additive ~ wt.gain | amount,
                 data=turkey, layout=c(1,3),
                 scales=list(cex=1.4, alternating=FALSE),
                 between=list(y=1),
                 xlab=list(cex=1.4),
                 ylab=list("additive",cex=1.4),
                 par.strip.text=list(cex=1.4),
                 strip=function(...)
                 strip.default(..., style = 1,
                               strip.names = c(TRUE, TRUE)))
          ))
})
## export.eps(hh("dsgn/figure/turkey.f2v.eps"))


if.R(r={
print(position=c(0,.3, 1,1),
  bwplot(wt.gain ~ additive.rev | amount,
         data=turkey, layout=c(3,1),
         scales=list(cex=1.4, alternating=FALSE),
         ylab=list(cex=1.4),
         xlab=list("additive",cex=1.4),
         par.strip.text=list(cex=1.4),
         strip=function(...)
         strip.default(..., style = 1,
                       strip.names = c(TRUE, TRUE)))

)
}, s={
  print(position=c(0,.3, 1,1),
        t(bwplot(additive.rev ~ wt.gain | amount,
                 data=turkey, layout=c(3,1),
                 scales=list(cex=1.4, alternating=FALSE),
                 xlab=list(cex=1.4),
                 ylab=list("additive",cex=1.4),
                 par.strip.text=list(cex=1.4),
                 strip=function(...)
                 strip.default(..., style = 1,
                               strip.names = c(TRUE, TRUE)))
          ))
})
## export.eps(hh("dsgn/figure/turkey.f2v-h.eps"))


print(position=c(.25,0, .75,1),
bwplot(additive ~ wt.gain | amount,
       data=turkey, layout=c(1,3),
       scales=list(cex=1.4, alternating=FALSE),
       between=list(y=1),
       xlab=list(cex=1.4),
       ylab=list("additive",cex=1.4),
       par.strip.text=list(cex=1.4),
       strip=function(...)
       strip.default(..., style = 1,
                     strip.names = c(TRUE, TRUE)))
)
## export.eps(hh("dsgn/figure/turkey.f2h-v.eps"))




## turkey.aov3.s
## follows dsgn/code/turkey.aov2.s

match(dimnames(coef(summary.lm(turkey3.aov)))[[1]],
      dimnames(turkey3.aov$x)[[2]])

tmp1 <- coef(summary.lm(turkey3.aov))
dimnames(tmp1)[[1]][3:5]
dimnames(tmp1)[[1]][3:5] <- c("additive","amount","additive:amount")
dimnames(tmp1)[[1]][3:5]

tmp2 <- turkey3.aov$x[,c(1,2,4,8,12)]
dimnames(tmp2)[[2]][3:5] <- c("additive","amount","additive:amount")

round(tmp1, 4)
tmp2[c(1,7,13,19,25),]


## turkey.aov3-match.s
## follows dsgn/code/turkey.aov3.s



                    summary.lm(turkey3.aov)
               coef(summary.lm(turkey3.aov))
      dimnames(coef(summary.lm(turkey3.aov)))
      dimnames(coef(summary.lm(turkey3.aov)))[[1]]

               turkey3.aov$x
      dimnames(turkey3.aov$x)
      dimnames(turkey3.aov$x)[[2]]

match(dimnames(coef(summary.lm(turkey3.aov)))[[1]],
      dimnames(turkey3.aov$x)[[2]])




## contrasts.s
A <- rep(c("r","s","t"), rep(4, 3))
B <- rep(c("w","x","y","z"), 3)
BwA <- letters[3:14]
AB <- paste(A, B, sep=".")
BwA <- ordered(BwA, levels=BwA)

## y <- round(rnorm(12), 2)
y <- c(-0.02, 1.19, -0.02, 0.23,
       0.67, 1.95, -0.71, -0.40,
       -0.56, 0.01, 0.13, 1.19)
y

abc <- data.frame(A, B, AB, BwA, y, row.names=AB)
abc

abc.oneway <- ## one-way
  if.R(s=
       design.table(abc[, c("y","A")])
       ,r= ## this also works in S-Plus
       matrix(abc$y, 4, 3,
              dimnames=list(1:4, abc$A[c(1,5,9)]))
       )
abc.oneway

abc.crossed <- ## crossed
  if.R(s=
       design.table(abc[, c("y","A","B")])
       ,r= ## this also works in S-Plus
       matrix(abc$y, 3, 4, byrow=TRUE,
              dimnames=list(abc$A[c(1,5,9)], abc$B[1:4]))
       )
abc.crossed

abc.nested <- ## nested
  if.R(s=
       design.table(abc[, c("y","A","BwA")])
       ,r= ## this also works in S-Plus
       matrix(c(abc$y[1:4],rep(NA,12),abc$y[5:8],rep(NA,12),abc$y[9:12]),
              3, 12, byrow=TRUE,
              dimnames=list(abc$A[c(1,5,9)],abc$BwA))
       )
abc.nested

abc.double.indexed <- ## doubly-indexed
  if.R(s=
       design.table(abc[, c("y","AB")])
       ,r= ## this also works in S-Plus
       data.matrix(abc)[,"y"]
       )
abc.double.indexed

## one-way
y.A.aov <- aov(y ~ A, data=abc)
anova(y.A.aov)
model.tables(y.A.aov)
resid(y.A.aov)

if.R(s=
     twoway(abc.crossed, trim=0) ## crossed
     ,r=
     {} ## no equivalent for means in R
)
y.ApB.aov <- aov(y ~ A+B, data=abc)
anova(y.ApB.aov)
model.tables(y.ApB.aov)
resid(y.ApB.aov)

if.R(s=
     twoway(design.table(abc[, c("y","A","BwA")]), trim=0) ## nested
     ,r=
     {} ## no equivalent for means in R
     )
y.ABwA.aov <- aov(y ~ A/BwA, data=abc)
anova(y.ABwA.aov)
model.tables(y.ABwA.aov)
resid(y.ABwA.aov)

## doubly-indexed
y.AB.aov <- aov(y ~ AB, data=abc)
anova(y.AB.aov)
model.tables(y.AB.aov)
resid(y.AB.aov)



## short function name for contrasts with redundant columns
ct <- function(f, contrasts=FALSE) contr.treatment(f, contrasts=FALSE)

## one-way
model.matrix(~A, data=abc, contrasts=list(A=ct))
model.matrix(~A, data=abc, contrasts=list(A=contr.treatment))
model.matrix(~A, data=abc, contrasts=list(A=contr.sum))
model.matrix(~A, data=abc, contrasts=list(A=contr.helmert))
model.matrix(~A, data=abc, contrasts=list(A=contr.poly))


A <- factor(A)
contrasts(A) <- contr.treatment(levels(A))
contrasts(A)
lm1 <- lm(y ~ A, x=TRUE)
lm1$x
summary(lm1)$coef
anova(lm1)

contrasts(A, how.many=3) <- contr.treatment(levels(A), contrasts=FALSE)
contrasts(A)
lm2 <- lm(y ~ A, x=TRUE, singular.ok=TRUE)
lm2$x
summary(lm2)$coef
anova(lm2)

aov2 <- aov(y ~ A, x=TRUE)
aov2$x
anova(aov2)
summary.lm(aov2)$coef

contrasts(A, how.many=3) <- contr.sum(levels(A), contrasts=FALSE)
contrasts(A)
aov3 <- aov(y ~ A, x=TRUE)
aov3$x
anova(aov3)
summary.lm(aov3)$coef


## main effects
A <- factor(A)
contrasts(abc$A, how.many=3) <- contr.sum(levels(A), FALSE)
contrasts(abc$A)
B <- factor(B)
contrasts(abc$B, how.many=4) <- contr.sum(levels(B), FALSE)
contrasts(abc$B)
model.matrix(~A+B, data=abc)


## main effects
model.matrix(~A+B, data=abc, contrasts=list(A=ct, B=ct))

## doubly indexed
model.matrix(~A:B, data=abc, contrasts=list(A=ct, B=ct))

if (if.R(r=TRUE,
         s=(version$major==6 && version$minor >= 2) || version$major>6))
  { ## crashes old S-Plus releases 4.5 and 6.0.3 rel 2
    model.matrix(~A*B, contrasts=list(A=ct, B=ct))
  } else {
    cbind(model.matrix(~A+B, data=abc, contrasts=list(A=ct, B=ct)),
          model.matrix(~-1+A:B, data=abc, contrasts=list(A=ct, B=ct)))
  }

## nesting
cbind(model.matrix(~A, data=abc, contrasts=list(A=ct)),
      model.matrix(~-1+A:BwA, data=abc, contrasts=list(A=ct, BwA=ct)))



## start over
A <- rep(c("r","s","t"), rep(4, 3))
B <- rep(c("w","x","y","z"), 3)
BwA <- letters[3:14]
AB <- paste(A, B, sep=".")
BwA <- ordered(BwA, levels=BwA)

A <- factor(A)
contrasts(A) <- contr.sum(levels(A))
contrasts(A)

B <- factor(B)
contrasts(B) <- contr.sum(levels(B))
contrasts(B)

BwA <- factor(BwA)
contrasts(BwA) <- contr.sum(levels(BwA))
contrasts(BwA)

AB <- factor(AB)
contrasts(AB) <- contr.sum(levels(AB))
contrasts(AB)


abc <- data.frame(A, B, AB, BwA, y, row.names=AB)
abc


## main effects
model.matrix(~A+B, data=abc)

## doubly indexed
model.matrix(~AB, data=abc)

## crossing
model.matrix(~A*B, data=abc)

## nesting
model.matrix(~A/B, data=abc)




## vulcan.s
## "A completely randomised 5x3x4 factorial design.
##    (Davies,O.L: Design and Analysis of Industrial Experiments
##     Oliver and Boyd 1954, page 291.)
##  Wear resistance of vulcanised rubber
##  Treatment factors:
##  A - 5 qualities of filler;
##  B - 3 methods of pretreatment of the rubber;
##  C - 4 qualities of raw rubber.
##  Only one replicate, thus 3 factor interaction used for error."

## We got the example from the Genstat Manual.

data(vulcan)
vulcan <- vulcan

vulcan.aov <- aov(wear ~ (filler + pretreat + raw)^2, data=vulcan)
anova(vulcan.aov)
model.tables(vulcan.aov, type="mean", se=TRUE)

print(position=c(0,0,.98,1),
interaction2wt(wear ~ filler+pretreat+raw, data=vulcan,
               main=list(cex=1.6,"wear: main effects and 2-way interactions"),
               par.strip.text=list(cex=.8),
               scales=list(
                 x=list(cex=1),
                 y=list(cex=1, alternating=FALSE)))
)
## export.eps(hh("dsgn/figure/vulcan.eps"))



## Version of the plot that appears in HH page 423.
vulcan$pretreat <- unpositioned(vulcan$pretreat)
vulcan$raw <- unpositioned(vulcan$raw)

tmp <-
interaction2wt(wear ~ filler+pretreat+raw, data=vulcan,
               main=list(cex=1.6,"wear: main effects and 2-way interactions"),
               par.strip.text=list(cex=1.1),
               x.factor.line=0,
               x.factor.cex=.8,
               scales=list(
                 x=list(cex=1),
                 y=list(cex=1, alternating=FALSE)))
## export.eps(hh("dsgn/figure/interaction.trellis.color.eps"))
## export.eps(hh("dsgn/figure/interaction.trellis.eps"))

position(vulcan$pretreat) <- c(2,3,4)
position(vulcan$raw) <- (1:4)+.5
tmp



## tires.latin.s
## tires.latin.2008.s
## follows dsgn/code/tires.s
##
## This is the assignment intended in HH Exercise 13.6
## At the time it was written S-Plus defaulted to Helmert contrasts
## The residual sum of squares in tr1.lm is 0.
contrasts(tires$car)      <- contr.helmert(4)
contrasts(tires$position) <- contr.helmert(4)
contrasts(tires$brand)    <- contr.helmert(4)

tires.aov <- aov(wear ~ car + position + brand, data=tires, x=TRUE)
tires.rc.aov <- aov(wear ~ car * position, data=tires, x=TRUE)

t(tires.aov$x[,8:10])
t(tires.rc.aov$x[,8:16])
tr1.lm <- lm(tires.aov$x[,8] ~ tires.rc.aov$x[,8:16])
anova(tr1.lm)


## R and S-Plus currently both default to treatment contrasts.
## The residual sum of squares in tr1t.lm is not 0.
contrasts(tires$car)      <- contr.treatment(4)
contrasts(tires$position) <- contr.treatment(4)
contrasts(tires$brand)    <- contr.treatment(4)


tirest.aov <- aov(wear ~ car + position + brand, data=tires, x=TRUE)
tirest.rc.aov <- aov(wear ~ car * position, data=tires, x=TRUE)

t(tirest.aov$x[,8:10])
t(tirest.rc.aov$x[,8:16])
tr1t.lm <- lm(tirest.aov$x[,8] ~ tirest.rc.aov$x[,8:16])
anova(tr1t.lm)




## -rwx------+   1 rmh None    80272 2004-06-08 16:18 dsgn.tex
