library(HH)

#### conc/code/simple-pdf.s
pp <- ppoints(101, a=0.5)
zz.norm <- qnorm(pp)
dd.norm <- dnorm(zz.norm)
xx.t2   <- qt(pp, df=2)
dd.t2   <- dt(xx.t2, df=2)

ratio <- dd.norm[51] / dd.t2[51]
xx.normr <- zz.norm * ratio
dd.normr <- dd.norm / ratio

xyplot(dd.norm ~ zz.norm)
xyplot(dd.normr ~ xx.normr)
xyplot(dd.t2 ~ xx.t2)

zz.norm.extra <- seq(3, 9, .5)
dd.norm.extra <- dnorm(zz.norm.extra)
xx.normr.extra <- zz.norm.extra * ratio
dd.normr.extra <- dd.norm.extra / ratio


xx.normr <- c(rev(-xx.normr.extra), xx.normr, xx.normr.extra)
dd.normr <- c(rev( dd.normr.extra), dd.normr, dd.normr.extra)


xx.negskew <- c(xx.t2[1:51],xx.normr[52:101],xx.normr.extra)
dd.negskew <- c(dd.t2[1:51],dd.normr[52:101],dd.normr.extra)

xyplot(dd.negskew ~ xx.negskew)

type <- c("negatively skewed",
          "symmetric",
          "positively skewed")
type <- ordered(type, levels=type)

conc <- data.frame(x=c(xx.negskew, xx.normr, rev(-xx.negskew)),
                   f.x=c(dd.negskew, dd.normr, rev(dd.negskew)),
                   type=rep(type, c(114,127,114)))

## these are scaled so they are all densities
print(position=c(0,0,1,.8),
      xyplot(f.x ~ x | type, data=conc,
             ylab=list("f(x)", cex=1.6),
             xlab=list(cex=1.6),
             par.strip.text=list(cex=1.4),
             scales=list(cex=1.2, alternating=FALSE),
             layout=c(3,1), type="l",
             between=list(x=1),
             panel=function(...) {
               panel.xyplot(...)
               panel.abline(h=0, v=0, lty=2)
             })
)
## export.eps(hh("conc/figure/skewdens.eps"))



xx <- seq(-3,6,.025)
dd <- (dnorm(xx, mean=0, sd=1) + dnorm(xx, mean=2.5, sd=1.1))/2

## bimodal density
xyplot(dd ~ xx, type="l",
       panel=function(...) {
         panel.xyplot(...)
         panel.abline(h=0, lty=2)
       },
       ylab=list("f(x)", cex=1.6, adj=.55),
       xlab=list("x", cex=1.6),
       scales=list(cex=1.4)
       )
## export.eps(hh("conc/figure/bimodal.eps"))


pr <- (pnorm(c(2,4), mean=0, sd=1) + pnorm(c(2,4), mean=2.5, sd=1.1))/2

## bimodal density shaded
xyplot(dd ~ xx, type="l",
       panel=function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(h=0, lty=2)
         if.R(s=polygon(x=x[c(201,201:281,281)],
                y=c(0,y[201:281],0),
                col=80)
              ,
              r=grid.polygon(x=x[c(201,201:281,281)],
                y=c(0,y[201:281],0),
                gp=gpar(fill="80"), default.units="native")
              )
       },
       ylab=list("f(x)", cex=1.6, adj=.55),
       xlab=list("x", cex=1.6),
       scales=list(cex=1.4),
       main=list(paste("Prob(2 < X < 4) =", round(diff(pr),3)), cex=1.5)
       )
## export.eps(hh("conc/figure/bimodal.shade.eps"))



## quartiles
pp <- ppoints(101, a=1)
pp[101] <- .999
q.pp <- qf(pp, df1=3, df2=36)
if.R(s={q.pp[1] <- 0},  ## bug in S-Plus 6.1.2
     r={})
dd <- df(q.pp, df1=3, df2=36)

if.R(s=
     xyplot(dd ~ q.pp, type="l",
            panel=function(x, y, ...) {
              polygon(x=x[c(1,1:26,26)],
                      y=c(0,y[1:26],0),
                      col=80)
              polygon(x=x[c(51,51:76,76)],
                      y=c(0,y[51:76],0),
                      col=80)
              axis(1, at=x[c(26,51,76)], labels=c("Q1","M","Q3"))
              axis(1, at=x[c(26,51,76)], line=2, ticks=FALSE,
                   labels=round(x[c(26,51,76)], 2))
              axis(1, at=x[51], line=2, ticks=FALSE,
                   labels=round(x[51], 2))
              panel.xyplot(x, y, ...)
              panel.abline(h=0, lty=2)
            },
            ylab=list("density", cex=1.6),
            xlab=list("quantile", cex=1.6),
            scales=list(cex=1.4),
            main=list("Quartiles of F(3,36)", cex=1.5)
            )
     ,r=
     xyplot(dd ~ q.pp, type="l",
            par.settings = list(clip = list(panel = "off")),
            panel=function(x, y, ...) {
              grid.polygon(x=x[c(1,1:26,26)],
                           y=c(0,y[1:26],0),
                           gp=gpar(fill="80"), default.units="native")
              grid.polygon(x=x[c(51,51:76,76)],
                           y=c(0,y[51:76],0),
                           gp=gpar(fill="80"), default.units="native")
              panel.axis("bottom", at=x[c(26,51,76)], labels=c("Q1","M","Q3"),
                         outside=TRUE,
                         half=FALSE,
                         rot=0)
              panel.axis("bottom", at=x[c(26,51,76)],
                         tck = 5, line.col = "transparent",
                         labels=round(x[c(26,51,76)], 2),
                         outside=TRUE,
                         half=FALSE,
                         rot=0)
              panel.xyplot(x, y, ...)
              panel.abline(h=0, lty=2)
            },
            ylab=list("density", cex=1.6),
            xlab=list("quantile", cex=1.6),
            scales=list(cex=1.4),
            main=list("Quartiles of F(3,36)", cex=1.5)
            )
     )
## export.eps(hh("conc/figure/quartiles.eps"))



## another skewness example

pp <- ppoints(101, a=0.5)
zz.norm <- qnorm(c(seq(.0010,.0045,.0005), pp, 1-rev(seq(.0010,.0045,.0005))))
dd.norm <- dnorm(zz.norm)

xx.f3.36 <- qf(pp, df1=3, df2=36)
dd.f3.36 <- df(xx.f3.36, df1=3, df2=36)
xx.f3.36 <- c(0,xx.f3.36)
dd.f3.36 <- c(0,dd.f3.36)

ratio <- max(dd.f3.36) / max(dd.norm)
xx.normr <- zz.norm / ratio
dd.normr <- dd.norm * ratio

type <- c("negatively skewed",
          "symmetric",
          "positively skewed")
type <- ordered(type, levels=type)

skew <- data.frame(x=c(rev(-(xx.f3.36-2.5)), xx.normr, xx.f3.36-2.5),
                   f.x=c(rev(dd.f3.36), dd.normr, dd.f3.36),
                   type=rep(type, c(102,117,102)))

## these are scaled so they are all densities
print(position=c(0,0,1,.8),
      xyplot(f.x ~ x | type, data=skew,
             ylab=list("f(x)", cex=1.6, adj=.4),
             xlab=list(cex=1.6),
             par.strip.text=list(cex=1.4),
             scales=list(cex=1.2, alternating=FALSE),
             layout=c(3,1), type="l",
             between=list(x=1),
             panel=function(...) {
               panel.xyplot(...)
               panel.abline(h=0, v=0, lty=2)
             })
)
## export.eps(hh("conc/figure/skewdens2.eps"))


#### conc/code/tv-graphs.s
#### conc/code/tv-graphs.r
data(tv)

## frequency table
tmp <- as.matrix(table(cut(tv$male.life.exp, breaks=seq(49.5,79.5,5))))
dimnames(tmp)[[2]] <- "frequency"
tmp

## old-style graphics
hist(tv$male.life.exp,
     breaks=c(49.5,54.5,59.5,64.5,69.5,74.5,79.5),
     plot=TRUE,
     probability=FALSE,
     include.lowest=FALSE,
     xlab = "male life expectancy")

## trellis graphics
histogram( ~ male.life.exp, data = tv,
          breaks=seq(49.5, 79.5, 5), type="count", col=80,
          scales=list(cex=1.5), xlab=list(cex=1.5), ylab=list(cex=1.5))
## export.eps(hh("conc/figure/tv-hist.eps"))

## stem-and-leaf
if.R(s=
stem(tv$male.life.exp, nl=5, scale=-1, twodig=FALSE, depth=TRUE)
,r=
stem(tv$male.life.exp)
)

#### conc/code/skew.s
sym <- rnorm(100)

neg <- sym[sym<0]
pos <- sym[sym>0]

neg.skew <- c(-(neg^2), pos^.5)

pos.skew <- c(-((-neg)^.5), pos^2)

skew.levels <- c("neg.skew", "sym", "pos.skew")
skew.df <- data.frame(y=c(neg.skew, sym, pos.skew),
                      dist=ordered(rep(skew.levels, c(100,100,100)),
                        levels=skew.levels))
                      
print(position=c(0,.3, 1,1),
bwplot(dist ~ y, data=skew.df,
       xlab="", scales=list(x=list(cex=1.5), y=list(cex=1.5)))
)
## export.eps(hh("conc/figure/skew.eps"))

#### conc/code/corrscat.s
## Bivariate Normal distribution---scatterplot at various correlations

x <- rnorm(100)
e <- rnorm(100)

## One per page, cycle through pages with Ctrl-PageUp and Ctrl-PageDown
old.par <- par(pty="s")
for (r in seq(-1, 1, .1)) {
  y <- r*x + (1-r^2)^.5 * e 
  plot(y ~ x, main=paste("correlation =",r), ylim=c(-3.5,3.5), lab=c(5,7,7))
}
par(old.par)

## trellis plot, at select correlations
corr.data <- data.frame(x=rep(x,7), y=rep(0,7*100), corr=rep(0,7*100))
r <- c(-1, -.9, -.5, 0, .5, .9, 1)
for (i in seq(along=r)) {
  corr.data$y[(i-1)*100 + 1:100] <- r[i]*x + (1-r[i]^2)^.5 * e 
  corr.data$corr[(i-1)*100 + 1:100] <- r[i]
}
corr.data$corr <- factor(corr.data$corr)

xyplot(y ~ x | corr, layout=c(7,1), data=corr.data,
       aspect=1, xlim=c(-3.5,3.5), ylim=c(-3.5,3.5),
       par.strip.text=list(cex=1.3), cex=.5,
       strip = function(...) strip.default(..., strip.names=c(TRUE,TRUE), style=1),
       scales=list(x=list(cex=1), y=list(cex=1), alternating=FALSE)
       )
## export.eps(hh("conc/figure/corr.eps"))

#### conc/code/bivnorm.s
## Bivariate Normal density in 3-D space with various viewpoints.
## Based on the function   example.draping2   in the trellis library.

example.bivariate.normal <-
  function(rho=0, layout.in=c(3,3),
           lwd.in=.1,
           col.regions.in=trellis.par.get("regions")$col)
{
  old.par <- par(lwd=lwd.in)
  on.exit(par(old.par))
  x <- seq(-2, 2, length=33)
  y <- x
  fxy <- 1/(2*pi*sqrt(1-rho^2)) *
    exp(-.5/(1-rho^2) * outer(x^2, y^2, "+") - 2*rho*outer(x,y,"*"))
  angle <- c(22.5, 67.5, 112.5, 337.5, 157.5, 292.5, 247.5, 202.5)
  Angle <- rep(angle, rep(length(fxy), 8))
  Viewing.Angle <- ordered(Angle, angle)
  wireframe(rep(fxy, 8) ~ rep(x[row(fxy)], 8) * rep(y[col(fxy)], 8) |
            Viewing.Angle,
            panel = function(x, y, subscripts, z, angle, ...)
            {
              w <- unique(angle[subscripts])
              panel.wireframe(x=x, y=y, subscripts=subscripts, z=z,
                              screen = list(z = w, x = -60, y = 0), ...)
            },
            angle = Angle, ## this is how to pass down external element
            strip = function(...)
            strip.default(..., strip.names = TRUE, style = "1"),
            skip = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
            drape = TRUE, layout = layout.in, distance = 0.3,
            main = paste("Bivariate Normal, rho=", rho),
            xlab = list("x", cex = 0.4),
            ylab = list("y", cex = 0.4),
            zlab = list("f(x,y)", cex = 0.4),
            col.regions=col.regions.in)
}

## example.bivariate.normal()                   # boring, with rho=0

example.bivariate.normal(.7)                 # all views on one page
### export.eps is not recommended for this example.
### The minimum lwd parameter is too thick on the graphsheet
##
### We recommend using a the postscript driver directly.
## black and white
## trellis.device(postscript, file=hh("conc/figure/bivnorm.eps"), color=FALSE)
## strip.background0()
## example.bivariate.normal(.7,                 # all views on one page
##                          col.regions.in=rep(82:106, rep(4,25)))
## dev.off()
##
## color
## trellis.device(postscript, file=hh("conc/figure/bivnorm-color.eps"), color=TRUE)
## strip.background0()
## example.bivariate.normal(.7)                 # all views on one page
## dev.off()


example.bivariate.normal(.7, layout=c(1,1))  # each view on its own page
## One per page, cycle through pages with Ctrl-PageUp and Ctrl-PageDown
##
## This is the plot from which the figure in the book is taken.
## We use just the Viewing.Angle=112.5 panel.
##
### We recommend using a the postscript driver directly.
## black and white
## trellis.device(postscript, onefile=FALSE, print.it=FALSE, color=FALSE)
## strip.background0()
## example.bivariate.normal(.7, layout=c(1,1),  # One panel per page
##                          col.regions.in=rep(82:106, rep(4,25)))
## dev.off()
## ## manually rename ps.out.003.ps to hh("conc/figure/bivnorm1125.eps")
##
## color
## trellis.device(postscript, onefile=FALSE, print.it=FALSE, color=TRUE)
## strip.background0()
## example.bivariate.normal(.7, layout=c(1,1))  # One panel per page
## dev.off()
## ## manually rename ps.out.003.ps to hh("conc/figure/bivnorm1125-color.eps")


## for (rho in seq(-.9,.9,.1))                  # one page for each rho
##   print(example.bivariate.normal(rho))

if.R(s={
  #### conc/code/bivnorm-rotate.s
  ## 3-D rotating Bivariate Normal density
  brush.bivariate.normal <- function(rho=0) {
    x <- seq(-2, 2, length=33)
    y <- x
    fxy <- 1/(2*pi*sqrt(1-rho^2)) *
      exp(-.5/(1-rho^2) * outer(x^2, y^2, "+") - 2*rho*outer(x,y,"*"))
    brush(cbind(x=rep(x,33), y=as.vector(t(matrix(y,33,33))), z=as.vector(fxy)),
          hist=TRUE)
  }
  brush.bivariate.normal(.7)
}, r={}
     )


     
#### conc/code/normpdf.s
norm.setup(xlim=c(75,125), mean=100, se=5)  ## setup region
norm.curve(100, 5, 100+5*(1.645), z=seq(-5, 5, .1), shade="left", col=80) ## draw plot
## export.eps(hh("conc/figure/normpdf.eps"))

#### conc/code/conc.oc.s
#### conc/code/conc.oc.r
### Operating Characteristics and Power Curves

oc.data <- data.frame(mu=seq(9.8, 11, .01))


if.R(r={
### 1

print(
xyplot(pnorm(41.644 - 4*mu) ~ mu, data=oc.data,
##       main=paste("1. Operating Characteristics Curve:  beta(mu) =",
##         "Probability of Retaining H0 for specified value of mu",
##         "H0: mu <= 10, alpha=.05, n=64, mu_c=10.411", sep="\n"),
       par.settings = list(clip = list(panel = "off")),
       xlab=list("H0 value of mu", cex=1.5),
       ylab=list("P(retain H0)", cex=1.5),
       scales=list(cex=1.5),
       type="l",
       panel=function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(h=.95, v=10)
         panel.abline(h=.5, v=10.411, lty=3)
         panel.abline(h=c(0,1), lty=2)
         panel.axis("bottom", at=c(10,10.411), labels=c("mu0","mu_c"),
                    half=FALSE,
                    tck = 4,
                    outside=TRUE,
                    tick=TRUE, text.cex=1, rot=0)
         panel.axis("left", at=c(0.5,0.95),
                    half=FALSE,
                    ## line=1.2,
                    outside=TRUE,
                    tick=TRUE, text.cex=1.45)
       })
)

### 2

print(
xyplot((1 - pnorm(41.644 - 4*mu)) ~ mu, data=oc.data,
##       main=paste("2. Power Curve:  1-beta(mu) =",
##         "Probability of Rejecting H0 for specified value of mu",
##         "Ho: mu <= 10, alpha=.05, n=64, mu_c=10.411", sep="\n"),
       par.settings = list(clip = list(panel = "off")),
       xlab=list("H0 value of mu", cex=1.5),
       ylab=list("P(reject H0)", cex=1.5),
       scales=list(cex=1.5),
       type="l",
       panel=function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(h=.05, v=10)
         panel.abline(h=.5, v=10.411, lty=3)
         panel.abline(h=c(0,1), lty=2)
         panel.axis("bottom", at=c(10,10.411), labels=c("mu0","mu_c"),
                    half=FALSE,
                    tck = 4,
                    outside=TRUE,
                    tick=TRUE, text.cex=1, rot=0)
         panel.axis("left", at=c(0.5,0.95),
                    half=FALSE,
                    ## line=1.2,
                    outside=TRUE,
                    tick=TRUE, text.cex=1.45)
       })
)
}, s={
### 1

print(
xyplot(pnorm(41.644 - 4*mu) ~ mu, data=oc.data,
##       main=paste("1. Operating Characteristics Curve:  beta(mu) =",
##         "Probability of Retaining H0 for specified value of mu",
##         "H0: mu <= 10, alpha=.05, n=64, mu_c=10.411", sep="\n"),
       xlab=list("H0 value of mu", cex=2.5),
       ylab=list("P(retain H0)", cex=2.5),
       scales=list(cex=1.5),
       type="l",
       panel=function(x, y, ...) {
         panel.xyplot(x, y, ...)
         abline(h=.95, v=10)
         abline(h=.5, v=10.411, lty=3)
         abline(h=c(0,1), lty=2)
         axis(1, at=c(10, 10.411), label=c("mu0","mu_c"),
              line=1.7, ticks=FALSE, cex=1.5)
         axis(1, at=c(10, 10.411), labels=FALSE, tck=-.05)
         axis(2, at=c(0.95, .5), line=1.2, ticks=FALSE, cex=1.45)
       })
)

### 2
print(
xyplot((1 - pnorm(41.644 - 4*mu)) ~ mu, data=oc.data,
##       main=paste("2. Power Curve:  1-beta(mu) =",
##         "Probability of Rejecting H0 for specified value of mu",
##         "Ho: mu <= 10, alpha=.05, n=64, mu_c=10.411", sep="\n"),
       xlab=list("H0 value of mu", cex=2.5),
       ylab=list("P(reject H0)", cex=2.5),
       scales=list(cex=1.5),
       type="l",
       panel=function(x, y, ...) {
         panel.xyplot(x, y, ...)
         abline(h=.05, v=10)
         abline(h=.5, v=10.411, lty=3)
         abline(h=c(0,1), lty=2)
         axis(1, at=c(10,10.411), label=c("mu0","mu_c"),
              line=1.7, ticks=FALSE, cex=1.5)
         axis(1, at=c(10, 10.411), labels=FALSE, tck=-.05)
         axis(2, at=c(0.95, .5), line=1.2, ticks=FALSE, cex=1.45)
       })
)}
)
