##   c:/HOME/rmh/hh/regc/:
##   \cref{rent4.s}
##   \cref{rent4b.s}
##   \cref{norm.prob.plot.s}
##   \cref[splus.library]{lm.case.s}
##   \cref[splus.library]{plotcasediagtrellis.s}
##   \cref{dfbeta.s}
## ## A simple executable version of this algorithm is in file \cref{dfbeta.s}.
## ## A more complete version (with protection against near-singularity)
## ## is in the \Splus\ function {\tt lm.influence}.
## -rwx------+   1 rmh None    51522 2004-05-23 16:08 regc.tex


## rent4.s
## Sanford Weisberg, Applied Linear Regression, 2nd Edition, Wiley (1985).
## Problem 6.5, page 163 with data in Table 6.10, page 162.

## S-Plus/R code by Heiberger and Holland

data(rent)

rent[1:3,]

if.R(s=
splom( ~ rent[,c(1,2,6,3,4)] | rent$lime, pch=16, cex=.6,
      superpanel=panel.pairs.hh, panel.cex=.8,
      subpanel.scales=list(cex=.5), pscales=4,
      par.strip.text=list(cex=1.4))
     ,r=
splom( ~ rent[,c(1,2,6,3,4)] | rent$lime, pch=16, cex=.6,
      superpanel=panel.pairs, varname.cex=.8,
      axis.text.cex=.5, pscales=4,
      par.strip.text=list(cex=1.4))
     )
##    main="scatterplot matrices | lime"
## export.eps(hh("regc/figure/rent1.eps"))

rent.lm3l <- lm(rnt.alf ~ rnt.till + cow.dens + prop.past + lime,
                data=rent)
summary(rent.lm3l, corr=FALSE)
anova(rent.lm3l)
## rent.lm31.st


## proportion in pasture is not needed
if.R(s=
splom(~rent[,c(1,2,3)] | rent$lime, pch=16,
      superpanel=panel.pairs.hh, panel.cex=1.2,
      subpanel.scales=list(cex=.7),
      par.strip.text=list(cex=1.4))
     ,r=
splom(~rent[,c(1,2,3)] | rent$lime, pch=16,
      superpanel=panel.pairs, varname.cex=1.2,
      axis.text.cex=.7,
      par.strip.text=list(cex=1.4))
)
##    main="scatterplot matrices"
## export.eps(hh("regc/figure/rent2.eps"))

rent.lm4ln <- lm(rnt.alf ~ rnt.till + cow.dens +
                 lime + cow.dens:lime, data=rent)
summary(rent.lm4ln, corr=FALSE)
anova(rent.lm4ln)
## rent.lm4ln.st


## residuals against x variables conditioned on lime
if.R(s=
print(position=c(0,0,1,.9),
xysplom(resid(rent.lm4ln) ~ rnt.till + cow.dens | lime, data=rent,
        layout=c(2,2),
        xlab="", ylab="",
        x.label="", y.label="",
        group.label.side="",
        par.strip.text=list(cex=1.2),
        panel=panel.cartesian,
        axis3.line=3.6,
        scales=list(
          relation="same",
          alternating=FALSE, labels=FALSE, ticks=FALSE),
        between=list(x=1, y=3))
)
,r=
print(position=c(0,0,1,.9),
xysplom(resid(rent.lm4ln) ~ rnt.till + cow.dens | lime, data=rent,
        layout=c(2,2),
        xlab="", ylab="",
        x.label="", y.label="",
        group.label.side="",
        par.strip.text=list(cex=1.2),
##        panel=panel.cartesian,
        axis3.line=3.6,
        scales=list(
          relation="same",
          alternating=FALSE, labels=FALSE, ticks=FALSE),
        between=list(x=1, y=3))
)
)
## export.eps(hh("regc/figure/rent4lnres.eps"))




## alf.till ratio
rent.lm12p <- lm(alf.till ~ lime * cow.dens + prop.past, data=rent)
summary(rent.lm12p, corr=FALSE)
anova(rent.lm12p)
## rentlm12p.st

rent.lm12m <- ancova(alf.till ~ lime * cow.dens, data=rent,
                     par.strip.text=list(cex=1.2))
print(position=c(0,0, 1,.6),
      attr(rent.lm12m,"trellis"))
## export.eps(hh("regc/figure/rent.lm12m.eps"))
anova(rent.lm12m)
## rentlm12.st
summary.lm(rent.lm12m)  ## coefficients illustrated in Figure 11.11


old.par <- par(pch=16, mar=c(5,6,4,2)+.1)
par(mfrow=c(2,3))
plot(rent.lm12m, which=1:6, cex=.8)
par(mfrow=c(1,1))
par(old.par)
## export.eps(hh("regc/figure/rent.plot.lm12m.eps"))


rent.case12m <-
  if.R(r=case(rent.lm12m),
       s=case.lm(rent.lm12m))
if.R(s=
print(position=c(0,0,1,1.1),  ## oversize for printing in book
      plot(rent.case12m, rent.lm12m, par.strip.text=list(cex=.9))
)
,r=
     plot(rent.case12m, rent.lm12m, par.strip.text=list(cex=.9),
          layout=c(2,3))
)
.lm.case.large ## needed to resolve overstriking in the graph
## export.eps(hh("regc/figure/rent.diag.lm12m.eps"))

## locate the outliers (19,33,60, 49) identified by the diagnostics
print(position=c(0,0,1,.7),
xyplot(alf.till ~ cow.dens | lime, data=rent, RowNames=row.names(rent),
       panel=function(x, y, subscripts, RowNames, ...) {
         panel.xyplot(x, y, ...)
         subs <- match(c(19,33,60, 49), subscripts, 0)
         if.R(r=panel.text, s=text)(x[subs], y[subs],
                              RowNames[subscripts][subs], adj=0, cex=1.5)
       },
       par.strip.text=list(cex=1.4),
       between=list(x=1),
       scales=list(alternating=FALSE, cex=1.2),
       xlab=list(cex=1.4), ylab=list(cex=1.4),
       pch=16, cex=.8,
       main=" ")
      )
## export.eps(hh("regc/figure/rent.text.lm12m.eps"))


## remove four observations
rent.lm12ms <- ancova(alf.till ~ lime * cow.dens,
                      data=rent[-c(19, 33, 60, 49),],
                      ylim=range(rent$alf.till),
                      par.strip.text=list(cex=1.2))
print(position=c(0,0, 1,.6),
      attr(rent.lm12ms,"trellis"))
## export.eps(hh("regc/figure/rent.lm12ms.eps"))
anova(rent.lm12ms)
## rent.lm12ms.st


## rent4b.s
## continue after rent4.s

## normal plot
qqmath( ~ resid(rent.lm12m))

## normal plot with straight line, identified outliers, ylim control
qqmath( ~ resid(rent.lm12m),
       ylim=c(-.4,.88),
       xlab=list(cex=1.4), ylab=list(cex=1.4),
       scales=list(cex=1.4),
       panel=function(x, y,...) {
         if.R(r={y <- x
                 x <- qnorm(ppoints(y))[order(order(y))]},
              s={})
         panel.qqmathline(y=y, distribution=qnorm)
         if.R(r=panel.qqmath(x=y, ...),
              s=panel.qqmath(x=x, y=y, ...))
         y.sort <- sort(y)
         y.len <- length(y)
         ii.label <- names(y.sort[c(1,65:67)])
         ii <- match(ii.label, names(y))
         if.R(s=text(x=x[ii], y=y[ii], names(y)[ii], adj=0, cex=1.4),
              r=panel.text(x=x[ii], y=y[ii], names(y)[ii], adj=0, cex=1.4))
       }
       )
## export.eps(hh("regc/figure/rent.residn.eps"))
## this graph looks not normal
## also do the test
shapiro.test(resid(rent.lm12m))



## normal plot after excluding four points identified from the
## diagnostics plots
qqmath( ~ resid(rent.lm12ms),
       panel=function(x,y,...) {
         if.R(r={y <- x
                 x <- qnorm(ppoints(y))[order(order(y))]},
              s={})
         panel.qqmathline(y=y, distribution=qnorm)
         if.R(r=panel.qqmath(x=y, ...),
              s=panel.qqmath(x=x, y=y, ...))
       }
       )
## this graph looks normal
## also do the test
shapiro.test(resid(rent.lm12ms))




rent.X <- update(rent.lm12m, x=TRUE)$x[,2:4] ## model has interactions
residual.plots(rent.lm12m, X=rent.X)

tmp <- residual.plots(rent.lm12m, X=rent.X)

## display on screen
if.R(s={
  print(position=c(-.025,-.050, 1.025,  .250), more=TRUE, tmp[[4]])
  print(position=c(-.025, .200, 1.025,  .500), more=TRUE, tmp[[3]])
  print(position=c(-.025, .450, 1.025,  .750), more=TRUE, tmp[[2]])
  print(position=c(-.025, .700, 1.025, 1.000), more=FALSE, tmp[[1]])
},r={
  ## 4*b - 3*a == 1
  a <- .05
  b <- (1+3*a)/4
  print(position=c(0,    .000, 1.000,   b    ), more=TRUE,  update(tmp[[4]], main=list(tmp[[4]]$main, cex=.9)))
  print(position=c(0,   b-  a, 1.000, 2*b-  a), more=TRUE,  update(tmp[[3]], main=list(tmp[[3]]$main, cex=.9)))
  print(position=c(0, 2*b-2*a, 1.000, 3*b-2*a), more=TRUE,  update(tmp[[2]], main=list(tmp[[2]]$main, cex=.9)))
  print(position=c(0, 3*b-3*a, 1.000,   1.000), more=FALSE, update(tmp[[1]], main=list(tmp[[1]]$main, cex=.9)))
}
)

## detailed look at each of the panels of the plot.case
tmp <- 
if.R(s=
     plot(rent.case12m, rent.lm12m, par.strip.text=list(cex=.9))
     ,r=
     plot(rent.case12m, rent.lm12m, par.strip.text=list(cex=.9),
          layout=c(2,3))
)
## change the printing details
tmp$layout <- c(1,1)
tmp$main <- list(tmp$main[[1]], cex=1.4)
tmp$par.strip.text <- list(cex=2)
tmp
##
rent.diag.names <- make.names(tmp$glist$group$levels)
rent.diag.names[c(1,2,4)] <- c("stu.res","si","cook")
rent.diag.names <- paste(rent.diag.names,  ".eps",sep="")
##
## export.eps(hh(paste("regc/figure", rent.diag.names[1], sep="/")))
## export.eps(hh(paste("regc/figure", rent.diag.names[2], sep="/")))
## export.eps(hh(paste("regc/figure", rent.diag.names[3], sep="/")))
## export.eps(hh(paste("regc/figure", rent.diag.names[4], sep="/")))
## export.eps(hh(paste("regc/figure", rent.diag.names[5], sep="/")))
## export.eps(hh(paste("regc/figure", rent.diag.names[6], sep="/")))
## export.eps(hh(paste("regc/figure", rent.diag.names[7], sep="/")))
## export.eps(hh(paste("regc/figure", rent.diag.names[8], sep="/")))
## export.eps(hh(paste("regc/figure", rent.diag.names[9], sep="/")))


## norm.prob.plot.s
## Calibrate your eye to reading a normal plot.
## Here are 6 normal plots fit to randomly generated normal data.
for (i in 1:3) for (j in 1:2)
  print(position=c((i-1)/3, (j-1)/2, i/3, j/2), more=TRUE,
        qqmath( ~ rnorm(67),
               ylab=list(cex=1.4),
               panel=function(x,y,...) {
                 if.R(r=panel.qqmathline(x, x, distribution=qnorm),
                      s=panel.qqmathline(y, distribution=qnorm))
                 panel.qqmath(x,y,...)}))
if.R(r=print(position=c(0, .95, 1, 1), more=FALSE,
       xyplot(0 ~ 0, panel=function(...){},
              xlab="", ylab="", scales=list(draw=FALSE),
              main=list("six randomly generated normal plots"),
              par.settings = list(axis.line = list(col = "transparent"))))
     ,s={
       mtext("six randomly generated normal plots", outer=TRUE, cex=2, line=-2)
       par(new=FALSE)
     })
## export.eps(hh("regc/figure/norm.prob.plot.eps"))


## Here are 6 normal plots fit to randomly generated nonnormal data.
tmp.r <- list()
tmp.r[["runif(67)"]] <- runif(67)
tmp.r[["rt(67, 5)"]] <- rt(67, 5)
tmp.r[["rf(67, 4, 24)"]] <- rf(67, 4, 24)
tmp.r[["rchisq(67, 15)"]] <- rchisq(67, 15)
tmp.r[["rbinom(67, 10, .2)"]] <- rbinom(67, 10, .2)
tmp.r[["rpois(67, 8)"]] <- rpois(67, 8)

for (i in 1:3) for (j in 1:2)
  print(position=c((i-1)/3, (j-1)/2, i/3, j/2), more=TRUE,
        qqmath( ~ tmp.r[[(i-1)*2+j]],
               ylab=list(names(tmp.r)[[(i-1)*2+j]], cex=1.4),
               panel=function(x,y,...) {
                 if.R(r=panel.qqmathline(x, x, distribution=qnorm),
                      s=panel.qqmathline(y, distribution=qnorm))
                 panel.qqmath(x,y,...)}))
if.R(r=print(position=c(0, .95, 1, 1), more=FALSE,
       xyplot(0 ~ 0, panel=function(...){},
              xlab="", ylab="", scales=list(draw=FALSE),
              main=list("six randomly generated normal plots"),
              par.settings = list(axis.line = list(col = "transparent"))))
     ,s={
       mtext("six randomly generated nonnormal plots", outer=TRUE, cex=2, line=-2)
       par(new=FALSE)
     })
## export.eps(hh("regc/figure/nonnorm.prob.plot.eps"))


## splus.library/lm.case.s
## functions in the HH package


## splus.library/plotcasediagtrellis.s
## functions in the HH package


## dfbeta.s
## ## A simple executable version of this algorithm is in file \cref{dfbeta.s}.
##    This produces a matrix of coefficients.
## ## A more complete version (with protection against near-singularity)
## ## is in the \Splus\ function {\tt lm.influence} with three components in the result:
##    coefficients, sigma, hat
## ## or in the \R\ function {\tt lm.influence} with four components in the result:
##    hat, coefficients, sigma, wt.res
dfbeta <- function(x,y) {
  xy.lm <- lm(y ~ x, qr=TRUE)
  R <- qr.R(xy.lm$qr)
  Q <- qr.Q(xy.lm$qr)
  h <- hat(x)
  e <- resid(xy.lm)
  b <- coef(xy.lm)

  z <- e/(1-h)

  Qz <- Q*z

  dbeta <- backsolve(R, t(Qz))

  dfbeta <- b - dbeta
  dimnames(dfbeta) <- list(c("Intercept",dimnames(x)[[2]]), dimnames(x)[[1]])
  t(dfbeta)
}

if.R(r={
  dfbeta(as.matrix(longley[,-7]), longley[,7])
  longley.lm <- lm( Employed ~ . , data=longley)
  lm.influence(longley.lm)
}, s={
  dfbeta(longley.x, longley.y)
  longley <- data.frame(longley.x, Employed = longley.y)
  longley.lm <- lm( Employed ~ . , data=longley)
  lm.influence(longley.lm)
})
