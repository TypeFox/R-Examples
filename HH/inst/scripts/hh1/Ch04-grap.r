library(HH)

#### grap/code/ecolo.s
## Ecological correlation

eco.corr <- data.frame(x=1:99,
                       g=factor(rep(1:3, c(33,33,33))),
                       eps=rnorm(99))
eco.corr$y <- 30*as.numeric(eco.corr$g) - .3*eco.corr$x + 8*eco.corr$eps

cor(eco.corr$x, eco.corr$y)

for (i in 1:3) print(cor(eco.corr$x[eco.corr$g==i], eco.corr$y[eco.corr$g==i]))
                     
xyplot(y ~ x, group=g, data=eco.corr, panel=panel.superpose)

summary(lm(y ~ x, data=eco.corr), corr=FALSE)

summary(lm(y ~ g + x, data=eco.corr), corr=FALSE)


## parallel lines
ecop.aov <- ancova(y ~ x + g, data=eco.corr,
                   par.strip.text=list(cex=1.2),
                   ignore.groups=TRUE,
                   layout=c(3,2), between=list(y=2),
                   main=list(cex=1.5)
                   )
ecop.aov


## top two panels on their own page
ecop.aov <- ancova(y ~ x + g, data=eco.corr,
                   par.strip.text=list(cex=1.2),
                   ignore.groups=TRUE,
                   layout=c(2,1), skip=c(FALSE,FALSE,FALSE,TRUE,FALSE,FALSE),
                   aspect=1.2,
                   main=list(cex=1.5, main="")
                   )
ecop.aov
## export.eps(hh("grap/figure/grap.ecop.eps"))  ## third page of three


#### grap/code/njgolf-read.s
data(njgolf)

#### grap/code/njgolf.graphs.s
## njgolf was defined in njgolf-read.s

print(position=c(0,.4, 1,1), more=FALSE,
      xyplot(sprice ~ lotsize, data=njgolf,
             scale=list(y=list(alternating=2)),
             ylab="selling price")
)
## export.eps(hh("grap/figure/grap.pric.lot.eps"))

print(position=c(0,.35, 1,1), more=FALSE,
      xyplot(sprice ~ lotsize,
             group=factor(njgolf$lotsizef=="house",
               labels=c("condominium","house")),
             data=njgolf,
             panel=panel.superpose,
             pch=16,
             scale=list(y=list(alternating=2)),
             ylab="selling price",
             key=list(
               border=TRUE,
               text=list(c("Condominium","House")),
               points=list(
                 col=trellis.par.get("superpose.symbol")$col[1:2],
                 pch=16),
               space="top"))
)
## export.eps(hh("grap/figure/grap.pric.lot.color.eps"))


xyplot(sprice ~ beds, data=njgolf)

xyplot(sprice ~ drarea, data=njgolf)

xyplot(sprice ~ kitarea, data=njgolf)

xysplom(sprice ~ beds + drarea + kitarea, data=njgolf)

print(position=c(0,.4, 1,1), more=FALSE,
xysplom(sprice ~ beds + drarea + kitarea, data=njgolf,
        ## above is required, below is for prettier formatting
        scales=list(x=list(relation="free"), cex=1),
        par.strip.text=list(cex=1.4),
        xlab=list("measures of dwelling size", cex=1.2), ylab="",
        between=list(x=1.5))
)
## export.eps(hh("grap/figure/grap.pric.bdk.eps"))


xysplom(sprice ~ beds + drarea + kitarea |
        factor(lotsizef=="house", labels=c("condominium","house")),
        data=njgolf,
        ## above is required, below is for prettier formatting
        scales=list(x=list(relation="free"), cex=1, y=list(alternating=FALSE)),
        par.strip.text=list(cex=1.4),
        xlab=list("measures of dwelling size", cex=1.2), ylab="",
        between=list(x=1.5), layout=c(3,2))


if.R(r={
  bdk <- with(njgolf,
              rbind(data.frame(sprice=sprice, value=beds,    var="beds",    lotsizef=lotsizef),
                    data.frame(sprice=sprice, value=drarea,  var="drarea",  lotsizef=lotsizef),
                    data.frame(sprice=sprice, value=kitarea, var="kitarea", lotsizef=lotsizef)))
  tmp <-
    xyplot(sprice ~ value | var * factor(lotsizef=="house", labels=c("condominium","house")),
           data=bdk,
           scales=list(
             cex=1,
             x=list(alternating=FALSE, relation="free"),
             y=list(alternating=FALSE, relation="same")
             ),
           par.strip.text=list(cex=1.2),
           xlab=list("measures of dwelling size", cex=1.2),
           between=list(x=1.5, y=1))
  combineLimits(tmp)
}, s=
     xysplom(sprice ~ beds + drarea + kitarea |
             factor(lotsizef=="house", labels=c("condominium","house")),
             data=njgolf,
             ## above is required, below is for prettier formatting
             panel="panel.cartesian",
             axis3.line=3.6,
             x.label="",
             y.label="",
             group.label.side="",
             scales=list(
               cex=1,
               x=list(alternating=TRUE, labels=FALSE, ticks=FALSE),
               y=list(alternating=TRUE, labels=FALSE, ticks=FALSE, relation="same")
               ),
             par.strip.text=list(cex=1.2),
             xlab=list("measures of dwelling size", cex=1.2), ylab=list(label=""),
             between=list(x=1.5,y=4), layout=c(3,2))
     )
## export.eps(hh("grap/figure/grap.pric.bdkc.eps"))


tmp <- cbind(njgolf[,c("sprice","lotsize","beds","drarea","kitarea")],
             cond.house=factor(njgolf$lotsizef=="house",
               labels=c("condominium","house")))

if.R(s=splom(~ tmp, cex=.6),
     r=splom(~ tmp, axis.text.cex=.5, varname.cex=.8, xlab=""))
## export.eps(hh("grap/figure/grap.pbdkcl.eps"))

if.R(s=
     splom(~ tmp, cex=.6, pch=16, group=tmp$cond.house, panel=panel.superpose,
           key=list(
             border=TRUE,
             text=list(c("Condominium","House")),
             points=list(
               col=trellis.par.get("superpose.symbol")$col[1:2],
               pch=16),
             space="right")),
     r=
     splom(~ tmp, axis.text.cex=.5, varname.cex=.8, xlab="",
           pch=16, group=tmp$cond.house, panel=panel.superpose,
           key=list(
             border=TRUE,
             text=list(c("Condominium","House")),
             points=list(
               col=trellis.par.get("superpose.symbol")$col[1:2],
               pch=16),
             space="right"))
     )
## export.eps(hh("grap/figure/grap.pbdkcl.color.eps"))

if.R(s=splom(~ tmp[,1:5] | tmp[,6], par.strip.text=list(cex=1.2), cex=.6),
     r=splom(~ tmp[,1:5] | tmp[,6], par.strip.text=list(cex=1.2),
             axis.text.cex=.5, varname.cex=.8, xlab=""))
## export.eps(hh("grap/figure/grap.pbdkc-l.eps"))


#### grap/code/grap.read.le.r
#### grap/code/grap.read.le.s
data(tv)

#### grap/code/grap.f1.s
par(mfrow=c(1,1))
par(pty="s")
plot(male.life.exp ~ fem.life.exp, data=tv,
     main="life expectancy",
     pch=16,
     xlim=c(50,85), ylim=c(50,85))
abline(a=0, b=1)
text(x=tv["Japan","fem.life.exp"],
     y=tv["Japan","male.life.exp"],
     "(82,76)", adj=0)
## export.eps(hh("grap/figure/grap.f1.eps"))


#### grap/code/grap.f2.s
par(mar=par("mar")+c(0,1,0,0))

split.screen(matrix(nrow=5, ncol=4, byrow=TRUE,
                    c(0.0, 0.5, 0.4, 1.0,
                      0.0, 0.5, 0.0, 0.4,
                      0.5, 1.0, 0.7, 1.0,
                      0.5, 1.0, 0.4, 0.7,
                      0.5, 1.0, 0.0, 0.35)))


screen(1)
par(pty="s")
plot(male.life.exp ~ fem.life.exp, data=tv,
     type="n",
     xlim=c(50,85), ylim=c(50,85),
     cex=1)
title("a. abbreviated names", cex=.8)
text(y=tv$male.life.exp, x=tv$fem.life.exp,
     abbreviate(row.names(tv)), cex=.8)
abline(a=0, b=1)

screen(2)
par(pty="s")
plot(male.life.exp ~ fem.life.exp, data=tv,
     type="n",
     xlim=c(50,85), ylim=c(50,85),
     cex=1)
points(y=tv$male.life.exp, x=tv$fem.life.exp, pch=16, cex=.6)
title("b. simulated interactive", cex=.8)
abline(a=0, b=1)
text(y=tv$male.life.exp[2], x=tv$fem.life.exp[2],
     row.names(tv)[2], cex=.8, adj=1)

screen(3)
par(pty="s")
plot(male.life.exp ~ fem.life.exp, data=tv,
     type="n",
     pch=16, cex=1)
title("c. square, unequal scale", cex=.8)
points(y=tv$male.life.exp, x=tv$fem.life.exp, pch=16, cex=.6)

screen(4)
par(pty="s")
plot(male.life.exp ~ fem.life.exp, data=tv,
     type="n",
     pch=16, cex=1)
title("d. x=y line, square, unequal scale", cex=.8)
points(y=tv$male.life.exp, x=tv$fem.life.exp, pch=16, cex=.6)
abline(a=0, b=1)

screen(5)
par(pty="m")
plot(male.life.exp ~ fem.life.exp, data=tv,
     type="n",
     pch=16, cex=1,
     xlim=c(50,85), ylim=c(50,85))
points(y=tv$male.life.exp, x=tv$fem.life.exp, pch=16, cex=.6)
mtext("e. same xlim and ylim, unequal scale,\nleast squares line",
      line=4, cex=1.2)
abline(a=0, b=1)
abline(lm(male.life.exp ~ fem.life.exp, data=tv))

close.screen(all=TRUE)
## export.eps(hh("grap/figure/grap.f2.eps"))


#### grap/code/grap.identify.s
par(pty="s")
plot(male.life.exp ~ fem.life.exp, data=tv,
     type="n",
     xlim=c(50,85), ylim=c(50,85),
     cex=1)
points(y=tv$male.life.exp, x=tv$fem.life.exp, pch=16, cex=.6)
title("Interactive, similar to Figure 4.8.b.", cex=.8)
abline(a=0, b=1)

if (FALSE) { ## identify() doesn't work inside source()
             ## For manual use, enter these lines directly in the console.
if.R(r={old.par <- par(xpd=TRUE)},
     s={})
identify(y=tv$male.life.exp, x=tv$fem.life.exp, labels=row.names(tv))
if.R(r=par(old.par),
     s={})
## Put the cursor on any point, then left-click.
## Repeat until satisfied, then right-click.
}

#### grap/code/grap.f3.s
#### grap/code/grap.f3.s
if.R(r={
splom( ~ tv[,c(4,5,1,2,3)],
      main=list("Televisions, Physicians, and Life Expectancy", cex=1.4),
      axis.text.cex=.5,
      cex=.7, varname.cex=1)
## export.eps(hh("grap/figure/grap.f3.eps"))
}, s={
  splom( ~ tv[,c(4,5,1,2,3)],
        main=list("Televisions, Physicians, and Life Expectancy", cex=1.4),
        superpanel=panel.pairs.hh, subpanel.scales=list(cex=.5),
        cex=.7, panel.cex=1)
})
## export.eps(hh("grap/figure/grap.f3.eps"))


#### grap/code/grap.f11.s,
par(pty="m")
old.par <- par(cex=1.7, pch=16, oma=par("mar")+c(0,4,0,0))
pairs(tv, cex=1)
mtext(outer=TRUE, 'pairs with NW--SE diagonal and rectangular panels',
      cex=1.5, line=1.5)

par(old.par)
par(pty="m")
## export.eps(hh("grap/figure/grap.f11.eps"))


#### grap/code/grap.f12.s
f12b <- matrix(c("Var1",1,2,4,
                 "1'","Var2",3,5,
                 "2'","3'","Var3",6,
                 "4'","5'","6'","Var4"), 4, 4)
dd <- c(1,6,11,16)

par(mfrow=c(1,2))
old.par <- par(pty="s", usr=c(.5, 4.5,  .5, 4.5), mar=par("mar")+c(0,1,0,0))

plot(x=c(.5, 4.5), y=c(.5, 4.5), type="n", las=1, yaxt="n",
     cex=1.2,
     main="a. multiple axes of symmetry.",
     xlab="variables", ylab="variables")
axis(2, at=1:4, labels=rep("",4), cex=1.2)
if.R(s=
     axis(2, at=1:4, labels=4:1, cex=1.2, line=.5, ticks=FALSE)
     ,r=
     axis(2, at=1:4, labels=4:1, cex=1.2, las=1)
)
text(f12b[-dd], y=5-row(f12b)[-dd], x=col(f12b)[-dd], cex=2)
text(f12b[ dd], y=5-row(f12b)[ dd], x=col(f12b)[ dd], cex=1.4)
abline(a=5, b=-1)
abline(a=-2, b=1, lty=3)
abline(a=-1, b=1, lty=2)
abline(a= 0, b=1, lty=2)
abline(a= 1, b=1, lty=2)
abline(a= 2, b=1, lty=2)
if.R(s={
  arrows(2.8, .75, 2.8,1.14, open=TRUE, size=.3, rel=TRUE) ## up arrow
  arrows(3.8,1.7, 4.15,1.7, open=TRUE, size=.3, rel=TRUE)  ## right arrow
},r={
  arrows(2.8, .75, 2.8,1.14, length=.1) ## up arrow
  arrows(3.8,1.7, 4.15,1.7, length=.1)  ## right arrow
})


plot(x=c(.5,4.5), y=c(.5,4.5), type="n", las=1,
     cex=1.2,
     main="b. single axis of symmetry.",
     xlab="variables", ylab="variables")
text(f12b[-dd], y=row(f12b)[-dd], x=col(f12b)[-dd], cex=2)
text(f12b[ dd], y=row(f12b)[ dd], x=col(f12b)[ dd], cex=1.4)
abline(a=0,b=1)
if.R(s={
  arrows(2.8,3.75, 2.8,4.14, open=TRUE, size=.3, rel=TRUE) ## up arrow
  arrows(3.8,2.7, 4.15,2.7, open=TRUE, size=.3, rel=TRUE)  ## right arrow
},r={
  arrows(2.8,3.75, 2.8,4.14, length=.1) ## up arrow
  arrows(3.8,2.7, 4.15,2.7, length=.1)  ## right arrow
})

par(old.par)

## export.eps(hh("grap/figure/grap.f12.eps"))


#### grap/code/grap.f5.s
if.R(r={
  splom( ~ tv[, 1:3],
        main=list("Televisions, Physicians, and Life Expectancy", cex=1.4),
        axis.text.cex=.7,
        cex=.9, varname.cex=1.3)
},s={
  splom( ~ tv[, 1:3],
        main=list("Televisions, Physicians, and Life Expectancy", cex=1.4),
        superpanel=panel.pairs.hh, subpanel.scales=list(cex=.7),
        cex=.9, panel.cex=1.3)
})
## export.eps(hh("grap/figure/grap.f5.eps"))


#### grap/code/grap.f6.s
if.R(r={splom( ~ cbind(tv[,1,drop=FALSE], log(tv[, 2:3])),
              main=list("log(Televisions, Physicians), and Life Expectancy", cex=1.4),
              axis.text.cex=1.1,
              cex=.9, varname.cex=1.3)
      },
     s={
       splom( ~ cbind(tv[,1,drop=FALSE], log(tv[, 2:3])),
             main=list("log(Televisions, Physicians), and Life Expectancy", cex=1.4),
             superpanel=panel.pairs.hh, subpanel.scales=list(cex=1.1),
             cex=.9, panel.cex=1.3)
     })
## export.eps(hh("grap/figure/grap.f6.eps"))


#### grap/code/grap.f8.s
## variants of the power transformations
## grap.f8.s
par(mfrow=c(1,1))
frame()
par(cex=.8)
dy2 <- format(c(-1, -0.5, 0, 0.5, 1, 2))

x <- seq(0, 2, length=101)
y <- data.matrix(ladder.f(x+.001))
if.R(s=
     par(fig=c(-.05,.31,.4,1))  ## left 1/3 and top .6 of the graph sheet
     ,r=
     par(fig=c( .00,.34,.4,1)      )  ## left 1/3 and top .6 of the graph sheet
)
matplot(x=x, y=y, type="l", lty=1, col=1, xlim=c(0,2), ylim=c(-2,4), err=-1,
        main="a. simple powers\nnegative reciprocals (correct)")
text(x[101]+.1, y[101,]+c(.2,-.1,0,0,0,-.5), dy2, adj=0)

y[,1:2] <- -y[,1:2]
if.R(s=
par(fig=c(.30,.66,.4,1))  ## middle 1/3 and top .6 of the graph sheet
     ,r=
par(fig=c(.33,.67,.4,1),new=TRUE)  ## middle 1/3 and top .6 of the graph sheet
)
matplot(x=x, y=y, type="l", lty=1, col=1, xlim=c(0,2), ylim=c(-2,4), err=-1,
        main="b. simple powers\n positive reciprocals (incorrect)")
text(x[101]+.1, y[101,]+c(-.2,.3,0,0,0,-.5), dy2, adj=0)

yy <- data.matrix(ladder.fstar(x+.001))
if.R(s=
par(fig=c(.65,1.01,.4,1))  ## right 1/3 and top .6 of the graph sheet
     ,r=
par(fig=c(.66,1.00,.4,1),new=TRUE)  ## right 1/3 and top .6 of the graph sheet
)
matplot(x=x, y=yy, type="l", lty=1, col=1, xlim=c(0,2), ylim=c(-2,2), err=-1,
        main="c. scaled powers")
text(x[101]+.1, yy[101,]+c(-.21,-.12,0,0.1,0.3,0.1), dy2, adj=0)
## export.eps(hh("grap/figure/grap.f8.eps"))


#### The HH text references filename splus.library/code/ladder.cartesian.s
#### which is a typo for file
#### splus.library/code/ladder.s
#### and means the ladder function in the HH package.


#### grap/code/grap.f7.s
## Default size labels
ladder(life.exp ~ ppl.per.phys, data=tv,
       main="Ladder of Powers for Life Expectancy and People per Physician")

## Larger labels, but smaller plotting area
ladder(life.exp ~ ppl.per.phys, data=tv,
       main="Ladder of Powers for Life Expectancy and People per Physician",
       par.strip.text=list(cex=1),
       axis3.line=if.R(r=.7, s=1),
       dsx="ppp", dsy="le",
       scales.in=list(alternating=TRUE, labels=FALSE, ticks=FALSE, cex=1), g.cex=1)

## Larger labels and larger plotting area.
## this opens the postscript device, so it will not happen automatically
if (FALSE) {
if.R(r={
  postscript("Fig416.eps", paper="special", width=11, height=11)
  ladder(life.exp ~ ppl.per.phys, data=tv,
         main="Ladder of Powers for Life Expectancy and People per Physician",
         par.strip.text=list(cex=1), axis3.line=.8,
         scales.in=list(alternating=TRUE, labels=FALSE, ticks=FALSE, cex=1), g.cex=1)
  dev.off()  
},
     s={
       ## Half the graph is off-screen.  All of the graph is on the exported eps file.
       print(position=c(-.1,0, 1.2, 1.9),
             ladder(life.exp ~ ppl.per.phys, data=tv,
                    main="Ladder of Powers for Life Expectancy and People per Physician",
                    par.strip.text=list(cex=1), axis3.line=2.0, dsx="ppp", dsy="le",
                    scales.in=list(alternating=TRUE, labels=FALSE, ticks=FALSE, cex=1), g.cex=1)
             )
     })
}
## export.eps(hh("grap/figure/grap.f7.eps"))


#### grap/code/grap.f9.s
par(mfrow=c(1,2))
frame()
par(cex=.8)

x <- sort(tv[,"life.exp"])

y <- data.matrix(ladder.fstar(x))
par(fig=c(0,.5,.4,1))  ## left 1/2 and top .6 of the graph sheet
matplot(x=x, y=y, type="l", lty=1, col=1,
        main="a. Powers of Life Expectancy")

y <- y*matrix(c(10, 6, 2.5, 1/1.4, 10/60, 1/200), byrow=TRUE, nrow=40, ncol=6)
par(fig=c(.5,1,.4,1), new=TRUE)  ## right 1/2 and top .6 of the graph sheet
matplot(x=x, y=y, type="l", lty=1, col=1,
        main="b. Powers of Life Expectancy\n matched by scaling")
par(mfrow=c(1,1))

## export.eps(hh("grap/figure/grap.f9.eps"))


#### grap/code/grap.f10.r
#### grap/code/grap.f10.s
if.R(r={
x <- sort(tv[,"ppl.per.phys"])
y <- ladder.f(x)
g <- paste("ppp ^",names(y))
g <- ordered(g, g)
g <- g[col(as.matrix(y))]
gg <- y
gg[] <- ''

print(position=c(0,.475,1,.95), more=TRUE,
bwplot(unlist(y) ~ unlist(gg) | g, layout=c(6,1), xlab="", ylab="",
       horizontal=FALSE,
       main=list("a. boxplots", cex=1),
       par.strip.text=list(cex=1),
       scales=list(relation="free", cex=.8))
)

print(position=c(0,0,1,.475), more=TRUE,
stripplot(unlist(y) ~ unlist(gg) | g, layout=c(6,1), xlab="", ylab="",
        horizontal=FALSE,
        main=list("b. stripplots", cex=1),
        par.strip.text=list(cex=1),
        scales=list(relation="free", cex=.8))
)

print(position=c(0,.95,1,1),
      xyplot(1 ~ .5, main='Powers of "People per Physician"',
             xlab="", ylab="",
             type="n", scales=list(draw=FALSE),
             par.settings=list(axis.line=list(col="transparent"))))
##
## this main title is not right yet
## mtext('Powers of "People per Physician"', outer=TRUE, cex=1.2)
##
## export.eps(hh("grap/figure/grap.f10.eps"))

stem(x, scale=2)
stem(x)

for (i in names(y)) {
    cat("\n(People per Physician)^", i,'\n')
    stem(y[[i]])
}
}, s={
x <- sort(tv[,"ppl.per.phys"])
y <- ladder.f(x)
g <- paste("ppp ^",names(y))
g <- ordered(g, g)
g <- g[col(y)]
gg <- y
gg[] <- ''

print.trellis(split=c(1,2, 1,2), more=TRUE,
t(bwplot(unlist(gg) ~ unlist(y) | g, layout=c(6,1), xlab='',
       main=list("a. boxplots", cex=1),
       par.strip.text=list(cex=1),
       scales=list(x=list(relation="free"), cex=.8)))
)

print.trellis(split=c(1,1, 1,2), more=FALSE,
t(stripplot(unlist(gg) ~ unlist(y) | g, layout=c(6,1), xlab='',
        main=list("b. stripplots", cex=1),
        par.strip.text=list(cex=1),
        scales=list(x=list(relation="free"), cex=.8)))
)

mtext('Powers of "People per Physician"', outer=TRUE, cex=1.2)
## export.eps(hh("grap/figure/grap.f10.eps"))

stem(x)
stem(x, fence=12)
stem(x, fence=12, nl=2)

for (i in names(y)) {
    cat("\n(People per Physician)^", i,'\n')
    stem(y[[i]])
}
})

#### grap/code/draft70mn-read.r
data(draft70mn)
