### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/grap.tex'

###################################################
### code chunk number 1: grap.tex:9-10
###################################################
library(HH)


###################################################
### code chunk number 2: grap.tex:13-18
###################################################
## the standard lattice color 2 is difficult for people with color deficient vision
data(col3x2)
## These colors look like a 3x2 color array when run through
## the vischeck simulator to see how they look for the three most
## common color vision deficiencies: Protanope, Deuteranope, Tritanope.


###################################################
### code chunk number 3: grap.tex:66-77
###################################################
data(njgolf)
## hhpdf("grap-pric-lot.pdf", width=7, height=4)
col.lot <- likertColor(2)[1:2]
xyplot(sprice ~ lotsize, data=njgolf, groups=(lotsize==0), pch=c("+","0"), cex=c(2, 1.8),
       ylab="Selling Price", xlab="Lot Size\n", col=col.lot[2:1],
       main="Single Family Homes",
       key=list(
         border=TRUE, space="bottom", columns=2,
         points=list(pch=c("0","+"), col=col.lot, cex=1.5),
         text=list(c("Zero","Positive"), col=col.lot)))
## hhdev.off()


###################################################
### code chunk number 4: grap.tex:150-178
###################################################
eco.corr <- data.frame(x=1:99,
                       g=factor(rep(1:3, c(33,33,33))),
                       eps=rnorm(99))
eco.corr$y <- 30*as.numeric(eco.corr$g) - .3*eco.corr$x + 8*eco.corr$eps

cor(eco.corr$x, eco.corr$y)

for (i in 1:3) print(cor(eco.corr$x[eco.corr$g==i], eco.corr$y[eco.corr$g==i]))

xyplot(y ~ x, group=g, data=eco.corr, panel=panel.superpose, pch=c(15,19,17), col=col3x2)

summary(lm(y ~ x, data=eco.corr))

summary(lm(y ~ g + x, data=eco.corr))
## hhpdf("grap-ecop%03d.pdf", width=7, height=4.3, onefile=FALSE)
## Top two panels of the cbind are printed on page 4 (grap-ecop004.pdf).
## Use only top two panels in book.
## The colors and lines are scrambled on the bottom 6 panels (grap-ecop00x.pdf for x=1,2,3).
A <- ancovaplot(y ~ x + g, data=eco.corr, col=col3x2)
B <- ancovaplot(y ~ x, groups=g, data=eco.corr, col=col3x2, col.line="black")
dimnames(B)[[1]][4] <- "Ignore Groups"

update(cbind(Superpose=A, "Ignore Groups"=B),
             strip.left=FALSE, between=list(x=1),
             scales=list(alternating=FALSE),
             layout=c(2,1),
             par.settings=list(layout.heights=list(strip=1)))
## hhdev.off()


###################################################
### code chunk number 5: grap.tex:242-257
###################################################
## hhpdf("grap-pric-bdk.pdf", width=7, height=3.35)
tmp <-
xyplot(sprice ~ beds + drarea + kitarea, data=njgolf, outer=TRUE,
       scales=list(x=list(relation="free"), y=list(alternating=2)),
       layout=c(3,1), between=list(x=1),
       groups=(lotsize==0), pch=c("+","0"), cex=c(2, 1.8), col=col.lot[2:1],
       ylab=NULL, ylab.right="Selling Price",
       xlab=NULL, xlab.top=list("Measures of Dwelling Size", cex=1.2),
       key=list(title="Lot Size", cex.title=1,
         border=TRUE, space="bottom", columns=2,
         points=list(pch=c("0","+"), col=col.lot, cex=c(1.5, 1.8)),
         text=list(c("Zero","Positive"), col=col.lot)))
dimnames(tmp)[[1]] <- c("Beds","Dining Room Area","Kitchen Area")
tmp
## hhdev.off()


###################################################
### code chunk number 6: grap.tex:278-300
###################################################
## hhpdf("grap-pric-bdkc.pdf", width=7.5, height=5.5)
tmp <-
combineLimits(useOuterStrips(
  update(transpose(
    xyplot(sprice ~ beds + drarea + kitarea |
           factor(lotsizef=="house", labels=c("condominium","house")),
           data=njgolf, outer=TRUE,
           groups=(lotsize==0),
           pch=c("+","0"), cex=c(2, 1.8), col=col.lot[2:1],
           ylab=NULL, ylab.right="Selling Price",
           xlab=NULL, xlab.top=list("Measures of Dwelling Size", cex=1.2),
           layout=c(2,3),
           key=list(title="Lot Size", cex.title=1,
              border=TRUE, space="bottom", columns=2,
              points=list(pch=c("0","+"), col=col.lot, cex=c(1.5, 1.8)),
              text=list(c("Zero","Positive"), col=col.lot)))
           ), between=list(x=1.5, y=1),
           scales=list(x=list(relation="free"), y=list(alternating=2)))
))
dimnames(tmp)[[1]] <- c("Beds","Dining Room Area","Kitchen Area")
tmp
## hhdev.off()


###################################################
### code chunk number 7: grap.tex:367-390
###################################################
## hhpdf("grap-pbdkcl-color.pdf", width=7, height=7)
tmp <- cbind(njgolf[,c("sprice","lotsize","beds","drarea","kitarea")],
             cond.house=factor(njgolf$lotsizef=="house",
               labels=c("condominium","house")))
tmp.splom <-
splom(~ tmp, axis.text.cex=.5, varname.cex=.8, xlab=NULL,
      group=tmp$cond.house, panel=panel.superpose,
      pch=c("0","+"), cex=c(1.8, 2), col=col.lot[1:2],
      key=list(
        border=TRUE,
        text=list(c("Condominium","House")),
        points=list(pch=c("0","+"), col=col.lot, cex=c(1.5, 1.8)),
        space="bottom", columns=2))
rrows <- c(1,1,1,1)
ccols <- c(2,3,4,5)
cols <- c("gray40","black","black","black")
## when I isolate the contents of the layer() as a function, it doesn't work.
tmp.splom + layer((function(...) {
  if (i %in% rrows && j %in% ccols)
     panel.abline(v=current.panel.limits()$xlim,
                  h=current.panel.limits()$ylim,
                  col=cols[which(i == rrows & j == ccols)], lwd=8)})())
## hhdev.off()


###################################################
### code chunk number 8: grap.tex:414-425
###################################################
## hhpdf("grap-pbdkc-l.pdf", width=9, height=5)
splom(~ tmp[,1:5] | tmp[,6], par.strip.text=list(cex=1.5),
      axis.text.cex=.5, varname.cex=.8, xlab=NULL,
      groups=tmp$cond.house,
      pch=c("0","+"), cex=c(1.3, 1.5), col=col.lot[1:2],
      key=list(
        border=TRUE,
        text=list(c("Condominium","House")),
        points=list(pch=c("0","+"), col=col.lot, cex=c(1.2, 1.6)),
        space="bottom", columns=2))
## hhdev.off()


###################################################
### code chunk number 9: grap.tex:510-513
###################################################
## hhcode("array3way.r", '
## This three-way array of scatterplots is constructed in file HHscriptnames("logi.R").
## ')


###################################################
### code chunk number 10: grap.tex:597-608
###################################################
## hhpdf("grap-f1.pdf", width=5.5, height=5.5)
data(tv)
xyplot(male.life.exp ~ fem.life.exp, data=tv,
     main="Life Expectancy", xlab="Female", ylab="Male",
     pch=19, aspect="iso", col=likertColor(2)[2],
     xlim=c(48,90), ylim=c(48,90)) +
  layer(panel.abline(a=0, b=1)) +
  layer(panel.text(x=tv["Japan","fem.life.exp"],
                   y=tv["Japan","male.life.exp"],
                   "(82,76)", pos=4))
## hhdev.off()


###################################################
### code chunk number 11: grap.tex:619-662
###################################################
## hhpdf("grap-f2.pdf", width=8, height=8)
AA <-
  xyplot(male.life.exp ~ fem.life.exp, data=tv,
         xlim=c(48,85), ylim=c(48,85), aspect="iso",
         panel=panel.text, labels=abbreviate(row.names(tv)),
         cex=.8, col=likertColor(2)[2],
         main="a. abbreviated names") +
  layer(panel.abline(a=0, b=1))


BB <-
  xyplot(male.life.exp ~ fem.life.exp, data=tv, pch=19,
         xlim=c(48,85), ylim=c(48,85), aspect="iso", col=likertColor(2)[2],
         main="b. simulated interactive") +
  layer(panel.abline(a=0, b=1)) +
  layer(panel.xyplot(x[2], y[2], cex=1.5, pch=19)) +
  layer(panel.text(x[2]+3, y[2]+3, abbreviate(row.names(tv)[2]), cex=1, pos=2))

CC <-
xyplot(male.life.exp ~ fem.life.exp, data=tv,
     pch=19, aspect=1, col=likertColor(2)[2],
     main="c. square, unequal scale")

DD <-
  xyplot(male.life.exp ~ fem.life.exp, data=tv,
       pch=19, aspect=1, col=likertColor(2)[2],
       main="d. x=y line, square, unequal scale") +
  layer(panel.abline(a=0, b=1))

EE <-
  xyplot(male.life.exp ~ fem.life.exp, data=tv,
         pch=19, aspect=.4, col=likertColor(2)[2],
         xlim=c(50,85), ylim=c(50,85),
         main="e. same xlim and ylim, unequal scale,\nleast squares line") +
  layer(panel.abline(a=0, b=1)) +
  layer(panel.abline(lm(male.life.exp ~ fem.life.exp, data=tv)))

print(position=c(0.0, 0.45, 0.55, 1.00), more=TRUE, AA)
print(position=c(0.0, 0.00, 0.55, 0.45), more=TRUE, BB)
print(position=c(0.5, 0.67, 1.0, 1.00), more=TRUE, CC)
print(position=c(0.5, 0.35, 1.0, 0.68), more=TRUE, DD)
print(position=c(0.5, 0.00, 1.0, 0.35), more=FALSE, EE)
## hhdev.off()


###################################################
### code chunk number 12: grap.tex:735-740
###################################################
## hhpdf("grap-f3.pdf", width=5.5, height=5.5)
splom( ~ tv[,c(4,5,1,2,3)],
      main=list("Televisions, Physicians, and Life Expectancy", cex=1.4),
      axis.text.cex=.5, pch=19, xlab=NULL, cex=.7, col=likertColor(2)[2])
## hhdev.off()


###################################################
### code chunk number 13: grap.tex:823-827
###################################################
## hhpdf("grap-f11a.pdf", width=9.9, height=8.58)
pairs(tv, pch=19, main="pairs with NW--SE diagonal and rectangular panels",
      cex.main=1.6)
## hhdev.off()


###################################################
### code chunk number 14: grap.tex:855-893
###################################################
## hhpdf("grap-f12.pdf", width=9, height=5)
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
     axis(2, at=1:4, labels=4:1, cex=1.2, las=1)
text(f12b[-dd], y=5-row(f12b)[-dd], x=col(f12b)[-dd], cex=2)
text(f12b[ dd], y=5-row(f12b)[ dd], x=col(f12b)[ dd], cex=1.4)
abline(a=5, b=-1)
abline(a=-2, b=1, lty=3, lwd=3)
abline(a=-1, b=1, lty=2)
abline(a= 0, b=1, lty=2)
abline(a= 1, b=1, lty=2)
abline(a= 2, b=1, lty=2)
  arrows(2.8, .75, 2.8,1.14, length=.1) ## up arrow
  arrows(3.8,1.7, 4.15,1.7, length=.1)  ## right arrow

plot(x=c(.5,4.5), y=c(.5,4.5), type="n", las=1,
     cex=1.2,
     main="b. single axis of symmetry.",
     xlab="variables", ylab="variables")
text(f12b[-dd], y=row(f12b)[-dd], x=col(f12b)[-dd], cex=2)
text(f12b[ dd], y=row(f12b)[ dd], x=col(f12b)[ dd], cex=1.4)
abline(a=0,b=1)
  arrows(2.8,3.75, 2.8,4.14, length=.1) ## up arrow
  arrows(3.8,2.7, 4.15,2.7, length=.1)  ## right arrow
par(old.par)
## hhdev.off()


###################################################
### code chunk number 15: grap.tex:972-977
###################################################
## hhpdf("grap-f5.pdf", width=6.5, height=6.5)
splom( ~ tv[, 1:3],
      main=list("Televisions, Physicians, and Life Expectancy", cex=1.4),
      xlab=NULL, pch=19, cex=.9, varname.cex=1.3, col=likertColor(2)[2])
## hhdev.off()


###################################################
### code chunk number 16: grap.tex:1013-1018
###################################################
## hhpdf("grap-f6.pdf", width=6.5, height=6.5)
splom( ~ cbind(tv[,1,drop=FALSE], log(tv[, 2:3])),
      main=list("log(Televisions, Physicians), and Life Expectancy", cex=1.4),
      xlab=NULL, pch=19, cex=.9, varname.cex=1.3, col=likertColor(2)[2])
## hhdev.off()


###################################################
### code chunk number 17: grap.tex:1094-1163
###################################################
## hhpdf("grap-f8.pdf", width=11, height=4.15, lwd=4) ## lwd is not an argument for grDevices:::pdf

x <- seq(0, 2.5, length=126)
xy <- cbind(x=x, ladder.f(x+.001))

xyA <- xy
outOfBoundsA <- ((xyA[,-1] > 6.40) | (xyA[,-1] < -2.25))
xyA[,-1][outOfBoundsA] <- NA
AA <-
xyplot(`-1` +  `-0.5` + `0` + `0.5` + `1` + `2` ~ x, data=xyA, col=col3x2,
       xlim=c(0, 2.5), ylim=c(-2.2, 6.4), type="l",
       main=list("       a. Simple Powers\n       with Negative Reciprocals\n       (monotonic, wrong order)", cex=1.2),
       ## ylab=list(expression(over(x^p, scriptstyle(sign(p)))), rot=0, cex=1.5),
       ylab=list(expression(T[p]^symbol(' ')), rot=0, cex=1.5),
       scales=list(x=list(at=seq(0, 2.5, .5))),
       par.settings=list(clip=list(panel=FALSE))) +
         layer(panel.axis("right", at=xyA[126, -1], labels=names(xyA)[-1],
                          outside=TRUE, text.col=col3x2, text.cex=1.2,
                          line.col=col3x2, tck=.5)) +
           layer(panel.axis("right", at=-1.5, labels="p", outside=TRUE, line.col="transparent", tck=1.5))

xyB <- xy
xyB[,2:3] <- -xyB[,2:3]
outOfBoundsB <- ((xyB[,-1] > 6.40) | (xyB[,-1] < -2.25))
xyB[,-1][outOfBoundsB] <- NA
BB <-
xyplot(`-1` +  `-0.5` + `0` + `0.5` + `1` + `2` ~ x, data=xyB, col=col3x2,
       xlim=c(0, 2.5), ylim=c(-2.2, 6.4), type="l",
       main=list("         b. Simple Powers\n         with Positive Reciprocals\n         (not monotonic, wrong order)", cex=1.2),
       ## ylab=list(expression(x^p), rot=0, cex=1.5),
       ylab=list(expression(W[p]^symbol(' ')), rot=0, cex=1.5),
       scales=list(x=list(at=seq(0, 2.5, .5))),
       par.settings=list(clip=list(panel=FALSE))) +
         layer(panel.axis("right", at=xyB[126, -1], labels=names(xyB)[-1],
                          outside=TRUE, text.col=col3x2,  text.cex=1.2,
                          line.col=col3x2, tck=.5)) +
           layer(panel.axis("right", at=-1.5, labels="p", outside=TRUE, line.col="transparent", tck=1.5))


xyC <- cbind(x=x, ladder.fstar(x+.001))
outOfBoundsC <- ((xyC[,-1] > 6.40) | (xyC[,-1] < -2.25))
xyC[,-1][outOfBoundsC] <- NA
CC <-
xyplot(`-1` +  `-0.5` + `0` + `0.5` + `1` + `2` ~ x, data=xyC, col=col3x2,
       xlim=c(0, 2.5), ylim=c(-2.2, 6.4), type="l",
       main=list("       c. Scaled Powers\n\n       (monotonic, right order)", cex=1.2),
       ## ylab=list(expression(over(x^p - 1, scriptstyle(p))), rot=0, cex=1.5),
       ylab=list(expression(T[p]^symbol('*')), rot=0, cex=1.5),
       scales=list(x=list(at=seq(0, 2.5, .5))),
       par.settings=list(clip=list(panel=FALSE))) +
         layer(panel.axis("right", at=xyC[126, -1], labels=names(xyC)[-1],
                          outside=TRUE, text.col=col3x2,  text.cex=1.2,
                          line.col=col3x2, tck=.5)) +
           layer(panel.axis("right", at=-1.5, labels="p", outside=TRUE, line.col="transparent", tck=1.5))

print(AA, panel.width=list(.6, "npc"), split=c(1,1,3,1), more=TRUE)
print(BB, panel.width=list(.6, "npc"), split=c(2,1,3,1), more=TRUE)
print(CC, panel.width=list(.6, "npc"), split=c(3,1,3,1), more=FALSE)


## hhdev.off()
## I am showing this in the script file.  This is not in the book.
## The interesting feature is the use of factor.levels in the strip.custom().
update(c(AA,BB,CC,layout=c(3,1)),
       between=list(x=3), scales=list(alternating=FALSE), ylab=NULL,
       par.strip.text=list(lines=3, cex=1.5),
       strip=strip.custom(factor.levels=
                          c(AA$ylab[[1]], BB$ylab[[1]], CC$ylab[[1]])),
       main=NULL)


###################################################
### code chunk number 18: grap.tex:1263-1309
###################################################
## hhpdf("grap-f7.pdf", width=7, height=7)
data(tv)
le <- ladder.fstar(tv$life.exp, "LE^")
ppp <- ladder.fstar(tv$ppl.per.phys, "PPP^")

form <-
  as.formula(paste(paste("`", names(le),  "`", sep="", collapse="+"), "~",
                   paste("`", names(ppp), "`", sep="", collapse="+")))

tmp <- xyplot(form, data=cbind(le, ppp), outer=TRUE, layout=c(6,6),
              pch=19, cex=.7, scales="free", col=likertColor(2)[2],
              between=list(x=.5, y=.5), aspect=1,
              xlab=NULL, xlab.top="People per Physician\n",
              ylab="Life Expectancy")
tmp2 <- useOuterStrips(combineLimits(matrix.trellis(tmp, byrow=TRUE, nrow=6, ncol=6)))
dimnames(tmp2) <- list(names(ppp), names(le))

tmp3 <- update(tmp2, scales=list(at=NULL, labels=NULL))

tmp3$strip
tmp3env <- environment(tmp3$strip)
tmp3$strip.left

tmp3$strip <-
function (which.given, which.panel, var.name, factor.levels, ...)
{
  if (which.given == 1)
    strip(which.given = 1, which.panel = which.panel[1],
          var.name = var.name[1],
          factor.levels=parse(text=dimnames(tmp3)[[1]]), ...)
}

tmp3$strip.left <-
  function (which.given, which.panel, var.name, factor.levels, ...)
{
  if (which.given == 2)
    strip.left(which.given = 1, which.panel = which.panel[2],
               var.name = var.name[2],
               factor.levels=parse(text=dimnames(tmp3)[[2]]), ...)
}

environment(tmp3$strip) <- tmp3env
environment(tmp3$strip.left) <- tmp3env

update(tmp3, par.strip.text=list(cex=1.3, lines=.7))
## hhdev.off()


###################################################
### code chunk number 19: grap.tex:1334-1345
###################################################
## hhpdf("grap-f9.pdf", width=7, height=4)
x <- sort(tv[,"life.exp"])
y <- ladder.fstar(x)
form <- as.formula(paste(paste("`", names(y), "`", sep="", collapse="+"), "~ x"))
xyplot(form, data=cbind(x=x, y), outer=TRUE, layout=c(3,2),
       strip=strip.custom(factor.levels=
           parse(text=paste("x^", names(y), sep=""))),
       par.strip.text=list(cex=1.5), col=likertColor(2)[2],
       type="l", scales=list(alternating=1, y="free", rot=0),
       between=list(y=1), ylab=NULL)
## hhdev.off()


###################################################
### code chunk number 20: grap.tex:1395-1423
###################################################
y <- ladder.fstar(tv$ppl.per.phys)
g <- paste("ppp ^",names(y))
g <- ordered(g, g)
g <- g[col(y)]
gg <- y
gg[] <- ''
## hhpdf("grap-f10bw.pdf", width=7, height=3)
tmpbw <-
bwplot(unlist(y) ~ unlist(gg) | g, groups=g, col=col3x2,
       layout=c(6,1), xlab=NULL, ylab=NULL,
       scales="free", panel=panel.bwplot.superpose)
update(tmpbw, scales=list(at=NULL, labels=NULL),
       strip=strip.custom(factor.levels=
           parse(text=paste("ppp^", names(y), sep=""))),
           par.strip.text=list(cex=1.5, lines=1.5),
           par.settings=list(plot.symbol=list(pch=19)))
## hhdev.off()
##
## hhpdf("grap-f10sp.pdf", width=7, height=3)
tmpsp <-
stripplot(unlist(y) ~ unlist(gg) | g, groups=g, col=col3x2,
          layout=c(6,1), xlab=NULL, ylab=NULL, pch=19,
          scales=list(relation="free", cex=.8))
update(tmpsp, scales=list(at=NULL, labels=NULL),
       strip=strip.custom(factor.levels=
           parse(text=paste("ppp^", names(y), sep=""))),
           par.strip.text=list(cex=1.5, lines=1.5))
## hhdev.off()


###################################################
### code chunk number 21: grap.tex:1444-1456
###################################################
scales <- c(2, 3, 3, 2, 2, 2)
names(scales) <- names(y)
for (i in names(y)) assign(i, y[[i]])
##
for (i in names(y)) {
  i.Rout <- paste("grap-f10-", i, ".Rout", sep="")
##  i.stem <- paste("stem(-y$`", i, "`, scale=", scales[i],")", sep="")
  i.stem <- paste("stem(-`", i, "`, scale=", scales[i],")", sep="")
  hhcapture(i.Rout, i.stem)
}
## I manually took these six files and edited them.
## The edited file is included as part c of the figure.


###################################################
### code chunk number 22: grap.tex:1514-1546
###################################################
## hhpng("ColorComparison.png", height=70, width=1060) ## pixels

data(col3x2)
ColTrellis <- trellis.par.get("superpose.symbol")$col[1:6]  ## $

grid.newpage()

pushViewport(viewport(width=.5, x=.5, just=1))
draw.key(list(rectangles=list(col=col3x2,
                              cex=1.5, size=8, height=2.8, border=FALSE),
              columns=6,
              between=1),
           draw=TRUE)
popViewport()

pushViewport(viewport(width=.5, x=1, just=1))
draw.key(list(rectangles=list(col=ColTrellis,
                              cex=1.5, size=8, height=2.8, border=FALSE),
              columns=6,
              between=1),
           draw=TRUE)
popViewport()

## hhdev.off()

## Open the file ColorComparison.png in ij.ImageJ.
## Open the "Vischeck Panel" in the ij.ImageJ plugin menu.
## Click the file and then the Vischeck button for one of
## the color vision deficieny types.
## Then click the original again and the next button.
## Then click the original again and the next button.



###################################################
### code chunk number 23: grap.tex:1953-1981
###################################################
bwdata <- data.frame(Y=(rt(80, df=5)*5 + rep(c(20,25,15,22, 22,28,16,14), 10)),
                     week=ordered(rep(c(1:4, 1:4), each=10)),
                     treatment= rep(c("A", "B"), each=40))

position(bwdata$week) <- c(1, 2, 4, 8)
levels(bwdata$week) <- c(1, 2, 4, 8)

bwdata$week.treatment <- with(bwdata, interaction(treatment, week))
position(bwdata$week.treatment) <-
   as.vector(t(outer(c(1, 2, 4, 8), c(-.18,.18), "+")))

BR <- likertColor(2, colorFunctionOption="default")[2:1]

## uses panel.bwplot.intermediate.hh to control position and colors
## hhpdf("bwplotposition.pdf", width=7, height=5)
bwplot(Y ~ week.treatment, data=bwdata,
       panel=panel.bwplot.intermediate.hh, xlim=c(0, 9),
       box.width=.25,
       pch=c(17, 16), col=BR,
       xlab="Week", ylab=list(rot=0),
       scales=list(x=list(at=position(bwdata$week), tck=1)),
       key=list(
          text=list(c("A","B"), col=BR),
          points=list(pch=c(17, 16), col=BR),
          space="top", columns=2, between=1, border=TRUE,
          title="Treatment", cex.title=.9)) +
   layer(panel.abline(h=37, col="gray60", lty=3, lwd=2))
## hhdev.off()


###################################################
### code chunk number 24: grap.tex:2118-2152
###################################################
tmp <- data.frame(x=1:6, y1=1:6, a=letters[c(1:3,1:3)], b=letters[c(4,4,4,5,5,5)],
                  y2=2:7, y3=3:8)

AA <- xyplot(y1 ~ x, data=tmp, pch=19, type="b", col=col3x2[1], main="AA: y1 ~ x")
BB <- xyplot(y2 ~ x, data=tmp, pch=19, type="b", col=col3x2[2], main="BB: y2 ~ x")
CC <- xyplot(y3 ~ x, data=tmp, pch=19, type="b", col=col3x2[3], main="CC: y3 ~ x")

## hhpdf("grap1.pdf", width=7, height=2)
print(AA, split=c(1,1,3,1), more=TRUE)
print(BB, split=c(2,1,3,1), more=TRUE)
print(CC, split=c(3,1,3,1), more=FALSE)
## hhdev.off()

## hhpdf("grap2.pdf", width=3, height=2)
update(AA + BB + CC, main="AA + BB + CC", ylim=c(.5, 8.5))
## hhdev.off()

## hhpdf("grap3.pdf", width=3, height=2)
update(c(y1=AA, y2=BB, y3=CC, layout=c(3,1), y.same=TRUE), main="c(AA, BB, CC)", scales=list(alternating=FALSE))
## hhdev.off()

## hhpdf("grap4.pdf", width=3, height=2)
xyplot(y1 + y2 + y3 ~ x, data=tmp, pch=19, type="b", col=col3x2,
       main="y1 + y2 + y3 ~ x")
## hhdev.off()

## hhpdf("grap5.pdf", width=3, height=2)
xyplot(y1 + y2 + y3 ~ x, data=tmp, pch=19, type="b", col=col3x2[1], outer=TRUE, layout=c(3,1),
       main="y1 + y2 + y3 ~ x, outer=TRUE", scales=list(alternating=FALSE))
## hhdev.off()

## hhpdf("grap6.pdf", width=4, height=2.5)
useOuterStrips(xyplot(y1 ~ x | a + b, data=tmp, pch=19, main="y1 ~ x | a + b", col=likertColor(2)[2]))
## hhdev.off()


