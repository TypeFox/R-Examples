### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/dsgntwo.tex'

###################################################
### code chunk number 1: dsgntwo.tex:9-10
###################################################
library(HH)


###################################################
### code chunk number 2: dsgntwo.tex:13-18
###################################################
## the standard lattice color 2 is difficult for people with color deficient vision
data(col3x2)
## These colors look like a 3x2 color array when run through
## the vischeck simulator to see how they look for the three most
## common color vision deficiencies: Protanope, Deuteranope, Tritanope.


###################################################
### code chunk number 3: dsgntwo.tex:230-266
###################################################
data(yatesppl)

tmpsppl <-
within(yatesppl, {
  levels(blocks  ) <- paste("B", levels(blocks  ), sep="")
  levels(plots   ) <- paste("P", levels(plots   ), sep="")
  levels(subplots) <- paste("S", levels(subplots), sep="")
  levels(variety ) <- paste("V", levels(variety ), sep="")
  levels(nitrogen) <- paste("N",
                            (as.numeric(levels(yatesppl$nitrogen))-1)*2,
                            sep=".")
})
tmpsppl$value <- interaction(tmpsppl$variety,
                             tmpsppl$nitrogen, sep=":")

tmp <-
reshape2::acast(tmpsppl, subplots ~ plots ~ blocks)

## hhcapture("yatesppl-layout.Rout", echo=FALSE, '
## blocks 1 2 3 4 5 6
tmp
## ')

tmp2 <- as.character(tmp)
dim(tmp2) <- dim(tmp)
dimnames(tmp2) <- dimnames(tmp)

## hhcapture("yatesppl-layout.123.Rout", echo=FALSE, '
## blocks 1 2 3
print(tmp2[,,1:3], quote=FALSE)
## ')

## hhcapture("yatesppl-layout.456.Rout", echo=FALSE, '
## blocks 4 5 6
print(tmp2[,,4:6], quote=FALSE)
## ')


###################################################
### code chunk number 4: dsgntwo.tex:292-300
###################################################
## hhpdf("yatesppl.pdf", height=7, width=7, col=col3x2) ## col is not an argument for grDevices:::pdf
position(tmpsppl$variety) <- (1:3) + c(.5, .6, .7)
tmpsppl$nit.lev <- factor(tmpsppl$nitlev)
levels(tmpsppl$nit.lev)[1] <- "0.0"
interaction2wt(y ~ variety + nit.lev, data=tmpsppl,
               par.strip.text=list(cex=1.4),
               main.cex=1.6)
## hhdev.off()


###################################################
### code chunk number 5: dsgntwo.tex:373-379
###################################################
## hhcapture("yatesppl-1.Rout", '
yatesppl.aov <-
   aov(y ~ variety*nitrogen + Error(blocks/plots/subplots),
       data=yatesppl)
summary(yatesppl.aov)
## ')


###################################################
### code chunk number 6: dsgntwo.tex:393-397
###################################################
## hhcapture("yatesppl-2.Rout", '
model.tables(yatesppl.aov, type="means")
model.tables(yatesppl.aov, type="effects", se=TRUE)
## ')


###################################################
### code chunk number 7: dsgntwo.tex:412-419
###################################################
## hhcapture("yatesppl-3.Rout", '
yatesppl.wrong.aov <-
  aov(terms(y ~ (blocks*variety) + (nitrogen*variety),
            keep.order=TRUE),
      data=yatesppl)
summary(yatesppl.wrong.aov)
## ')


###################################################
### code chunk number 8: dsgntwo.tex:481-487
###################################################
## hhcapture("yatesppl-alt.Rout", '
yatesppl2.anova <-
  aov(y ~ variety*nitrogen + Error(blocks/variety/nitrogen),
      data=yatesppl)
summary(yatesppl2.anova)
## ')


###################################################
### code chunk number 9: dsgntwo.tex:558-572
###################################################
## hhcapture("yatesppl-st6.Rout", '
## polynomial contrasts in nitrogen
contrasts(yatesppl$nitrogen)
contrasts(yatesppl$nitrogen) <- contr.poly(4)
contrasts(yatesppl$nitrogen)

## split plot analysis with polynomial contrasts
yatespplp.aov <-
  aov(y ~ variety*nitrogen + Error(blocks/plots/subplots),
      data=yatesppl)
summary(yatespplp.aov,
        split=list(nitrogen=list(linear=1, quad=2, cub=3)),
        expand.split=FALSE)
## ')


###################################################
### code chunk number 10: dsgntwo.tex:597-605
###################################################
## hhcapture("yatesppl-mmc.Rout", '
yatesppl.mmc <- mmc(yatesppl.wrong.aov, focus="nitrogen")
nitrogen.lmat <- contr.poly(4)
rownames(nitrogen.lmat) <- levels(yatesppl$nitrogen)
yatesppl.mmc <- mmc(yatesppl.wrong.aov, focus="nitrogen",
                    focus.lmat=nitrogen.lmat)
yatesppl.mmc
## ')


###################################################
### code chunk number 11: dsgntwo.tex:624-632
###################################################
## hhpdf("yatesppl-mmc.pdf", width=9, height=4.5)
yatesmmcpair <- mmcplot(yatesppl.mmc)
yatesmmcpoly <- mmcplot(yatesppl.mmc, type="lmat")
update(c("Pairwise Contrasts"=yatesmmcpair,
         "Polynomial Contrasts"=yatesmmcpoly,
         layout=c(2,1)),
       between=list(x=3), scales=list(x=list(alternating=FALSE)))
## hhdev.off()


###################################################
### code chunk number 12: dsgntwo.tex:760-764
###################################################
## hhcode("Design_2.8-2.R", '
data(Design_2.8_2)
Design_2.8_2
## ')


###################################################
### code chunk number 13: dsgntwo.tex:857-873
###################################################
## hhcode("2.8-2.R", '
R282 <- t(sapply(strsplit(Design_2.8_2$trt,""),
                function(trtcomb)
                   as.numeric(letters[1:8] %in% trtcomb)))
dimnames(R282) <- list(Design_2.8_2$trt, letters[1:8])
R282 <- data.frame(blocks=Design_2.8_2$blocks, R282)
R282
data(R282.y) ## R282.y was randomly generated
R282.aov <- aov(R282.y ~ blocks + (a+b+c+d+e+f+g+h)^2, data=R282)
anova(R282.aov)
model.matrix(R282.aov)
## confirm aliasing
R282E.aov <- aov(R282.y ~ Error(blocks) + (a+b+c+d+e+f+g+h)^2,
                 data=R282)
summary(R282E.aov)
## ')


###################################################
### code chunk number 14: dsgntwo.tex:928-932
###################################################
## hhpdf("circuit.pdf", col=likertColor(2)[2:1]) ## col is not an argument for grDevices:::pdf
data(circuit)
interaction2wt(yield ~ A*B*C, data=circuit, main.cex=1.6, box.width=.3)
## hhdev.off()


###################################################
### code chunk number 15: dsgntwo.tex:946-953
###################################################
##   hhcapture("circuit.Rout", '
circuit.aov <- aov( yield ~ A + B + C + A:B, data=circuit)
summary(circuit.aov)
model.tables(circuit.aov, type="means")
tapply(circuit[,"yield"], circuit[,"D"], mean)
tapply(circuit[,"yield"], circuit[,c("A","C")], mean)
## ')


###################################################
### code chunk number 16: dsgntwo.tex:1121-1131
###################################################
## hhpdf("cc135-fig.pdf", width=7, height=2.5)
data(cc135)
LLL <- xyplot(yield ~ period | cow, data=cc135,
              type="l", lty=3, col="gray50")
TTT <- xyplot(yield ~ period | cow, data=cc135,
              group=treat, pch=levels(cc135$treat), cex=2.5, col=col3x2)
update(LLL+TTT,
       layout=c(6,1), between=list(x=c(1)), scales=list(alternating=FALSE),
       strip=strip.custom(strip.names=c(TRUE,TRUE)), ylim=c(-15, 145))
## hhdev.off()


###################################################
### code chunk number 17: dsgntwo.tex:1206-1220
###################################################
##   hhcapture("cc135.Rout", '
data(cc135)
a1c <-  aov(terms(yield ~ cow + square:period + treat + res.treat,
                  keep.order=TRUE), data=cc135)
summary(a1c)
model.tables(a1c, type="means")
apply(summary(a1c)[[1]][,1:2], 2, sum)

a1cr <- aov(terms(yield ~ cow + square:period + res.treat + treat,
                  keep.order=TRUE), data=cc135)
summary(a1cr)
model.tables(a1cr, type="means")
apply(summary(a1cr)[[1]][,1:2], 2, sum)
## ')


###################################################
### code chunk number 18: dsgntwo.tex:1222-1232
###################################################
##  hhcode("cc135.R", '
data(cc135)
a1c <-  aov(terms(yield ~ cow + square:period + treat + res.treat,
                  keep.order=TRUE), data=cc135)
summary(a1c)

a1cr <- aov(terms(yield ~ cow + square:period + res.treat + treat,
                  keep.order=TRUE), data=cc135)
summary(a1cr)
## ')


###################################################
### code chunk number 19: dsgntwo.tex:1300-1334
###################################################
## hhpdf("cc135fbwplot.pdf", width=8, height=5, col=col3x2) ## col is not an argument for grDevices:::pdf
## construct the yield adjusted for the blocking factors and res.treat
cc135.block.res.aov <- aov(terms(yield ~ square/(cow + period) + res.treat,
                                 keep.order=TRUE), data=cc135)
summary(cc135.block.res.aov)
cc135$y.adj.res <- mean(cc135$yield) + resid(cc135.block.res.aov)

a1cr.adj.res <- aov(terms(y.adj.res ~ square/(cow + period) + res.treat + treat,
                          keep.order=TRUE), data=cc135)
summary(a1cr.adj.res)

print(position = c(0,0,.47,1), more = TRUE,  # left
      bwplot(y.adj.res ~ treat, data=cc135, ylim=c(35, 85),
             scales=list(x=list(cex=.7), y=list(cex=1.4)), ylab=list(cex=1.4),
             main="treatment means adjusted\n for blocks and residual treatments",
             panel=panel.bwplot.superpose, groups=treat))


## construct the yield adjusted for the blocking factors and treat
cc135.block.treat.aov <- aov(terms(yield ~ square/(cow + period) + treat,
                                   keep.order=TRUE), data=cc135)
summary(cc135.block.treat.aov)
cc135$y.adj.treat <- mean(cc135$yield) + resid(cc135.block.treat.aov)

a1cr.adj.treat <- aov(terms(y.adj.treat ~ cow + square/period + treat + res.treat,
                            keep.order=TRUE), data=cc135)
summary(a1cr.adj.treat)

print(position = c(.47,0,1,1), more = FALSE,  # right
      bwplot(y.adj.treat ~ res.treat, data=cc135, ylim=c(35, 85),
             scales=list(x=list(cex=.7), y=list(cex=1.4)), ylab=list(cex=1.4),
             main="residual treatment means\n adjusted for blocks and treatments",
             panel=panel.bwplot.superpose, groups=res.treat))
## hhdev.off()


###################################################
### code chunk number 20: dsgntwo.tex:1436-1474
###################################################
## Set up the individual ancova models.
## The graphs will be displayed in the next Scode section.
## The anova tables will be displayed in the following Scode section.

data(apple)

apple$yield.block <- apple$yield - proj(lm(yield ~ block, data=apple))[,'block']
apple$pre.block   <- apple$pre   - proj(lm(  pre ~ block, data=apple))[,'block']

`y.b~p.b*t`<-  ancovaplot(yield.block ~ pre.block*treat, data=apple,
                          groups=block, pch=letters[1:4], cex=1.6,
                          col=col3x2,
                          col.by.groups=FALSE)
## `y.b~p.b*t`

`y.b~p.b+t` <-  ancovaplot(yield.block ~ pre.block + treat, data=apple, col=col3x2)
## `y.b~p.b+t`

apple.aov.4 <- aov(yield.block ~ pre.block + treat, data=apple)

apple$yield.block.pre <-
  apple$yield.block - predict.lm(apple.aov.4, type="terms", terms="pre.block")

`y.b~t` <-  ancovaplot(yield.block ~ treat, x=pre.block, data=apple, col=col3x2)
## `y.b~t`

`y.b.p~t` <-  ancovaplot(yield.block.pre ~ treat, x=pre.block, data=apple, col=col3x2)
## `y.b.p~t`

`y.b~p.b` <-  ancovaplot(yield.block ~ pre.block, groups=treat, data=apple, col=col3x2)
## `y.b~p.b`

`y~b+p*t` <-  ancovaplot(yield ~ pre * treat, data=apple,
                         groups=block, pch=letters[1:4], cex=1.6,
                         col=col3x2,
                         condition=apple$treat,
                         col.by.groups=FALSE)
## `y~b+p*t`


###################################################
### code chunk number 21: dsgntwo.tex:1480-1538
###################################################

## convert each vector of panels into a row matrix
FF <- transpose(useOuterStrips(as.matrix(`y.b.p~t`  )))
EE <- transpose(useOuterStrips(as.matrix(`y.b~p.b+t`)))
DD <- transpose(useOuterStrips(as.matrix(`y.b~t`    )))
CC <- transpose(useOuterStrips(as.matrix(`y.b~p.b`  )))
BB <- transpose(useOuterStrips(as.matrix(`y.b~p.b*t`)))
AA <- transpose(useOuterStrips(as.matrix(`y~b+p*t`  )))
AA.strip <- AA$strip ## save top strip with one line after useOuterStrips

## with these row names
dimnames(FF)[[2]] <- "y.b.p~t"
dimnames(EE)[[2]] <- "y.b~p.b+t"
dimnames(DD)[[2]] <- "y.b~t"
dimnames(CC)[[2]] <- "y.b~p.b"
dimnames(BB)[[2]] <- "y.b~p.b*t"
dimnames(AA)[[2]] <- "y~b+p*t"

## align panel heights and axes
AA <- update(AA, ylim=c(161.1, 402.9), xlim=c(5.15, 11.23), scales=list(y=list(at=c(200,250,300,350)), x=list(at=c(5,7,9,11), alternating=0, tck=0     )), par.strip.text=list(lines=.7), between=list(x=c(.5,.5,.5,.5,.5,1.2), y=1), strip=TRUE , ylab=" ", xlab=NULL, xlab.top="Treatment", cex=2.2)
BB <- update(BB, ylim=c(214.4, 343.1), xlim=c(5.15, 11.23), scales=list(y=list(at=c(200,250,300,350)), x=list(at=c(5,7,9,11), alternating=0, tck=0     )), par.strip.text=list(lines=.7), between=list(x=c(.5,.5,.5,.5,.5,1.2), y=1), strip=FALSE, ylab=" ", xlab=NULL, cex=2.2)
CC <- update(CC, ylim=c(214.4, 343.1), xlim=c(5.15, 11.23), scales=list(y=list(at=c(200,250,300,350)), x=list(at=c(5,7,9,11), alternating=0, tck=0     )), par.strip.text=list(lines=.7), between=list(x=c(.5,.5,.5,.5,.5,1.2), y=1), strip=FALSE, ylab="yield (adjusted as indicated", xlab=NULL, cex=1.3)
DD <- update(DD, ylim=c(214.4, 343.1), xlim=c(5.15, 11.23), scales=list(y=list(at=c(200,250,300,350)), x=list(at=c(5,7,9,11), alternating=0, tck=0     )), par.strip.text=list(lines=.7), between=list(x=c(.5,.5,.5,.5,.5,1.2), y=1), strip=FALSE, ylab=" ", xlab=NULL, cex=1.3)
EE <- update(EE, ylim=c(214.4, 343.1), xlim=c(5.15, 11.23), scales=list(y=list(at=c(200,250,300,350)), x=list(at=c(5,7,9,11), alternating=0, tck=0     )), par.strip.text=list(lines=.7), between=list(x=c(.5,.5,.5,.5,.5,1.2), y=1), strip=FALSE, ylab=" ", xlab=NULL, cex=1.3)
FF <- update(FF, ylim=c(214.4, 343.1), xlim=c(5.15, 11.23), scales=list(y=list(at=c(200,250,300,350)), x=list(at=c(5,7,9,11), alternating=1, tck=c(1,0))), par.strip.text=list(lines=.7), between=list(x=c(.5,.5,.5,.5,.5,1.2), y=1), strip=FALSE, ylab=" ", xlab="pre (adjusted as indicated)", cex=1.3)
AA$strip <- AA.strip  ## restore top strip with one line after useOuterStrips
environment(AA$strip) <- environment(AA$strip.left)

## These numbers make the figure look right.  I don't totally understand them.
## The 1.879 makes the inches between tickmarks 250 and 300 be constant over all six rows.
## > diff(c(161.1, 402.9))
## [1] 241.8
## > diff(c(214.4, 343.1))
## [1] 128.7
## > 241.8/128.7
## [1] 1.879

aa <- 2.0
bb <- 1
cc <- .3
dd <- 1.879
ee <- 3.0

bottom <- cumsum(c(0, aa+bb-3.5*cc, rep(bb+1.9*cc, 3), bb-3.0*cc))
top    <- cumsum(c(aa+bb, rep(bb, 4), dd+ee-3.0*cc))
bottom <- bottom/top[6]
top    <- top/top[6]

## print each row separately
## hhpdf("appleAncovaPlots.pdf", width=8.2, height=10)
print(AA, position=c(0, bottom[6], 1, top[6]), panel.height=list(1.879*5, "lines"), more=TRUE )
print(BB, position=c(0, bottom[5], 1, top[5]), panel.height=list(1.000*5, "lines"), more=TRUE )
print(CC, position=c(0, bottom[4], 1, top[4]), panel.height=list(1.000*5, "lines"), more=TRUE )
print(DD, position=c(0, bottom[3], 1, top[3]), panel.height=list(1.000*5, "lines"), more=TRUE )
print(EE, position=c(0, bottom[2], 1, top[2]), panel.height=list(1.000*5, "lines"), more=TRUE )
print(FF, position=c(0, bottom[1], 1, top[1]), panel.height=list(1.000*5, "lines"), more=FALSE)
## hhdev.off()



###################################################
### code chunk number 22: dsgntwo.tex:1641-1644
###################################################
##    hhcapture("appleAncovaTable1.Rout", echo=FALSE, '
  anova(aov(yield ~ block + pre * treat, data=apple))
##  ')


###################################################
### code chunk number 23: dsgntwo.tex:1695-1714
###################################################
## hhpdf("apple-y-p.pdf", height=4, width=6, col=likertColor(2)[2]) ## col is not an argument for grDevices:::pdf
y.b <- bwplot(yield + yield.block ~ block, data=apple, outer=TRUE,
       par.settings=list(box.dot=list(
          col=trellis.par.get()$superpose.symbol$col[1])))
p.b <- bwplot(pre + pre.block ~ block, data=apple, outer=TRUE,
       par.settings=list(box.dot=list(
          col=trellis.par.get()$superpose.symbol$col[1])))

p.y.b <- matrix.trellis(c(p.b, y.b), nrow=2, ncol=2, byrow=TRUE)
pybos <- useOuterStrips(combineLimits(p.y.b))

pybos <- update(pybos,
                xlab="Blocks",
                xlab.top=c("Observed", "Adjusted for Blocks"),
                ylab=c("Pre", "Yield"),
                strip=FALSE, strip.left=FALSE,
                between=list(x=1, y=1))
pybos
## hhdev.off()


###################################################
### code chunk number 24: dsgntwo.tex:1763-1773
###################################################
  anova(aov(yield.block ~ pre.block * treat, data=apple))
##    hhcapture("appleAncovaTable2.Rout", '
apple$yield.block <-
  apple$yield - proj(lm(yield ~ block, data=apple))[,"block"]
apple$pre.block   <-
  apple$pre   - proj(lm(  pre ~ block, data=apple))[,"block"]
## wrong interaction sum of squares
anova(aov(yield.block ~ block + pre.block * treat,
      data = apple))
##  ')


###################################################
### code chunk number 25: dsgntwo.tex:1796-1805
###################################################
##    hhcapture("appleAncovaTable2qr.Rout", '
applebpst.aov <- aov(yield ~ block + pre * treat, data=apple,
                     x=TRUE)
appleQ <- qr.Q(qr(applebpst.aov$x))
`(pre.block:treat).block` <- appleQ[,11:15]
## correct anova for `y.b~p.b*t`
anova(aov(yield.block ~ pre.block + treat +
                        `(pre.block:treat).block`, data=apple))
##  ')


###################################################
### code chunk number 26: dsgntwo.tex:1838-1841
###################################################
##    hhcapture("appleAncovaTable3.Rout", echo=FALSE, '
  anova(aov(yield.block ~ pre.block, data=apple))
##  ')


###################################################
### code chunk number 27: dsgntwo.tex:1863-1866
###################################################
##    hhcapture("appleAncovaTable4.Rout", echo=FALSE, '
  anova(aov(yield.block ~ treat, data=apple))
##  ')


###################################################
### code chunk number 28: dsgntwo.tex:1910-1913
###################################################
##    hhcapture("appleAncovaTable5.Rout", echo=FALSE, '
  anova(aov(yield.block ~ pre.block + treat, data=apple))
##  ')


###################################################
### code chunk number 29: dsgntwo.tex:1930-1933
###################################################
##    hhcapture("appleAncovaTable5b.Rout", '
anova(aov(yield ~ block + pre + treat, data=apple))
##  ')


###################################################
### code chunk number 30: dsgntwo.tex:1950-1953
###################################################
##    hhcapture("appleAncovaTable5c.Rout", '
anova(aov(yield.block ~ block + pre.block + treat, data=apple))
##  ')


###################################################
### code chunk number 31: dsgntwo.tex:2026-2046
###################################################
## hhpdf("appleshift5.pdf", height=3, width=10)
update(`y.b~p.b+t`[1:6],
       cex=2,
       xlim=BB$x.limits,
       ylim=BB$y.limits,
       panel=panel.ancova.superpose,
       plot.resids=TRUE, print.resids=TRUE,
       mean.x.line=mean(unlist(sapply(`y.b~p.b+t`$panel.args, `[[`, "x"))),
       col.mean.x.line="gray40")
## hhdev.off()
## hhpdf("appleshift6.pdf", height=3, width=10)
update(`y.b.p~t`[1:6],
       cex=2,
       xlim=BB$x.limits,
       ylim=BB$y.limits,
       panel=panel.ancova.superpose,
       plot.resids=TRUE, print.resids=TRUE,
       mean.x.line=mean(unlist(sapply(`y.b.p~t`$panel.args, `[[`, "x"))),
       col.mean.x.line="gray40")
## hhdev.off()


###################################################
### code chunk number 32: dsgntwo.tex:2062-2069
###################################################
  anova(aov(yield.block.pre ~ treat, data=apple))
##    hhcapture("appleAncovaTable6.Rout", '
apple.aov.4 <-  aov(yield.block ~ pre.block + treat, data=apple)
apple$yield.block.pre <- apple$yield.block -
        predict.lm(apple.aov.4, type="terms", terms="pre.block")
anova(aov(yield.block.pre ~ treat, data = apple))
##  ')


###################################################
### code chunk number 33: dsgntwo.tex:2108-2116
###################################################
##    hhcapture("apple2.2mc.Rout", '
apple5.aov <- aov(yield ~ block + pre + treat, data=apple)
anova(apple5.aov)
apple5d.mmc <- mmc(apple5.aov,
    linfct=mcp(treat=contrMat(table(apple$treat),
                              type="Dunnett", base=6)))
apple5d.mmc
##  ')


###################################################
### code chunk number 34: dsgntwo.tex:2135-2139
###################################################
## hhpdf("appleMMC.pdf", height=6, width=6)
mmcplot(apple5d.mmc, main="Dunnett comparisons against Control=6", style="both",
        sub=list("\n             The MMC panel shows informative overprinting.  Please see Tiebreaker panel and caption.", cex=.75))
## hhdev.off()


###################################################
### code chunk number 35: dsgntwo.tex:2445-2463
###################################################
## hhpdf("testscorefsplom1.pdf", width=8, height=9)
## testscore.s
## testscore data:
## P. O. Johnson and F. Tsao, 1945
## R. L. Anderson and T. A. Bancroft, 1952

data(testscore)
## shorter level names for splom
levels(testscore$standing)[2] <- "avg"
levels(testscore$order)[2] <- "med"

splom( ~ testscore, pch=16,
      main="Original ordering of factor values", axis.text.cex=.6, xlab=NULL,
      col=likertColor(2)[2])
## restore longer level names
levels(testscore$standing)[2] <- "average"
levels(testscore$order)[2] <- "medium"
## hhdev.off()


###################################################
### code chunk number 36: dsgntwo.tex:2480-2497
###################################################
## hhpdf("testscorefsplom2.pdf", width=8, height=9)
## reorder the levels of the factors, and the order of the variables
## to improve the simplicity of the graph.
testscore$standing  <- ordered(testscore$standing,
                              levels=c("poor", "average", "good"))
testscore$order     <- ordered(testscore$order,
                              levels=c("low", "medium", "high"))
## shorter level names for splom
levels(testscore$standing)[2] <- "avg"
levels(testscore$order)[2] <- "med"
splom( ~ testscore[,c(1,6,7,4,5,3,2)], pch=16,
      main="Revised ordering of factor values", axis.text.cex=.6, xlab=NULL,
      col=likertColor(2)[2])
## restore longer level names
levels(testscore$standing)[2] <- "average"
levels(testscore$order)[2] <- "medium"
## hhdev.off()


###################################################
### code chunk number 37: dsgntwo.tex:2539-2578
###################################################
##   hhcapture("testscoreComp.Rout", '
## factors only
testscore1.aov <- aov(final ~ (sex + grade + standing + order)^2,
                      data=testscore)
summary(testscore1.aov)


## continuous first
testscore2.aov <- aov(final ~ initial + mental.age +
                      (sex + grade + standing + order)^2,
                      data=testscore)
summary(testscore2.aov)


## continuous second
testscore3s.aov <-
  aov(terms(final ~
            sex + grade + standing + order +
            sex:grade + sex:standing + sex:order +
            grade:standing + grade:order + standing:order +
            initial + mental.age,
            keep.order=TRUE),
      data=testscore)
summary(testscore3s.aov)


## continuous only
testscore4.aov <- aov(final ~ initial + mental.age,
                       data=testscore)
summary(testscore4.aov)


## comparisons
## factor vs both
anova(testscore1.aov, testscore3s.aov)

## continuous vs both
anova(testscore4.aov, testscore3s.aov)
## ')


###################################################
### code chunk number 38: dsgntwo.tex:2640-2664
###################################################
old.stars <- options(show.signif.stars=FALSE)
##   hhcapture("testscore5.Rout", '
## after looking at all of above
## Total Sum of Squares
var(testscore$final) * (length(testscore$final)-1)

testscore5.aov <- aov(final ~ initial + mental.age +
                      grade + standing + order + sex
                      + sex:order,
                      data=testscore)
summary(testscore5.aov)


testscore6.aov <- aov(final ~ initial + mental.age +
                      standing + order + grade + sex,
                      data=testscore)
summary(testscore6.aov)

testscore7.aov <- aov(final ~ initial + mental.age +
                      standing + order + sex,
                      data=testscore)
summary(testscore7.aov)
## ')
options(old.stars)


###################################################
### code chunk number 39: dsgntwo.tex:2711-2742
###################################################
testscore$final.adj <-
  testscore$final - as.vector(apply(proj(testscore7.aov, onedf=TRUE)[,2:3],1,sum))

yso.testscore.plot <- function(yname, main="Main Title Goes Here", cex=1.2, ...) {
  jj <- jitter(as.numeric(testscore$standing), amount=.1)
  formula <- as.formula(paste(yname, "~ jj | order"))
  xyplot(formula,
         groups=sex,
         main=main, cex=cex,
         xlab="standing", xlab.top="order",
         auto.key=list(border=TRUE, reverse=TRUE, space="right", title="sex"),
         par.settings=list(superpose.symbol=list(pch=c(19,17), cex=cex,
                             col=trellis.par.get()$superpose.symbol$col[2:1])),
         scales=list(x=list(at=1:3, labels=levels(testscore$standing),
                            alternating=1)),
         layout=c(3,1), between=list(x=1),
         ...)
}

## hhpdf("initial.pdf", width=7, height=3.4, col=likertColor(2)[2:1]) ## col is not an argument for grDevices:::pdf
yso.testscore.plot("initial", "Initial Score", ylim=c(5,35), data=testscore)
## hhdev.off()

## hhpdf("observed.pdf", width=7, height=3.4, col=likertColor(2)[2:1]) ## col is not an argument for grDevices:::pdf
yso.testscore.plot("final", "Observed Data", ylim=c(5,35), data=testscore)
## hhdev.off()

## hhpdf("adjusted1.pdf", width=7, height=3.4, col=likertColor(2)[2:1]) ## col is not an argument for grDevices:::pdf
yso.testscore.plot("final.adj", "Adjusted for Covariates", ylim=c(5,35), data=testscore)
## hhdev.off()



###################################################
### code chunk number 40: dsgntwo.tex:2779-2793
###################################################
## hhcapture("testscore7.Rout", '
newdata <- cbind(initial=mean(testscore$initial),
                 mental.age=mean(testscore$mental.age),
                 testscore[c(1:9,28:36),
                           c("standing","order","sex")])
newdata[c(1,2,3,4,18),]
final.pred <- predict(testscore7.aov, newdata=newdata)
final.pred.table <- tapply(final.pred, newdata[,3:5], c)
final.pred.table

## now summarize this over each factor to get predicted values
apply(final.pred.table, 1, mean) ## each scholastic standing
apply(final.pred.table, 2, mean) ## each individual order
## ')


###################################################
### code chunk number 41: dsgntwo.tex:2808-2815
###################################################
## hhpdf("adjusted2.pdf", width=7, height=3.5, col=likertColor(2)[2:1]) ## col is not an argument for grDevices:::pdf
yso.testscore.plot("final.adj", "Adjusted for Covariates", ylim=c(15,25), data=testscore)
## hhdev.off()

## hhpdf("prediction.pdf", width=7, height=3.5, col=likertColor(2)[2:1]) ## col is not an argument for grDevices:::pdf
yso.testscore.plot("final.pred", "Final Predictions", ylim=c(15,25), data=newdata)
## hhdev.off()


###################################################
### code chunk number 42: dsgntwo.tex:2900-2911
###################################################
data(crash)
## hhpdf("crash-bar.pdf", height=3, width=8, col=likertColor(2)[2]) ## col is not an argument for grDevices:::pdf
barchart(crashrate ~ passengers | agerange, data=crash,
         reference=TRUE, origin=0, between=list(x=1),
         ## xlab.top="Age Range",
         xlab="Number of Passengers",
         main="Crash Rates by Driver Age and Passenger Presence per 10,000 Trips",
         layout=c(3,1),
         col=trellis.par.get()$superpose.symbol$col[1],
         border=trellis.par.get()$superpose.symbol$col[1])
## hhdev.off()


###################################################
### code chunk number 43: dsgntwo.tex:2921-2926
###################################################
## hhpdf("crash-interaction.pdf", height=5, width=8, col=col3x2) ## col is not an argument for grDevices:::pdf
interaction2wt(crashrate ~ agerange + passengers, data=crash,
               main="Crash Rates by Driver Age and Passenger Presence per 10,000 Trips",
               strip=FALSE)
## hhdev.off()


###################################################
### code chunk number 44: dsgntwo.tex:2952-2956
###################################################
##   hhcapture("crash.aov.Rout", '
crash.aov <- aov(crashrate ~ agerange + passengers, data=crash)
summary(crash.aov)
## ')


###################################################
### code chunk number 45: dsgntwo.tex:2977-3014
###################################################
##   hhcapture("crash.cv.Rout", '
## means polish of crash data
tapply(crash$crashrate, crash[1:2], c)

mte <- model.tables(crash.aov, type="effects")
mte

mtm <- model.tables(crash.aov, type="means")
mtm

crash.resid <- tapply(resid(crash.aov), crash[,1:2], c)
crash.resid

## remove columns
cbind(rbind(crash.resid + as.vector(mte$tables$agerange),
            col=as.vector(mtm$tables$passengers)),
      row=0)


## remove rows
tmp <- rbind(cbind(t(t(crash.resid) + as.vector(mte$tables$passengers)),
                   row=as.vector(mtm$tables$agerange)),
             col=0)
zapsmall(tmp)


## removed both
cbind(rbind(crash.resid,
            col=as.vector(mte$tables$passengers)),
      row=c(as.vector(mte$tables$agerange), mtm[[1]][[1]][[1]]))

## cv: comparison value
cv <- outer(mte$tables$agerange,
            mte$tables$passengers) / mtm$tables$"Grand mean"[1]
dimnames(cv) <- list(levels(crash$agerange), levels(crash$passengers))
cv
## ')


###################################################
### code chunk number 46: dsgntwo.tex:3200-3206
###################################################
##   hhcapture("crash2.aov.Rout", '
crash2.aov <- aov(crashrate ~ agerange + passengers +
                  as.vector(cv), data=crash)
summary(crash2.aov)
coef(crash2.aov)
## ')


###################################################
### code chunk number 47: dsgntwo.tex:3227-3238
###################################################
## hhpdf("crash-diag.pdf", height=5, width=8)
## diagnostic plot. UREDA page 200--204
## crash-diag.eps.gz
crashr.lm <- lm(resid(crash.aov) ~ as.vector(cv))
anova(crashr.lm)
coef(crashr.lm)
plot(resid(crash.aov) ~ as.vector(cv), pch=16,
     main=paste("resid(crash.aov) ~ as.vector(cv)\nslope =",
       round(coef(crashr.lm)[2],4)))
abline(crashr.lm)
## hhdev.off()


###################################################
### code chunk number 48: dsgntwo.tex:3271-3318
###################################################
## hhpdf("crash-original.pdf", height=5, width=9, col=col3x2, lwd=2)
crash.a <-
interaction2wt(crashrate ~ agerange + passengers, data=crash,
               main="                  a. original scale,  k = 1",
               strip=FALSE)
crash.a$lattice.options$axis.options$right$labels2 <- "" ## remove right tick labels
crash.a$lattice.options$axis.options$right$labels3 <- c("","") ## remove right axis labels
crash.a
## hhdev.off()

## The regression coefficient of the cv term is 1.551647, from
## crashr.lm.  The recommended power for a transformation
## is 1-1.551647 = -0.551647
## We illustrate power = 0 -.5 -1

## hhpdf("crash-log.pdf", height=5, width=9, col=col3x2, lwd=2)
crash.b <-
interaction2wt(log(crashrate) ~ agerange + passengers,
               data=crash,
               main="                 b. log scale,  k = 0",
               strip=FALSE)
crash.b$lattice.options$axis.options$right$labels2 <- "" ## remove right tick labels
crash.b$lattice.options$axis.options$right$labels3 <- c("","") ## remove right axis labels
crash.b
## hhdev.off()

## hhpdf("crash-neg-rec-sqrt.pdf", height=5, width=9, col=col3x2, lwd=2)
crash.c <-
interaction2wt( I(-1/sqrt(crashrate)) ~ agerange + passengers,
               data=crash,
               main="             c. negative reciprocal square root scale,  k = -.5",
               strip=FALSE)
crash.c$lattice.options$axis.options$right$labels2 <- "" ## remove right tick labels
crash.c$lattice.options$axis.options$right$labels3 <- c("","") ## remove right axis labels
crash.c
## hhdev.off()

## hhpdf("crash-neg-rec.pdf", height=5, width=9, col=col3x2, lwd=2)
crash.d <-
interaction2wt( I(-1/(crashrate)) ~ agerange + passengers,
               data=crash,
               main="             d. negative reciprocal scale,  k = -1",
               strip=FALSE)
crash.d$lattice.options$axis.options$right$labels2 <- "" ## remove right tick labels
crash.d$lattice.options$axis.options$right$labels3 <- c("","") ## remove right axis labels
crash.d
## hhdev.off()


###################################################
### code chunk number 49: dsgntwo.tex:3352-3358
###################################################
## hhpdf("crash-interaction-rec.pdf", height=5, width=8, col=col3x2)
crash$crashrate.rec <- 10000/crash$crashrate
interaction2wt(crashrate.rec ~ agerange + passengers, data=crash,
               main="Trips per Crash by Driver Age and Passenger Presence",
               strip=FALSE, ylim=c(31000, 0))
## hhdev.off()


###################################################
### code chunk number 50: dsgntwo.tex:3376-3388
###################################################
## The appearance of the -.5 and -1 transformations are similar.
## We choose the reciprocal (power = -1) because it is easy to
## explain.  The units are "crashes per mile".  I dropped the
## negative in the anova table and table of means.  I left the
## negative in the graph so it would go in the same order as
## the other graphs.
##   hhcapture("crashi.aov.Rout", '
crashi.aov <-
   aov(10000/crashrate ~ agerange + passengers, data=crash)
summary(crashi.aov)
## ')
model.tables(crashi.aov, type="means")


###################################################
### code chunk number 51: dsgntwo.tex:3440-3450
###################################################
## hhpdf("crash-bar-rec.pdf", height=3, width=8, col=likertColor(2)[2]) ## col is not an argument for grDevices:::pdf
barchart(10000/crashrate ~ passengers | agerange, data=crash,
         reference=TRUE, origin=0, between=list(x=1),
         ## xlab.top="Age Range",
         xlab="Number of Passengers",
         main="Trips per Crash by Driver Age and Passenger Presence",
         layout=c(3,1),
         col=trellis.par.get()$superpose.symbol$col[1],
         border=trellis.par.get()$superpose.symbol$col[1])
## hhdev.off()


###################################################
### code chunk number 52: dsgntwo.tex:3495-3522
###################################################
## exploration of dummy variables
crashin.aov <- aov(1/crashrate ~ agerange/passengers, data=crash)
summary(crashin.aov)
## coef(summary.lm(crashin.aov))
summary(crashin.aov,
        split=list("agerange:passengers"=
          list(teens=1:2, adults=3, rest=4:9)))
## model.tables(crashin.aov, type="means")

## Now that we can see the different behavior for the
## passengers conditional on the agerange, let us make the
## anova table show it.
## Selection of just the linear contrasts.
old.width <- options(width=66)
##   hhcapture("crashinlin.aov.Rout", '
pass <- as.numeric(crash$passengers)
crashinlin.aov <- aov(10000/crashrate ~ agerange/pass,
                      data=crash)
## summary(crashinlin.aov)
## coef(summary.lm(crashinlin.aov))
print(coef(summary.lm(crashinlin.aov))[4:6,], digits=4)

print(digits=3,
summary(crashinlin.aov,
    split=list("agerange:pass"=list(teens=1:2, adults=3))))
## ')
options(old.width)


###################################################
### code chunk number 53: dsgntwo.tex:3767-3815
###################################################
##  hhcode("yatesppl.ex.R", '
### a. The whole plot column space is defined by the
###       plots %in% blocks
### dummy variables generated by the
       ## alternate residuals formula: orthogonal contrasts are critical
       data(yatesppl)
       yatesppl.resida.aov <- aov(y ~ blocks/plots,
                                  data=yatesppl, x=TRUE,
                                  contrasts=list(blocks=contr.helmert,
                                                 plots=contr.helmert))
       summary(yatesppl.resida.aov)
       t(yatesppl.resida.aov$x)
###
### b. This is the same column space defined by the
###       variety + blocks:variety
### dummy variables generated by the
       ## computational shortcut
       yatesppl.short.aov <-
         aov(terms(y ~ blocks + variety + blocks*variety +
                   nitrogen + variety*nitrogen,
                   keep.order=TRUE),  ## try it without keep.order=TRUE
             data=yatesppl, x=TRUE)
       summary(yatesppl.short.aov)
       t(yatesppl.short.aov$x)
###
### c. We illustrate this by regressing the response variable y on
### the  variety + blocks:variety dummy variables
       ## project y onto blocks/plots dummy variables
       plots.aov <- lm(y ~ yatesppl.resida.aov$x[,7:18], data=yatesppl)
       summary.aov(plots.aov)
       y.bp <- predict(plots.aov)
       variety.aov <- aov(y.bp ~ blocks*variety, data=yatesppl)
       summary(variety.aov)
### and seeing that we reproduce the plots %in% blocks
### stratum of the ANOVA table
###     Error: plots %in% blocks
###               Df Sum of Sq  Mean Sq F Value     Pr(F)
###       variety  2  1786.361 893.1806 1.48534 0.2723869
###     Residuals 10  6013.306 601.3306
### obtained from the complete five-factor specification.
###
       ## split plot analysis
       yatesppl.anova <- aov(y ~ variety*nitrogen +
                                 Error(blocks/plots/subplots),
                             data=yatesppl)
       summary(yatesppl.anova)
###
## ')


