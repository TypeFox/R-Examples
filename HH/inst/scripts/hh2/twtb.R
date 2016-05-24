### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/twtb.tex'

###################################################
### code chunk number 1: twtb.tex:7-8
###################################################
library(HH)


###################################################
### code chunk number 2: twtb.tex:11-16
###################################################
## the standard lattice color 2 is difficult for people with color deficient vision
data(col3x2)
## These colors look like a 3x2 color array when run through
## the vischeck simulator to see how they look for the three most
## common color vision deficiencies: Protanope, Deuteranope, Tritanope.


###################################################
### code chunk number 3: twtb.tex:80-81
###################################################
require(vcd)


###################################################
### code chunk number 4: twtb.tex:103-107
###################################################
## hhcapture("drunk.Rout", '
data(drunk)
drunk
## ')


###################################################
### code chunk number 5: twtb.tex:127-137
###################################################
BlyCol <- likertColor(12)[c(4,12)]
prop.female <- drunk["females",]/colSums(drunk)
ages <- ordered(dimnames(drunk)$age, levels=dimnames(drunk)$age)
## hhpdf("drunk-prop-fem.pdf", width=5.5, height=3)
barchart(prop.female ~ ages,
         horizontal=FALSE, origin=0,
         ylab="", main="proportion female",
         col=BlyCol[1],
         border=BlyCol[1])
## hhdev.off()


###################################################
### code chunk number 6: twtb.tex:150-160
###################################################
## hhpdf("drunk-mosaic.pdf", width=8, height=4)
mosaic(t(drunk), direction=c("v","h"),
       gp=gpar(fill=BlyCol[2:1], col="transparent"),
       rot_labels=c(0,0,0,0),  ## zero is horizontal
       rot_varnames=c(0,0,0,0),
       offset_labels=c(0, -0.6, 0, 1),  ## top, right, bottom, left ## positive means outward
       offset_varnames=c(0, -0.6, 0, 2.4),
       margins=c(left=6.5),
       keep_aspect_ratio=FALSE)
## hhdev.off()


###################################################
### code chunk number 7: twtb.tex:180-188
###################################################
## hhcapture("drunk2.Rout", '
drunk.chisq <- chisq.test(drunk)
drunk.chisq
drunk.chisq$observed
drunk.chisq$expected
drunk.chisq$residuals   ## cell chi values
drunk.chisq$residuals^2 ## cell chi-square values
## ')


###################################################
### code chunk number 8: twtb.tex:214-223
###################################################
## hhpdf("drunk-chi.pdf", width=5.5, height=3)
barchart(Freq ~ age | sex, as.data.frame(drunk.chisq$residuals),
         origin=0, layout=c(1,2), as.table=TRUE,
         scales=list(alternating=2), ##between=list(y=1),
         ylab=list("sex", rot=0),
         ylab.right=list("Chi values", rot=0), xlab="Age",
         strip=FALSE, strip.left=TRUE,
         col=BlyCol[2:1], border=BlyCol[2:1], groups=sex)
## hhdev.off()


###################################################
### code chunk number 9: twtb.tex:317-326
###################################################
## hhpdf("drunk-assoc.pdf", width=5.5, height=3.5)
assoc(drunk, gp=gpar(fill=BlyCol[rep(2:1, each=5)], col=0),
       margins=c(left=5),
      rot_labels=c(0,0,0,0), rot_varnames=c(0,0,0,0),
      just_labels=c("center","right","center","right"),
      just_varnames=c("center","right","center","left"),
      offset_labels=c(0, 0, 0, 0),
      offset_varnames=c(0, 0, 0, 2.5))
## hhdev.off()


###################################################
### code chunk number 10: twtb.tex:451-455
###################################################
## hhcapture("glasses.Rout", '
data(glasses)
glasses
## ')


###################################################
### code chunk number 11: twtb.tex:484-488
###################################################
## hhcapture("glasses2.Rout", '
fisher.test(glasses)
## ')
chisq.test(glasses, corr=FALSE)


###################################################
### code chunk number 12: twtb.tex:512-537
###################################################
## hhcapture("glasses-all.Rout", '
## construct all possible two-way tables with the same margins as the
## initial table
all.tables <- function(x) {

  xx <- x

  r.margin <- rowSums(x)
  c.margin <- colSums(x)

  result <- array(0, dim=c(r.margin[1]+1, 2, 2),
                  dimnames=c(table=list(0:r.margin[1]), rev(dimnames(xx))))
  for (x11 in 0:r.margin[1]) {
    xx[1,1] <- x11
    xx[1,2] <- r.margin[1] - xx[1,1]
    xx[2,1] <- c.margin[1] - xx[1,1]
    xx[2,2] <- sum(x) - (xx[1,1] + xx[1,2] + xx[2,1])
    if (min(xx) >= 0) result[as.character(x11),,] <- t(xx)
  }
  result
}

glasses.all <- all.tables(glasses)
aperm(glasses.all, c(3,2,1))
## ')


###################################################
### code chunk number 13: twtb.tex:558-564
###################################################
## hhpdf("glasses-exact.pdf", width=11, height=3)
(mosaic(glasses.all, direction=c("v","v","h"),
        gp=gpar(fill=BlyCol, col="transparent"),
        highlighting=3, highlighting_fill=likertColor(2),
        spacing=spacing_increase(rate=c(.2, 3, 3.5))))
## hhdev.off()


###################################################
### code chunk number 14: twtb.tex:584-618
###################################################
## hhcapture("glasses-hypergeometric.Rout", '
g.p <- apply(glasses.all, 1,
             function(x)
               c(prob=dhyper(x[1,1], sum(x[1,]), sum(x[2,]), sum(x[,1])),
                 min=min(x)))
names(dimnames(g.p))[1] <- ""
dimnames(g.p)[[1]] <- c("prob", "min") ## names in the apply FUN not picked up
g.p2 <- data.frame(t(g.p), which=I(""))
## initial table
g.p2[as.character(glasses[1,1]),"which"] <- "*"
## more extreme tables (min value is smaller)
g.p2[g.p2[as.character(glasses[1,1]),"min"] > g.p2[,"min"],"which"] <- "<"
g.p2
g.p2$cumsum <- cumsum(g.p2$prob)
g.p2$rev.cumsum <- rev(cumsum(rev(g.p2$prob)))
g.p2
## ')
##
## hhpdf("glasses-exact-prob.pdf", width=7, height=2.5)
barchart(g.p2$prob ~ factor(0:6), horizontal=FALSE,
         ylab=NULL, origin=0,
         main=list(labels=
           "probability of table with specified [1,1] position"),
         scales=list(x=list(at=1+0:6, labels=paste(0:6,g.p2$which))),
         xlab.top=list(format(round(g.p[1,], digits=4)), cex=.8),
         col=BlyCol[2],
         border=BlyCol[2],
         key=list(
           text=list(c("observed","more extreme")),
           text=list(c("*","<")),
           columns=2,
           border=TRUE,
           space="bottom"))
## hhdev.off()


###################################################
### code chunk number 15: twtb.tex:692-735
###################################################
## hhcapture("blyth.Rout", '
require(vcd)
require(reshape2)

data(blyth)

## rearrange as 3-way array
blyth3 <- blyth
dim(blyth3) <- c(2,2,2)
dimnames(blyth3) <- list(survival=c("not","survive"),
                         treatment=c("standard","new"),
                         location=c("A","B"))
blyth3x <- abind::abind(blyth3,
                        "A&B combined"=apply(blyth3, 1:2, sum))
names(dimnames(blyth3x)) <- names(dimnames(blyth3))
## blyth3x
structable(aperm(blyth3x, c(3,2,1)), direction=c("v","v","h"))


blyth3x.pct <- 100 * blyth3x /
  abind::abind(apply(blyth3x, c(2,3), sum),
               apply(blyth3x, c(2,3), sum), along=.5)
round(
structable(aperm(blyth3x.pct, c(3,2,1)),
           direction=c("v","v","h"))
)

blyth3xdf <- cbind(as.data.frame.table(blyth3x),
                   Pct=as.vector(blyth3x.pct))
blyth3xdf$Survival <-
  factor(blyth3xdf$survival,
         levels=rev(levels(blyth3xdf$survival)))
blyth3xdf

blyth3xdf.1.8 <- blyth3xdf[1:8,]
blyth3xdf.1.8$location <- factor(blyth3xdf.1.8$location)

blyth3xdf.9.12 <- blyth3xdf[9:12,]
blyth3xdf.9.12$location <- factor(blyth3xdf.9.12$location)

blyth3xc <- dcast(location + treatment ~ survival,
                  value.var="Freq", data=blyth3xdf)
## ')


###################################################
### code chunk number 16: twtb.tex:769-781
###################################################
## BlyCol <- likertColor(8)[c(2,4)]
BlyCol <- rainbow(12)[c(12,9)]
## hhpdf("bC3r.pdf", width=7, height=2.5)
resizePanels(w=c(.31,.31,.38),
barchart(Freq ~ treatment | location, groups=Survival, data=blyth3xdf,
         stack=TRUE,
         horizontal=FALSE, ylab="Count",
         ylab.right=list(c("survive","not"), rot=0),
         col=BlyCol[2:1], border=BlyCol[2:1],
         layout=c(3,1), between=list(x=c(0,2)))
             )
## hhdev.off()


###################################################
### code chunk number 17: twtb.tex:784-794
###################################################
## hhpdf("bP3r.pdf", width=7, height=2.5)
resizePanels(w=c(.31,.31,.38),
barchart(Pct ~ treatment | location, groups=Survival, data=blyth3xdf,
         stack=TRUE,
         horizontal=FALSE, ylab="Percent",
         ylab.right=list(c("survive","not"), rot=0),
         col=BlyCol[2:1], border=BlyCol[2:1],
         layout=c(3,1), between=list(x=c(0,2)))
             )
## hhdev.off()


###################################################
### code chunk number 18: twtb.tex:797-808
###################################################
## hhpdf("bP3s.pdf", width=7, height=2.5)
print(position=c(0, 0, .93, 1),
resizePanels(w=c(.31,.31,.38),
barchart(Pct ~ treatment | location, data=blyth3xdf,
         subset=(survival=="survive"),
         origin=0, ylim=c(-7, 107), col=BlyCol[2], border=BlyCol[2],#"white",
         ylab="Percent Survive",
         layout=c(3,1), between=list(x=c(0,2)))
             )
)
## hhdev.off()


###################################################
### code chunk number 19: twtb.tex:811-824
###################################################
## Figures mc2.pdf and mc1.pdf have the panel borders drawn by strucplot
## and the labeling by mosaic inside each panel.  They are positioned by LaTeX.
## hhpdf("mc2.pdf", width=6, height=3)
cotabplot(~ treatment + survival | location, data=blyth3xdf.1.8,
          layout = c(1, 2), direction=c("v","h"),
          rot_labels=c(0,0,0,0), rot_varnames=c(0,0,0,0), ## zero is horizontal
          offset_labels=c(0, -0.6, 0, .5),  ## top, right, bottom, left ## positive means outward
          offset_varnames=c(0, -0.6, 0, 2),
          panel_args=list(margins=c(4,1,2,6)),  ## need room for horizontal variable name
          keep_aspect_ratio=FALSE,
          spacing=spacing_highlighting(rate=6),
          gp=gpar(fill=BlyCol, col=0))
## hhdev.off()


###################################################
### code chunk number 20: twtb.tex:826-837
###################################################
## hhpdf("mc1.pdf", width=4.2, height=3)
cotabplot(~ treatment + survival | location, data=blyth3xdf.9.12,
          layout = c(1, 2), direction=c("v","h"),
          rot_labels=c(0,0,0,0), rot_varnames=c(0,0,0,0), ## zero is horizontal
          offset_labels=c(0, -0.6, 0, .5),  ## top, right, bottom, left ## positive means outward
          offset_varnames=c(0, -0.6, 0, 2),
          panel_args=list(margins=c(4,1,2,6)),  ## need room for horizontal variable name
          keep_aspect_ratio=FALSE,
          spacing=spacing_highlighting(rate=6),
          gp=gpar(fill=BlyCol, col=0))
## hhdev.off()


###################################################
### code chunk number 21: twtb.tex:840-863
###################################################
## Figures mc3a.pdf and mc3b.pdf have the panel borders drawn by lattice
## and the labeling outside all panels.  The figures are manually superposed using LaTeX.
## The offsets, margins and such are carefully tailored to these pdf settings.
mosaic.labels <- TRUE  ## to see a completely labeled mosaic plot
mosaic.labels <- FALSE ## to see a mosaic plot with all labels suppressed
## {if (mosaic.labels)
##    hhpdf("mc3.pdf", width=7, height=2.25)
## else
##   hhpdf("mc3a.pdf", width=7, height=2.25)
## }
## hhpdf(if (mosaic.labels) "mc3.pdf" else "mc3a.pdf", width=7, height=2.25)
mosaic(~ treatment + survival | location, data=blyth3xdf,
       layout = c(1, 3), direction=c("v","v","h"),
       rot_labels=c(0,0,0,0), rot_varnames=c(0,0,0,0), ## zero is horizontal
       offset_labels=c(0, -0.6, 0, 1.5),  ## top, right, bottom, left ## positive means outward
       offset_varnames=c(0, -0.6, 0, 2.8),
       varnames=mosaic.labels,
       labels=mosaic.labels,
       margins=c(left=6.5),
       keep_aspect_ratio=FALSE,
       spacing=spacing_highlighting(rate=3.5),
       gp=gpar(fill=BlyCol, col=0))
## hhdev.off()


###################################################
### code chunk number 22: twtb.tex:865-877
###################################################
## hhpdf("mc3b.pdf", width=8.5, height=2.5)
print(position=c(0, 0, .93, 1),
resizePanels(w=c(.28, .26, .46),
barchart(Pct ~ treatment | location, data=blyth3xdf,
         subset=(survival=="survive"),
         origin=0, ylim=c(-7, 107), col=0, border=0,
         ylab="Percent",
         ylab.right=list(c("survive","not"), rot=0),
         layout=c(3,1), between=list(x=c(0,2)))
             )
)
## hhdev.off()


###################################################
### code chunk number 23: twtb.tex:880-908
###################################################
## Figure mc3pdf.pdf is merged at the R level.
## The offsets, margins and such are carefully tailored to these pdf settings.
## hhpdf("mc3pdf.pdf", width=7, height=2.5)
##
mosaic(~ treatment + survival | location, data=blyth3xdf,
       layout = c(1, 3), direction=c("v","v","h"),
       rot_labels=c(0,0,0,0), rot_varnames=c(0,0,0,0), ## zero is horizontal
       offset_labels=c(0, -0.6, 0, 1.5),  ## top, right, bottom, left ## positive means outward
       offset_varnames=c(0, -0.6, 0, 2.8),
       varnames=mosaic.labels,
       labels=mosaic.labels,
       margins=c(top=3.45, right=6, bottom=2.935, left=5.5),
       keep_aspect_ratio=FALSE,
       spacing=spacing_highlighting(rate=3.5),
       gp=gpar(fill=BlyCol, col=0))
##
print(more=TRUE,
resizePanels(w=c(.305, .26, .435),
barchart(Pct ~ treatment | location, data=blyth3xdf,
         subset=(survival=="survive"),
         origin=0, ylim=c(-7, 107), col=0, border=0,
         ylab="Percent",
         ylab.right=list(c("survive","not"), rot=0),
         layout=c(3,1), between=list(x=c(0, 1)))
             )
)
##
## hhdev.off()


###################################################
### code chunk number 24: twtb.tex:911-917
###################################################
## hhpdf("lC3r.pdf", width=7, height=2.5)
likert(treatment ~ .| location, blyth3xc, horizontal=FALSE,
       main=NULL,  xlab=NULL, col=BlyCol,
       layout=c(3,1), between=list(x=c(0,2)), w.resizePanels=c(.31,.31,.38),
       ylab.right=list(c("not","survive"), rot=0), auto.key=FALSE)
## hhdev.off()


###################################################
### code chunk number 25: twtb.tex:920-926
###################################################
## hhpdf("lP3r.pdf", width=7, height=2.5)
likert(treatment ~ .| location, blyth3xc, horizontal=FALSE, as.percent=TRUE,
       main=NULL, xlab=NULL, col=BlyCol,
       layout=c(3,1), between=list(x=c(0,2)), w.resizePanels=c(.31,.31,.38),
       ylab.right=list(c("not","survive"), rot=0), auto.key=FALSE)
## hhdev.off()


###################################################
### code chunk number 26: twtb.tex:1352-1361
###################################################
## hhcapture("hypothermia.Rout", '
hypothermia <-
    matrix(c(75,54,61,83),
           nrow=2,
           dimnames=list(
             Treatment=c("treated","control"),
             Outcome=c("favorable","not.favorable")))
hypothermia
## ')


###################################################
### code chunk number 27: twtb.tex:1375-1379
###################################################
## hhpdf("hypothermiamosaic.pdf", width=4, height=3.5)
mosaic(Outcome ~ Treatment, data=as.data.frame.table(hypothermia), direction=c("v","h"),
       gp=gpar(fill=BlyCol[2:1], col="white"), keep_aspect_ratio=FALSE)
## hhdev.off()


###################################################
### code chunk number 28: twtb.tex:1390-1402
###################################################
## hhpdf("hypothermiacount.pdf", width=5, height=3)
## not included in book
hypothermia.df <- as.data.frame.table(hypothermia)
hypothermia.df$Outcome <-
  factor(hypothermia.df$Outcome, levels=rev(levels(hypothermia.df$Outcome)))

barchart(Freq ~ Treatment, groups=Outcome, stack=TRUE, hypothermia.df,
         horizontal=FALSE, origin=0, ylab="Count", col=BlyCol,
         auto.key=list(
           border=TRUE, space="right", rectangles=FALSE,
           reverse=TRUE, rect=list(col=BlyCol)))
## hhdev.off()


###################################################
### code chunk number 29: twtb.tex:1414-1426
###################################################
## hhpdf("hypothermiaproportion.pdf", width=5, height=3)
## not included in book
hypothermiaProportion.df <- as.data.frame.table(hypothermia / rep(colSums(hypothermia), each=2))
hypothermiaProportion.df$Outcome <-
  factor(hypothermiaProportion.df$Outcome, levels=rev(levels(hypothermiaProportion.df$Outcome)))

barchart(Freq ~ Treatment, groups=Outcome, stack=TRUE, hypothermiaProportion.df,
         horizontal=FALSE, origin=0, ylab="Proportion", col=BlyCol,
         auto.key=list(
           border=TRUE, space="right", rectangles=FALSE,
           reverse=TRUE, rect=list(col=BlyCol)))
## hhdev.off()


###################################################
### code chunk number 30: twtb.tex:1438-1451
###################################################
## hhpdf("hypothermiaodds.pdf", width=3.5, height=2.5)
barchart(hypothermia[,1] / hypothermia[,2] ~ dimnames(hypothermia)[[1]],
          horizontal=FALSE, origin=0, ylab="odds favorable",
         col=BlyCol[2],
         border=BlyCol[2])
## hhdev.off()
##
## hhpdf("hypothermialogit.pdf", width=3.5, height=2.5)
barchart(log(hypothermia[,1] / hypothermia[,2]) ~ dimnames(hypothermia)[[1]],
          horizontal=FALSE, origin=0, ylab="logit favorable",
         col=BlyCol[2],
         border=BlyCol[2])
## hhdev.off()


###################################################
### code chunk number 31: twtb.tex:1567-1580
###################################################
## hhpdf("hypothermiaplotOddsRatio.pdf", width=7.5, height=5, lwd=4)
tmp <- plotOddsRatio(t(hypothermia[2:1, 2:1]), col=col3x2)
## tmp

update(tmp, par.settings=list(clip=list(panel=FALSE), layout.widths=list(axis.right=.7, key.right=.9, axis.key.padding=0)),
        xlab.top=list("Given for Control", rot=0, just=.8), ylab.right=list(c(" ", " ", "Confidence\nInterval\non Predicted\nfor Treated"), rot=0, adjust=0)) +
layer(panel.abline(v=.3942, h=c(.5515,.4318,.6655), lty=2, lwd=1.75, col="gray50")) + ## works, use this one!
layer(panel.axis("right",  at=c(.5515, .4318, .6655), labels=FALSE,              line.col=col3x2[c(1,3,3)],               tck=1.3, outside=TRUE       )) + ## ticks
layer(panel.axis("right",  at=c(.5515, .4318, .6655),
                       labels=c(.5515, .4318, .6655), text.col=col3x2[c(1,3,3)], line.col="transparent",    text.cex=1.2, tck=1.3, outside=TRUE, rot=0)) + ## labels
layer(panel.axis("top",    at=.3942,                  text.col=col3x2[1],        line.col=col3x2[1],        text.cex=1.2, tck=0.9, outside=TRUE))

## hhdev.off()


###################################################
### code chunk number 32: twtb.tex:1724-1736
###################################################
data(salk)
salk.adjlabels <- salk
levels(salk.adjlabels$vaccine)
## [1] "no.vac" "vac"
levels(salk.adjlabels$vaccine) <- c("no.vac    ", "   vac") ## avoid overprinting x-axis labels

## hhpdf("salkMosaic.pdf", width=12, height=5.5) ## height=3.5 for portrait, height=5.5 for landscape
(mosaic(Freq ~ vaccine + paralyze | age, data=salk.adjlabels, direction=c("v","v","h"),
        main="Observed number of observations in each age group",
        gp=gpar(fill=BlyCol[2:1], col=0),
        spacing=spacing_increase(rate=c(.4, 1.4, 3.5))))
## hhdev.off()


###################################################
### code chunk number 33: twtb.tex:1755-1764
###################################################
## hhcapture("MHsalk.Rout", '
## Code for calculation of the Cochran--Mantel--Haenszel test of the polio example
salk2 <- tapply(salk$Freq, salk[c(2,3,1)], c)
class(salk2) <- "table"
salk2

mantelhaen.test(salk2)
mantelhaen.test(salk2, correct=FALSE)
## ')


###################################################
### code chunk number 34: twtb.tex:1768-1848
###################################################
## hhcapture("arithMHsalk.Rout", '
## Code for "Detail for calculation of the Cochran--Mantel--Haenszel test of the polio example."
## counts
salk2

## proportion without paralysis
pp <- apply(salk2, c(3,1),
            function(x) x[1]/(x[1]+x[2]))
pp

## binomial variance for proportion without paralysis
apply(salk2, c(3,1),
      function(x) (x[1]/(x[1]+x[2]))*(x[1]/(x[1]+x[2])) / (x[1]+x[2]))


## average proportion without paralysis
p <- apply(salk2, 3,
      function(x) sum(x[,1])/sum(x))
p

## weight per table
w <- apply(salk2, 3,
           function(x) 1/sum(1/(x[,1]+x[,2])))
w

## diff of proportion without paralysis
dp <- pp[,1] - pp[,2]
dp

## binomial variance for difference of proportions without paralysis
p*(1-p)

sum(w*dp) / sqrt(sum(w*p*(1-p)))


## chi-square for each table
chisq.table <-
t(apply(salk2, 3,
      function(x) {
        e <- (x[,1]+x[,2]) %o% (x[1,]+x[2,]) / sum(x)
        chisq <- sum((x-e)^2/e)
        p <- 1-pchisq(chisq,1)
        c(chisq=chisq, p.chisq=p)
      }))
chisq.table

## expected counts under independence for each table
E <- apply(salk2, 3,
           function(x) {
             (x[,1]+x[,2]) %o% (x[1,]+x[2,]) / sum(x)
           })
dimnames(E) <- NULL
dim(E) <- dim(salk2)
dimnames(E) <- dimnames(salk2)
E

## mh chi-square for each table (hypergeometric assumption)
apply(salk2, 3,
      function(x) {
        e <- (x[,1]+x[,2]) %o% (x[1,]+x[2,]) / sum(x)
        v <- prod(x[,1]+x[,2], x[1,]+x[2,]) / (sum(x)^2 * (sum(x)-1))
        (x-e)[1,1]^2 / v
      })


## Mantel-Haenszel chi-square components for each table
## (hypergeometric assumption)
mh.c <-
t(apply(salk2, 3,
      function(x) {
        e <- (x[,1]+x[,2]) %o% (x[1,]+x[2,]) / sum(x)
        v <- prod(x[,1]+x[,2], x[1,]+x[2,]) / (sum(x)^2 * (sum(x)-1))
        c(O=x[1,1], E=e[1,1], O.E=(x-e)[1,1], v=v, n=sum(x),
          dev=(x-e)[1,1]/sqrt(v), mh=(x-e)[1,1]^2 / v)
      }))
mh.c

## Cochran-Mantel-Haenszel test statistic
sum(mh.c[,"O.E"])^2 / sum(mh.c[,"v"])
## ')


###################################################
### code chunk number 35: twtb.tex:1974-1987
###################################################
## hhpdf("salk-dev.pdf", width=7.5, height=4)
ages <- ordered(dimnames(mh.c)[[1]], levels=dimnames(mh.c)[[1]])
barchart(mh.c[,"dev"] ~ ages, origin=0, horizontal=FALSE,
         xlab="Age Group", xlab.top="Number of Observations",
         scales=list(cex=1), ylab="standardized table deviations",
         col=BlyCol[2],
         border=BlyCol[2],
         par.settings=list(clip=list(panel=FALSE)),
         panel=function(...) {
           panel.barchart(...)
           panel.axis("top", labels=mh.c[,"n"], outside=TRUE, ticks=FALSE, rot=0)
         })
## hhdev.off()


###################################################
### code chunk number 36: twtb.tex:2043-2052
###################################################
## hhcapture("salkFisher.Rout", '
data(salk)
salk2 <- tapply(salk$Freq, salk[c(2,3,1)], c)
class(salk2) <- "table"
## salk2  ## salk2 is structured as a set of 2x2 tables
lt <- apply(salk2, 3, fisher.test, alternative="less")
## odds ratio and p-value
sapply(lt, `[`, c("estimate","p.value"))
## ')


###################################################
### code chunk number 37: twtb.tex:2105-2112
###################################################
## hhpdf("AEdotplot.pdf", width=6, height=6)
data(AEdata)
head(AEdata)
AEdotplot(AE ~ nAE/nTRT, groups = TRT, data = AEdata,
          col.AB=likertColor(2, colorFunctionArgs=list(h = c(260, 0), c = 150, l = c(30, 90))),
          panel.widths=c(.7, .3, 0))
## hhdev.off()


###################################################
### code chunk number 38: twtb.tex:2269-2293
###################################################
LikCol5 <- likertColor(5, colorFunctionOption="default")
## hhpdf("ProfChal.pdf", width=8.5, height=11)
data(ProfChal)

levels(ProfChal$Question)[4] <- "Federal, state, and local\ngovernment" ## insert line break
levels(ProfChal$Question)[5] <- "Private consultant\nself-employed" ## insert line break
levels(ProfChal$Question)[6] <- "Other (including retired,\nstudents, not employed, etc.)" ## insert line break
attributes(ProfChal)$names.dimnames <- c("Characteristic","ResponseLevel")

likert(Question ~ . | Subtable, ProfChal,
       as.percent=TRUE,    ## implies display Row Count Totals
       col=LikCol5,
       ylab=NULL, xlab=list("Percent", cex=1.5),
       main=list("Is your job professionally challenging?", x=unit(.65, "npc"), cex=1.5),
       strip.left=strip.custom(bg="gray95"),
       strip=FALSE,
       par.strip.text=list(cex=1.1, lines=5),
       auto.key=list(cex=1, size=2, between=.2),
       positive.order=TRUE,
       layout=c(1,6),
       par.settings=list(layout.widths=list(axis.key.padding=3)),
       scales=list(y=list(relation="free", cex=1.2), cex=1.2)  ## implies resizePanels
       )
## hhdev.off()


###################################################
### code chunk number 39: twtb.tex:2416-2425
###################################################
## hhpdf("PC2C.pdf", width=9.6, height=5)
EmpRows <- ProfChal$Subtable == "Employment sector"
likert(Question ~ . , data=ProfChal[EmpRows,],
       ylab=NULL, xlab=list("Count", cex=1.5),
       scales=list(x=list(cex=1.2), y=list(cex=1.2)),
       auto.key=list(cex=1, size=3, between=.2),
       col=LikCol5,
       main=list("Is your job professionally challenging?", cex=1.5))
## hhdev.off()


###################################################
### code chunk number 40: twtb.tex:2441-2450
###################################################
## hhpdf("PC2Cpct.pdf", width=9.6, height=5)
likert(Question ~ . , data=ProfChal[EmpRows,],
       as.percent=TRUE,
       ylab=NULL, xlab=list("Percent", cex=1.5),
       scales=list(x=list(cex=1.2), y=list(cex=1.2)),
       auto.key=list(cex=1, size=3, between=.2),
       col=LikCol5,
       main=list("Is your job professionally challenging?", cex=1.5))
## hhdev.off()


###################################################
### code chunk number 41: twtb.tex:2464-2474
###################################################
## hhpdf("PC2Cpctpo.pdf", width=9.6, height=5)
likert(Question ~ . , data=ProfChal[EmpRows,],
       as.percent=TRUE,
       ylab=NULL, xlab=list("Percent", cex=1.5),
       scales=list(x=list(cex=1.2), y=list(cex=1.2)),
       auto.key=list(cex=1, size=3, between=.2),
       col=LikCol5,
       main=list("Is your job professionally challenging?", cex=1.5),
       positive.order=TRUE)
## hhdev.off()


###################################################
### code chunk number 42: twtb.tex:2533-2562
###################################################
data(NZScienceTeaching)

## insert line breaks in long lines
levels(NZScienceTeaching$Question) <- c(
 "Science is really interesting",
 "I'm glad I decided to\ntake science subjects this year",
 "High school has increased\nmy interest in science",
 "I was interested in science even\nbefore I started high school",
 "My science classes are\noften taught in a boring way",
 "Science is mostly just about\nlearning facts",
 "I feel overwhelmed by\nall the options",
 "I wish there were more people\nI could talk to",
 "I worry that I'm not making\ngood choices",
 "I've been given advice or\ninformation that wasn't helpful"
)

## hhpdf("NZscienceteaching.pdf", width=8, height=6)
likert(Question ~ . | Subtable, NZScienceTeaching,
       main="New Zealand Students Still Taking Science in Year 13",
       layout=c(1,2),
       xlab="Percent",
       ylab=NULL,
       col=LikCol5,
       scales=list(y=list(relation="free")),
       strip=FALSE,
       strip.left=strip.custom(bg="gray97"),
       par.strip.text=list(cex=1.2, lines=1.2)
       )
## hhdev.off()


###################################################
### code chunk number 43: twtb.tex:2591-2618
###################################################
## hhpdf("PL5.pdf", width=17, height=9)
LikCol2 <- likertColor(2, colorFunctionOption="default")
data(USAge.table) ## from package:latticeExtra
tmp.df <- as.likertDataFrame(USAge.table[1:75, 2:1, seq(40,80,10)]/1000000)
names(tmp.df)[3] <- "Age"
names(tmp.df)[4] <- "Year"
likert(Age ~ . | Year, tmp.df,
       col=LikCol2,
       main=list("Population of United States (ages 0-74)", cex=2.5),
       xlab=list("Count in Millions\nLook for the Baby Boom\n", cex=2.5),
       ylab=list("Age", cex=2),
       ## sub=list("\nLook for the Baby Boom", cex=2.5),
       reverse=FALSE,
       scales=list(cex=1.7,
         y=list(
           limits=c(0,77),
           at=seq(1,76,5),
           labels=seq(0,75,5),
           alternating=3,
           tck=1),
         x=list(alternating=FALSE, at=-2:2)),
       auto.key=list(title=NULL, cex=2, size=4),
       strip=strip.custom(bg="gray97"),
       par.strip.text=list(cex=2.5),
       layout=c(5,1), between=list(x=.5),
       yscale.components=yscale.components.default)
## hhdev.off()


