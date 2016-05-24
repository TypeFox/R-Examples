### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/npar.tex'

###################################################
### code chunk number 1: npar.tex:9-10
###################################################
library(HH)


###################################################
### code chunk number 2: npar.tex:13-18
###################################################
## the standard lattice color 2 is difficult for people with color deficient vision
data(col3x2)
## These colors look like a 3x2 color array when run through
## the vischeck simulator to see how they look for the three most
## common color vision deficiencies: Protanope, Deuteranope, Tritanope.


###################################################
### code chunk number 3: npar.tex:123-131
###################################################
## hhcapture("vocabnpar.Rout", '
data(vocab)
table(vocab)
table(vocab$score - 10)
table( sign(vocab$score-10) )
pbinom(48.5, 50, .5, lower=FALSE)
1 - pbinom(48.5, 50, .5)
## ')


###################################################
### code chunk number 4: npar.tex:174-214
###################################################
## hhpdf("vocab-sign.pdf", width=7, height=7.5)
tmp <- dbinom(0:50, 50, .5)
names(tmp) <- 0:50

colnpar <- likertColor(2)[2]

print(more=TRUE, position=c(0,.67,1,1),
barchart(tmp ~ factor(0:50), origin=0,
         horizontal=FALSE,
         ylab="",
         border=colnpar, col=colnpar,
         scales=list(
           y=list(axs="s", at=seq(0,.1,.05)),
           x=list(at=seq(1,51,5), labels=seq(0,50,5), tick=TRUE),
           cex=1.2),
         main="a. Binomial(n=50, p=.5)")
)
print(more=TRUE, position=c(0,.335,1,.665),
barchart(tmp ~ factor(0:50),
         horizontal=FALSE,
         ylab="",
         border=colnpar, col=colnpar,
           scales=list(
             y=list(axs="s", log=10, at=10^-c(1,5,10,15)),
             x=list(at=seq(1,51,5), labels=seq(0,50,5), tick=TRUE),
             cex=1.2),
         main="b. Binomial(n=50, p=.5) on log scale")
)
print(more=FALSE, position=c(0,0,1,.33),
barchart(tmp ~ factor(0:50), origin=0,
         horizontal=FALSE,
         ylab="",
         border=colnpar, col=colnpar,
         scales=list(
           y=list(axs="s", limits=c(-.01e-13,.6e-13), at=c(0,3,6)*1e-14),
           x=list(at=seq(1,51,5), labels=seq(0,50,5), tick=TRUE),
           cex=1.2),
         main="c. magnified Binomial(n=50, p=.5)\nto emphasize x=49 and x=50")
)
## hhdev.off()


###################################################
### code chunk number 5: npar.tex:306-324
###################################################
## hhpdf("har1.pdf", width=7, height=5)
data(har1)
LL <-
histogram( ~ Post + Pre + (Post-Pre), data=har1, layout=c(1,3), breaks=seq(-16,6,1),
          between=list(y=1), xlab=NULL, scales=list(alternating=1),
          col=colnpar,
          xlim=c(-17,10)) +
  layer(panel.abline(v=0, lty=2, col="gray40"))
RR <-
densityplot( ~ Post + Pre + (Post-Pre), data=har1, layout=c(1,3), breaks=seq(-16,6,1),
            between=list(y=1), xlab=NULL, outer=TRUE, scales=list(alternating=1),
            col=colnpar,
            xlim=c(-17,10)) +
  layer(panel.abline(v=0, lty=2, col="gray40"))
LLRR <- matrix.trellis(c(LL, RR), byrow=TRUE, nrow=2, ncol=3)
dimnames(LLRR) <- list((dimnames(LL)[[1]]), c(LL$ylab.default, RR$ylab.default))
useOuterStrips(update(LLRR, between=list(x=1.5, y=2), ylab=NULL))
## hhdev.off()


###################################################
### code chunk number 6: npar.tex:381-385
###################################################
## hhcapture("har1a.Rout", '
table(sign(har1$Post - har1$Pre))
pbinom(15, 46, .5, lower=TRUE)
## ')


###################################################
### code chunk number 7: npar.tex:441-451
###################################################
## hhcapture("har2.Rout", '
har <- data.frame(diff=har1$Pre - har1$Post)
har$abs <- abs(har$diff)
har$abs[har$abs==0] <- NA
har$rank <- rank(har$abs)
har$rank[har$diff == 0] <- 0
har$prnk <- har$rank    ## rank for positive differences
har$prnk[har$diff < 0] <- 0
har[order(har$abs),]    ## manually edit into columns
## ')


###################################################
### code chunk number 8: npar.tex:480-501
###################################################
## hhcapture("har2b.Rout", '
## calculate the statistic
sum.prnk <- sum(har$prnk)
n <- sum(har$diff != 0)
sum.prnk
n

## normal approximation
mean.prnk <- n*(n+1)/4
numerator <- sum.prnk - mean.prnk
var.prnk1 <- n * (n + 1) * (2 * n + 1)/24

NTIES <- table(har$abs[1:46]) ## non-zero differences
TIEadjustment <- sum(NTIES^3 - NTIES)/48
var.afterTIES <- var.prnk1 - TIEadjustment

z <- (numerator - .5) / sqrt(var.afterTIES)
z
p.val <- pnorm(z, lower.tail = FALSE)
p.val
## ')


###################################################
### code chunk number 9: npar.tex:532-545
###################################################
## hhpdf("har1-rank-diff.pdf", width=7.5, height=7.5)
har.order <- har[order(har$abs)[1:46],]
xyplot(rank ~ jitter(diff, factor=10), data=har.order, pch=19, cex=1.5, col=colnpar,
       xlim=c(-6, 10),
       sub="Displayed numbers are the signed-ranks.\nTied ranks are shown by average rank and overstrikes.") +
  layer({rank <- har.order$rank
         diff.neg <- har.order$diff < 0
         rank[diff.neg] <- -rank[diff.neg]
         x[diff.neg] <- x[diff.neg] - 1.3
         panel.text(x+.2, y, labels=rank, adj=0)
       }) +
  layer(panel.abline(v=0, h=0, lty=2, col="gray50"))
## hhdev.off()


###################################################
### code chunk number 10: npar.tex:566-570
###################################################
## hhcapture("har1b.Rout", '
wilcox.test(har1$Pre, har1$Post, alternative="greater",
            paired=TRUE, exact=FALSE)
## ')


###################################################
### code chunk number 11: npar.tex:671-678
###################################################
## hhpdf("balance.pdf", width=5.5, height=3, col=likertColor(2)) ## col is not an argument for grDevices:::pdf
data(balance)
bwplot(sway ~ age, data=balance,
       panel=panel.bwplot.superpose, groups=age,
       par.settings=list(plot.symbol=list(pch=19)),
       xlab="age")
## hhdev.off()


###################################################
### code chunk number 12: npar.tex:705-710
###################################################
## hhcapture("balancesway.Rout", '
data(balance)
cbind(balance, rank(balance$sway))
sum(rank(balance$sway)[10:17])
## ')


###################################################
### code chunk number 13: npar.tex:741-745
###################################################
## hhpdf("balancedots.pdf", width=5.5, height=3, lwd=4) ## lwd is not an argument of grDevices::pdf
dotplot(rank(sway) ~ age, group=age, data=balance, col=likertColor(2),
       par.settings=list(superpose.symbol=list(pch=19)))
## hhdev.off()


###################################################
### code chunk number 14: npar.tex:754-759
###################################################
## hhcapture("balance.Rout", '
wilcox.test(balance$sway[balance$age=="young"],
            balance$sway[balance$age=="old"],
            alternative="less", exact=FALSE)
## ')


###################################################
### code chunk number 15: npar.tex:774-802
###################################################
## hhcapture("balancebyhand.Rout", '
data(balance)

unname(n.old <- table(balance$age)["old"])
unname(n.young <-  table(balance$age)["young"])

(r <- rank(balance$sway))
(NTIES <- table(r))
unname(tie.adjustment <-
         sum(NTIES^3 - NTIES) /
         ((n.old + n.young) * (n.old + n.young - 1)))

unname(SIGMA <-
         sqrt((n.old * n.young/12) *
              ((n.old + n.young + 1) - tie.adjustment)))

unname(STATISTIC.old <-
         c(W = sum(r[1:9]) - n.old * (n.old + 1)/2))
unname(z.old <- STATISTIC.old - n.old * n.young/2)
unname(z.old <- (z.old - .5)/SIGMA)
unname(pnorm(z.old, lower.tail = FALSE))

unname(STATISTIC.young <-
         c(W = sum(r[10:17]) - n.young * (n.young + 1)/2))
unname(z.young <- STATISTIC.young - n.old * n.young/2)
unname(z.young <- (z.young + .5)/SIGMA)
unname(pnorm(z.young))
## ')


###################################################
### code chunk number 16: npar.tex:875-882
###################################################
## hhpdf("pulse.pdf", width=7, height=3.5, lwd=4) ## lwd is not an argument of grDevices::pdf
data(pulse)
dotplot(pulse ~ task, data=pulse,
        horizontal=FALSE, col=col3x2[3],
        scales=list(cex=1.4),
          ylab=list(cex=1.4), xlab=list("task", cex=1.4))
## hhdev.off()


###################################################
### code chunk number 17: npar.tex:898-903
###################################################
old.width <- options(width=70)
## hhcapture("pulse.Rout", '
kruskal.test(pulse ~ task, data=pulse)
## ')
options(old.width)


