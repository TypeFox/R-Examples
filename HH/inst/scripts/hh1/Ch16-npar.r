## vocab2.s
## har1.s
## har2.s
## balance.s
## pulse.s
##  -rwx------+   1 rmh None    31080 2004-06-08 19:45 npar/npar.tex


## vocab2.s
## vocab2.r
## continue from iinf/code/vocab.s
## continue from Ch05-iinf.r
data(vocab)

## trellis.device(color=FALSE)
## dark gray for printing the book

old.bar.fill <- trellis.par.get("bar.fill")
trellis.par.set("bar.fill", list(col=75))

stem(vocab$score-10)

tmp <- table(sign(vocab$score-10))
tmp

1-pbinom(q    = tmp["1"]-1,
         size = sum(tmp[c("-1","1")]),
         prob = .5)

tmp <- dbinom(0:50, 50, .5)
names(tmp) <- 0:50
if.R(s=
     t(barchart( ~ tmp, ref=0,
                scales=list(
                  x=list(axs="s"),
                  y=list(at=seq(1,51,5), labels=seq(0,50,5), tick=TRUE),
                  cex=1.4)))
     ,r=
     barchart(tmp ~ factor(0:50), origin=0,
              horizontal=FALSE,
              ylab="",
              scales=list(
                y=list(axs="s"),
                x=list(at=seq(1,51,5), labels=seq(0,50,5), tick=TRUE),
                cex=1.4))
     )

if.R(s=
print(position=c(0,.3,1,.7), more=TRUE,
t(barchart( ~ tmp, ref=0,
           panel=function(...) {
             panel.barchartt(...)
             axis(1, at=50, labels=FALSE, cex=1.4, tick=TRUE)
             axis(1, at=50, labels=49, cex=1.4, line=-.5, tick=FALSE)
           },
           scales=list(
             x=list(axs="s"),
             y=list(at=seq(1,51,5), labels=seq(0,50,5), tick=TRUE),
             cex=1.4),
           main="a. Binomial(n=50, p=.5)"))
)
,r=
barchart(tmp ~ factor(0:50), origin=0,
         horizontal=FALSE,
         ylab="",
         scales=list(
           y=list(axs="s"),
           x=list(at=seq(1,51,5), labels=seq(0,50,5), tick=TRUE),
           cex=1.4),
         main="a. Binomial(n=50, p=.5)")
     )

if.R(s=
print(position=c(0,0,1,.3), more=FALSE,
t(barchart( ~ tmp, ref=0,
           panel=function(...) {
             panel.barchartt(...)
             axis(1, at=50, labels=FALSE, cex=1.4, tick=TRUE)
             axis(1, at=50, labels=49, cex=1.4, line=-.5, tick=FALSE)
           },
           scales=list(
             x=list(axs="s", limits=c(0,.6e-13), tick.number=2),
             y=list(at=seq(1,51,5), labels=seq(0,50,5), tick=TRUE),
             cex=1.4),
           main="b. magnified Binomial(n=50, p=.5) to emphasize x=49 and x=50"))
)
,r=
     barchart(tmp ~ factor(0:50), origin=0,
              horizontal=FALSE,
              ylab="",
              scales=list(
                y=list(axs="s", limits=c(0,.6e-13), tick.number=2),
                x=list(at=seq(1,51,5), labels=seq(0,50,5), tick=TRUE),
                cex=1.4),
              main="b. magnified Binomial(n=50, p=.5)\nto emphasize x=49 and x=50")
)
## export.eps(hh("npar/figure/vocab.sign.eps"))

if.R(s=
t(barchart( ~ tmp,
           scales=list(
             x=list(axs="s", log=10),
             y=list(at=seq(1,51,5), labels=seq(0,50,5), tick=TRUE),
             cex=1.4)))
,r=
barchart(tmp ~ factor(0:50),
         ylab="",
         horizontal=FALSE,
           scales=list(
             y=list(axs="s", log=10),
             x=list(at=seq(1,51,5), labels=seq(0,50,5), tick=TRUE),
             cex=1.4))
 )

trellis.par.set("bar.fill", old.bar.fill)  ## restore original color


## har1.s
data(har1)

stem(har1$Post)

## sign test
tmp <- table(sign(har1$Post - har1$Pre))
tmp

pbinom(q    = tmp["1"],
       size = sum(tmp[c("-1","1")]),
       prob = .5)


## Wilcoxon signed-rank test
wilcox.test(har1$Pre, har1$Post, alternative="greater", paired=TRUE, exact=FALSE)



## har2.s
## follows har1.s

## manual construction of Wilcoxon signed-rank test

har <- data.frame(diff=har1$Pre - har1$Post)
har$abs <- abs(har$diff)
har$abs[har$abs==0] <- NA
har$rank <- rank(har$abs)
har$rank[har$diff == 0] <- 0
har$prnk <- har$rank  ## rank for positive differences
har$prnk[har$diff < 0] <- 0

old.opt <- options(length=13)  ## control titles for printing
har[order(har$abs),]           ## manually edit into columns
options(old.opt)


## calculate the statistic
har.wilk <- sum(har$prnk)
n <- sum(har$diff != 0)


## illustrate the statistics
xyplot(rank ~ diff, data=har) ## simple plot

## control font sizes
xyplot(rank ~ diff, data=har,
       scales=list(cex=1.4),
       xlab=list(cex=1.5), ylab=list(cex=1.5))
## export.eps(hh("npar/figure/har1.rank.diff.eps"))



## balance.s
data(balance)

## look at a few cases and the variable names

balance[c(1,9,10,17),]

if.R(s=
t(bwplot(age ~ sway, data=balance,
         cex=2,
         xlab=list(cex=2), ylab=list("age", cex=2),
         scales=list(cex=2))
  )
,r=
bwplot(sway ~ age, data=balance,
         cex=2,
         ylab=list(cex=2), xlab=list("age", cex=2),
         scales=list(cex=2))
)

wilcox.test(balance$sway[balance$age=="old"],
            balance$sway[balance$age=="young"],
            alternative="greater", exact=FALSE)



## pulse.s
data(pulse)
if.R(s=
print(position=c(0,0,1,.6),
t(dotplot(task ~ pulse, data=pulse,
          scales=list(cex=1.4),
          xlab=list(cex=1.4), ylab=list("task", cex=1.4)))
)
,r=
dotplot(pulse ~ task, data=pulse,
        horizontal=FALSE,
        scales=list(cex=1.4),
          ylab=list(cex=1.4), xlab=list("task", cex=1.4))
)
## export.eps(hh("npar/figure/pulse.eps"))

kruskal.test(pulse$pulse, pulse$task)
