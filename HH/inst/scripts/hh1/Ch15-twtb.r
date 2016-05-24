## Ch15-twtb.r
## drunk.s
## glasses.s
## glasses.exact.s
## blyth.s
## hypothermia.s
## -rwx------+   1 rmh None    56332 2004-05-25 01:03 twtb/twtb.tex
## salk-mh.s
## salk.s
## -rwx------+   1 rmh None     3818 2004-04-03 23:07 twtb/twtb-cmh.tex



## drunk.s
data(drunk)
prop.female <- drunk["females",]/apply(drunk,2,sum)
ages <- ordered(dimnames(drunk)$age, levels=dimnames(drunk)$age)
if.R(s=
     print(position=c(.2,0,.8,.4),
           t(barchart(ages ~ prop.female, ref=0, xlim=c(0,.13),
                      xlab="", main="proportion female",
                      col=55))
           )
     ,r=
     barchart(as.numeric(prop.female) ~ ages,
              horizontal=FALSE, origin=0,
              ylab="", main="proportion female",
              col=55)
     )
## export.eps(hh("twtb/figure/drunk.prop.fem.eps"))

drunk.er <- apply(drunk,1,sum)
drunk.ec <- apply(drunk,2,sum)
drunk.n <- sum(drunk.er)
drunk.e <- outer(drunk.er, drunk.ec) / drunk.n

chi2.ij <- (drunk -drunk.e)^2/drunk.e
chi2 <- sum(chi2.ij)

chi1.ij <- sqrt(chi2.ij) * sign(drunk -drunk.e)
chi1.ij
if.R(s=
     chi1.ij.vec <- as.vector(chi1.ij)
     ,r=
     chi1.ij.vec <- as.vector(data.matrix(chi1.ij))
     )
chi1 <- data.frame(chi=chi1.ij.vec,
                   sex=rep(dimnames(drunk)[[1]], 5),
                   age=rep(ages, rep(2,5)))
if.R(s=
print(position=c(.2,0,.8,.7),
t(barchart(age ~ chi | sex, data=chi1, ref=0,
           par.strip.text=list(cex=1.2),
           scales=list(cex=1.2, alternating=FALSE),
           layout=c(1,2),
           between=list(y=1),
           xlab=list(cex=1.2),
           ylab=list("age", cex=1.2),
           main=list("chi deviations for drunk.dat", cex=1.4),
           col=55))
)
,r=
barchart(chi ~ age | sex, data=chi1, origin=0,
         horizontal=FALSE,
         par.strip.text=list(cex=1.2),
         scales=list(cex=1.2, alternating=FALSE),
         layout=c(1,2),
         between=list(y=1),
         ylab=list(cex=1.2),
         xlab=list(cex=1.2),
         main=list("chi deviations for drunk.dat", cex=1.4),
         col=55)
)
## export.eps(hh("twtb/figure/drunk.chi.eps"))

chisq.test(data.matrix(drunk))  ## warning from cell [2,5] having expected value less than 5.




## glasses.s
data(glasses)

## draw the iso-odds ratio plot with 50% CI and 95% CI,
plotOddsRatio(glasses)
## export.eps(hh("twtb/figure/glasses.ci.eps"))

## display the 95% CI and supporting values
OddsRatio(glasses)


## calculate values in text
tmp <- log(20) + c(-1,1) * 1.96 * sqrt(1/1 + 1/5 + 1/8 + 1/2)
exp(tmp)
c(1.416116, 282.462689) / 8
2.5/3.5
c(0.1770145, 35.3078361)/(1+c(0.1770145, 35.3078361))




## glasses.exact.s

fisher.test(glasses)


## study the construction of the Fisher Exact test

## construct all possible two-way tables with the same margins as the
## initial table
all.tables <- function(x) {

  xx <- x

  result <- list()
  r.margin <- apply(x,1,sum)
  c.margin <- apply(x,2,sum)

  for (x11 in 0:r.margin[1]) {
    xx[1,1] <- x11
    xx[1,2] <- r.margin[1] - xx[1,1]
    xx[2,1] <- c.margin[1] - xx[1,1]
    xx[2,2] <- sum(x) - (xx[1,1] + xx[1,2] + xx[2,1])
    if (min(xx) >= 0) result[[as.character(x11)]] <- xx
  }
  result
}

glasses.all <- all.tables(glasses)
glasses.all


## Format the glasses.all table in LaTeX as it will appear in the book.
## We did additional editing of the file after it was automatically constructed.
library("Hmisc")  ## access latex() function
glasses.latex <-
  latex(do.call("cbind", glasses.all),
        cgroup=names(glasses.all),
        rowlabel="glasses",
        file="glasses.exact.tex")
if.R(s=
     detach("Hmisc")
     ,r=
     detach("package:Hmisc")
     )
## format end

## find hypergeometric probability for each table
g.p <-
  do.call("rbind",
          lapply(glasses.all,
                 function(x)
                 c(prob=dhyper(x[1,1], sum(x[1,]), sum(x[2,]), sum(x[,1])),
                   min=min(x))))
g.p

g.p2 <- if.R(s=
             data.frame(g.p, which="", stringsAsFactors=FALSE)
             ,r=
             data.frame(g.p, which=I(""))
             )

## intitial table
g.p2[as.character(glasses[1,1]),"which"] <- "*"
## more extreme tables (min value is smaller)
g.p2[g.p2[as.character(glasses[1,1]),"min"] > g.p2[,"min"],"which"] <- "<"
g.p2

g.p2$cumsum <- cumsum(g.p2$prob)
g.p2$rev.cumsum <- rev(cumsum(rev(g.p2$prob)))
g.p2

if.R(s=
     print(position=c(.1,0,.9,.4),
           t(barchart(0:6 ~ g.p2$prob, col=55))
           )
     ,r=
     barchart(g.p2$prob ~ factor(0:6), col=55, horizontal=FALSE)
     )

if.R(s=
print(position=c(0,0,1,.4),
      t(barchart(0:6 ~ g.p2$prob, ref=0,
                 col=55,
                 xlab="",
                 main="probability of table with specified [1,1] position",
                 panel=function(x, y, ...) {
                   panel.barchartt(x, y, ...)
                   axis(1, line=1.5, at=x, label=g.p2$which, ticks=FALSE)
               },
                 key=list(
                   text=c("observed","more extreme"),
                   text=c("*","<"),
                   border=TRUE,
                   space="right")
                 ))
)
     ,r=
     barchart(g.p2$prob ~ factor(0:6), col=55, horizontal=FALSE,
              ylab=NULL, origin=0,
              main=list(labels=
                "probability of table with specified [1,1] position"),
              scales=list(x=list(at=0:6, labels=paste(0:6,g.p2$which))),
              key=list(
                text=list(c("observed","more extreme")),
                text=list(c("*","<")),
                border=TRUE,
                space="right")
              )
     )
## export.eps(hh("twtb/figure/glasses.exact.eps"))




## blyth.s
data(blyth)

## append the location margin
blyth.t <- cbind(blyth, standard=blyth[,1]+blyth[,3], new=blyth[,2]+blyth[,4])
blyth.t

## proportion
print(format(
   round(blyth.t[2,] / apply(blyth.t, 2, sum), 2)
 ), quote=FALSE)


{
## count
blyth.t.tmp <- blyth.t
## simplify the column names by suppressing the location
dimnames(blyth.t.tmp)[[2]] <- dimnames(blyth.t.tmp)[[2]][c(5,6,5,6,5,6)]
tmp <-
  barplot(blyth.t.tmp,
          ##angle=c(45, 135),
          col=c(55,80),
          legend=dimnames(blyth.t)[[1]],
          names=as.character(dimnames(blyth.t.tmp)[[2]]),
          ylab="count", ylim=c(0,13000),
          main="Blyth's example of Simpson's paradox")
abline(v=sum(tmp[4:5])/2)
mtext(side=1, line=2, at=sum(tmp[1:2])/2, "location A", adj=.5)
mtext(side=1, line=2, at=sum(tmp[3:4])/2, "location B", adj=.5)
mtext(side=1, line=2, at=sum(tmp[5:6])/2, "combined locations", adj=.5)
## export.eps(hh("twtb/figure/blyth.count.eps"))
}

{
## proportion
blyth.t.tmp.p <- blyth.t.tmp[2,] / apply(blyth.t.tmp, 2, sum)
blyth.t.tmp.p
tmp <-
barplot(blyth.t.tmp.p,
        ##angle=135,
        col=55,
        legend=dimnames(blyth.t)[[1]][2],
        names=as.character(dimnames(blyth.t.tmp)[[2]]),
        ylab="proportion surviving",  ylim=c(0,1),
        main="Blyth's example of Simpson's paradox")
abline(v=sum(tmp[4:5])/2)
mtext(side=1, line=2, at=sum(tmp[1:2])/2, "location A", adj=.5)
mtext(side=1, line=2, at=sum(tmp[3:4])/2, "location B", adj=.5)
mtext(side=1, line=2, at=sum(tmp[5:6])/2, "combined locations", adj=.5)
mtext(side=3, at=tmp, as.character(apply(blyth.t.tmp, 2, sum)))
mtext(side=3, at=0, "count", adj=1)
## export.eps(hh("twtb/figure/blyth.proportion.eps"))
}

{
## odds
blyth.t.tmp.o <- blyth.t.tmp[2,] / blyth.t.tmp[1,]
tmp <-
barplot(blyth.t.tmp.o,
        ##angle=135,
        col=55,
        legend=dimnames(blyth.t)[[1]][2],
        names=as.character(dimnames(blyth.t.tmp)[[2]]),
        ylab="odds in favor of surviving", ylim=c(0,20),
        main="Blyth's example of Simpson's paradox")
abline(v=sum(tmp[4:5])/2)
mtext(side=1, line=2, at=sum(tmp[1:2])/2, "location A", adj=.5)
mtext(side=1, line=2, at=sum(tmp[3:4])/2, "location B", adj=.5)
mtext(side=1, line=2, at=sum(tmp[5:6])/2, "combined locations", adj=.5)
mtext(side=3, at=tmp, as.character(apply(blyth.t.tmp, 2, sum)))
mtext(side=3, at=0, "count", adj=1)
## export.eps(hh("twtb/figure/blyth.odds.eps"))
}

{
## logit
blyth.t.tmp.l <- log(blyth.t.tmp.o)
tmp <-
barplot(blyth.t.tmp.l,
        ##angle=135,
        col=55,
        legend=dimnames(blyth.t)[[1]][2],
        names=as.character(dimnames(blyth.t.tmp)[[2]]),
        ylab="logit in favor of surviving",
        main="Blyth's example of Simpson's paradox")
abline(v=sum(tmp[4:5])/2)
mtext(side=1, line=2, at=sum(tmp[1:2])/2, "location A", adj=.5)
mtext(side=1, line=2, at=sum(tmp[3:4])/2, "location B", adj=.5)
mtext(side=1, line=2, at=sum(tmp[5:6])/2, "combined locations", adj=.5)
mtext(side=3, at=tmp, as.character(apply(blyth.t.tmp, 2, sum)))
mtext(side=3, at=0, "count", adj=1)
## export.eps(hh("twtb/figure/blyth.logit.eps"))
}




## hypothermia.s
hypothermia <- matrix(c(75,54,61,83),
                      nrow=2,
                      dimnames=list(
                        c("treated","control"),
                        c("favorable","not.favorable")))

plotOddsRatio(hypothermia)
## export.eps(hh("twtb/figure/hypothermia.ci.eps"))

## count
barplot(t(hypothermia),
        ylim=c(0,160), ylab="count",
        legend=dimnames(hypothermia)[[2]],
        names=as.character(dimnames(hypothermia)[[1]]))

## proportion
barplot(hypothermia[,1] / apply(hypothermia, 1, sum),
        ylim=c(0,.6), ylab="proportion favorable",
        legend=dimnames(hypothermia)[[2]][1],
        names=as.character(dimnames(hypothermia)[[1]]))

## odds
barplot(hypothermia[,1] / hypothermia[,2],
        ylim=c(0,1.4), ylab="odds favorable",
        legend=dimnames(hypothermia)[[2]][1],
        names=as.character(dimnames(hypothermia)[[1]]))

## logit
barplot(log(hypothermia[,1] / hypothermia[,2]),
        ylim=c(-.6,.3), ylab="logit favorable",
        legend=dimnames(hypothermia)[[2]][1],
        names=as.character(dimnames(hypothermia)[[1]]))




## salk.s
data(salk)
if.R(r={
  require(vcd)
  structable(~ age + vaccine + paralyze, data=salk,  direction=c("v","v","h"),
             main="Recommended display---specified, paralyze split last")
  detach("package:vcd")
}, s={
    tmp <- tapply(salk$Freq, salk[, c(2,3,1)], c)
    names(dimnames(tmp)) <- c("vaccine","paralyze", "age")
    class(tmp) <- "table"
    print(tmp)
    rm(tmp)
})

salk2 <- tapply(salk$Freq, salk[c(2,3,1)], c)
if.R(r=class(salk2) <- "table",
     s={})

salk2
mantelhaen.test(salk2)
mantelhaen.test(salk2, correct=FALSE)


## salk-mh.s


## Use latex() to generate the table in latex format.
## Then lots of editing to smooth out the appearance.
library("Hmisc")
for (i in 1:6) {
  tmp <- salk2[,,i]
  latex(tmp, file=paste("salk",i,"tex", sep="."),
        caption=dimnames(salk2)[[3]][i],
        label=dimnames(salk2)[[3]][i])
}
if.R(s=
     detach("Hmisc")
     ,r=
     detach("package:Hmisc")
     )
##


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
## latex
library("Hmisc")
for (i in 1:6) {
  tmp <- E[,,i]
  latex(tmp, file=paste("E",i,"tex", sep="."),
        caption=dimnames(E)[[3]][i],
        label=dimnames(E)[[3]][i])
}
if.R(s=
     detach("Hmisc")
     ,r=
     detach("package:Hmisc")
     )
##


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

## plot the dev values for each table.
if.R(r={
  ages <- ordered(dimnames(mh.c)[[1]], levels=dimnames(mh.c)[[1]])
  barchart(mh.c[,"dev"] ~ ages, origin=0, horizontal=FALSE,
           scales=list(y=list(axs="s"), cex=1.4),
           ylab="standardized table deviations", col=55)
}, s=
     print(position=c(0,0,1,.7),
           t(barchart( ~ mh.c[,"dev"], ref=0,
                      scales=list(x=list(axs="s"), cex=1.4),
                      xlab="", col=55,
                      panel=function(x,y,...) {
                        panel.barchartt(x,y,...)
                        axis(3, at=x, labels=mh.c[,"n"],
                             ticks=FALSE, cex=1.4)
                        axis(3, at=.3, labels="count",
                             ticks=FALSE, xpd=TRUE, adj=1, cex=1.4)
                        mtext(side=2, "standardized table deviations",
                              cex=1.4, line=4.5)
                      }
                      ))
           )
     )
## export.eps(hh("twtb/figure/salk.dev.eps"))



## prepare to print the table with LaTeX
abind::abind(salk2, E, along=2)

tmp <-
abind::abind(salk2, round(E,2), "p(no.par)"=round(t(pp),3),
      abind::abind(round(t(chisq.table),3),
            array(NA, dim=dim(t(chisq.table))),
            along=.1),
      abind::abind(round(t(mh.c),2), array(NA, dim=dim(t(mh.c))), along=.1),
      along=2)
tmp


tmp4 <- aperm(tmp, c(1,3,2))
dn <- dimnames(tmp4)
dimnames(tmp4) <- NULL
dim(tmp4) <- c(prod(dim(tmp4)[1:2]), dim(tmp4)[3])
dimnames(tmp4) <- list( as.vector(t(outer(dn[[2]], dn[[1]], paste))), dn[[3]] )
dimnames(tmp4)[[1]] <- rep(dn[[1]],6)
tmp4

##
library("Hmisc")
salk.latex <-
latex(tmp4, rgroup=dn[[2]], n.rgroup=rep(2,6), rowlabel="", na.blank=TRUE,
      cdec=c(0,0,2,2,3,3,3,0,2,2,2,0,2,2),
      cgroup=c("Observed","Expected","","Chi-Square",
        "[1,1] position for Mantel--Haenszel test"),
      n.cgroup=c(2,2,1,2,7),
      file="salk.tex")
if.R(s=
     detach("Hmisc")
     ,r=
     detach("package:Hmisc")
     )
## Edit the latex output into Table
