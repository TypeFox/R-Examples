## Ch14-dsgntwo.r
## yatesppl.s
## yatesppl-layout.s
## yatesppl-bwplot.s
## yatesppl-alt.s
## yatesppl-mmc2.s
## yatesppl.ex.s
## 2.8-2-full.s
## 2.8-2.s
## circuit.s
## cc135.s
## cc135-bwplot.s
## apple3.r
## apple3.s
## testscore.s
## crash.s
## crash-cover.s
##  -rwx------+   1 rmh None   109190 2004-06-25 18:15 dsgntwo/dsgntwo.tex

library(HH)


## yatesppl.s
## Split Plot design.  This example is from Yates 1937.
data(yatesppl)
yatesppl <- yatesppl

## split plot analysis
yatesppl.anova <- aov(y ~ variety*nitrogen +
                          Error(blocks/plots/subplots),
                      data=yatesppl)
summary(yatesppl.anova)
model.tables(yatesppl.anova, type="means", se=TRUE)


## incorrect analysis that ignores split plot
yatesppl.wrong.anova <-
  aov(y ~ (blocks*variety)+(nitrogen*variety),
      data=yatesppl)
summary(yatesppl.wrong.anova)
model.tables(yatesppl.wrong.anova, type="means", se=TRUE)


## Whole plot analysis, compare to top half of split plot

yates.whole   <-
  yatesppl[seq(1,72,4), c("blocks","plots","variety")]
yates.whole$y <- apply(matrix(yatesppl$y, 4,18), 2, sum)
yates.whole$y
yatesppl.whole.anova <- aov(y ~ blocks + variety,
                            data=yates.whole)
summary(yatesppl.whole.anova)

par(mfrow=c(2,3))
plot(y ~ blocks+variety+nitrogen, data=yatesppl, ask=FALSE)
plot(y ~ blocks + variety, data=yates.whole, ask=FALSE)
par(mfrow=c(1,1))

## polynomial contrasts in nitrogen
contrasts(yatesppl$nitrogen)
contrasts(yatesppl$nitrogen) <- contr.poly(4)
contrasts(yatesppl$nitrogen)

## split plot analysis with polynomial contrasts
yatespplp.anova <- aov(y ~ variety*nitrogen +
                          Error(blocks/plots/subplots),
                      data=yatesppl)
summary(yatespplp.anova,
        split=list(nitrogen=list(linear=1, quad=2, cub=3)))
model.tables(yatespplp.anova, type="means", se=TRUE)



## yatesppl-layout.s
## display design set up by yatesppl.s

yatesppl$blocks   <- factor(yatesppl$blocks,
                            labels=paste("B",1:6,sep=""))
yatesppl$plots    <- factor(yatesppl$plots,
                            labels=paste("P",1:3,sep=""))
yatesppl$subplots <- factor(yatesppl$subplots,
                            labels=paste("S",1:4,sep=""))
yatesppl$variety  <- factor(yatesppl$variety,
                            labels=paste("V",1:3,sep=""))
levels(yatesppl$nitrogen) <- paste("N",c(0,2,4,6),sep=".")
contrasts(yatesppl$nitrogen) <- contr.poly(4)  ## the contrasts were changed
                                               ## by the levels statement
contrasts(yatesppl$nitrogen)

tmp <-
tapply(paste(as.character(yatesppl$variety),
             as.character(yatesppl$nitrogen), sep=":"),
       yatesppl[c("subplots","plots","blocks")],
       FUN=c)
print(tmp, quote=FALSE)


## yatesppl-bwplot.s
position(yatesppl$variety) <- (1:3)+.5
yatesppl$nit.lev <- factor(yatesppl$nitlev)
interaction2wt(y ~ variety * nit.lev, data=yatesppl,
               par.strip.text=list(cex=1.4),
               scale=list(x=list(cex=1), y=list(cex=1, alternating=FALSE)),
               main.cex=1.6)
## export.eps(hh("dsgntwo/figure/yatesppl.eps"))


## yatesppl-alt.s
## statement in book and recommended way to think
## about split plot designs
yatesppl.anova <- aov(y ~ variety*nitrogen +
                          Error(blocks/plots/subplots),
                      data=yatesppl)
summary(yatesppl.anova)


## alternate specification
yatesppl2.anova <- aov(y ~ variety*nitrogen +
                          Error(blocks/variety/nitrogen),
                      data=yatesppl)
summary(yatesppl2.anova)



## yatesppl-mmc2.s
## follows yatesppl.s and yatesppl-layout.s

yatespplp.anova <- aov(y ~ variety*nitrogen +
                           Error(blocks/plots/subplots),
                       data=yatesppl)
summary(yatespplp.anova,
        split=list(nitrogen=list(linear=1, quad=2, cub=3)),
        expand.split=FALSE)
##
try.out <- if.R(r=try(glht(yatespplp.anova)),            ## glht can't handle aovlist
                s=try(multicomp(yatespplp.anova)$lmat))  ## multicomp can't handle aovlist

## incorrect analysis that ignores split plot (hence variety test is wrong),
## but gives all the right numbers for the nitrogen effect.
contrasts(yatesppl$nitrogen) <- contr.treatment(4)
contrasts(yatesppl$variety)  <- contr.treatment(3)
yatesppl.wrong.anova <-
  aov(terms(y ~ (blocks*variety)+(variety*nitrogen), keep.order=TRUE),
      data=yatesppl)
summary(yatesppl.wrong.anova)
##
if.R(r={tmp <- glht(yatesppl.wrong.anova,
                    linfct=mcp(nitrogen="Tukey",
                      `interaction_average`=TRUE, `covariate_average`=TRUE))},
     s={
       tmp <- multicomp(yatesppl.wrong.anova, focus="nitrogen")
       print(tmp)
       print(zapsmall(tmp$lmat))
       print(zapsmall(tmp$lmat[29:32,]))
       tmp
     })

if.R(r={
 ywa.mca <- glhtWithMCP.993(yatesppl.wrong.anova, linfct=mcp(nitrogen="Tukey"))
 ywa.mmc <- mmc(yatesppl.wrong.anova, linfct=ywa.mca$linfct, focus="nitrogen")
 old.omd <- par(omd=c(0,.95,0,1))
 plot(ywa.mmc, ry=c(73,127))
 par(old.omd)
 print(ywa.mmc)
}, s={
  old.mar <- par(mar=c(5,4,4,4)+.1)
  multicomp.mmc(yatesppl.wrong.anova, focus="nitrogen")
  par(old.mar)
})
## export.eps(hh("dsgntwo/figure/yatesppl-mmc-pairs.eps"))

contrasts(yatesppl$nitrogen)

## Pedagogical construction, see below for application construction of lmat
if.R(r={
  tl <- t(tmp$linfct)
  rows.nitr <- 19:21
}, s={
  tl <- tmp$lmat
  rows.nitr <- 29:32
})
yatesppl.lmat <- tl[,1:3]
dimnames(yatesppl.lmat)[[2]] <- c("lin","quad","cub")
yatesppl.lmat[,"lin"]  <- -(tl[,1]+tl[,2]+tl[,3]+tl[,4]+tl[,5]+tl[,6])
yatesppl.lmat[,"quad"] <-   tl[,1]-tl[,6]
yatesppl.lmat[,"cub"]  <-  -tl[,1]+tl[,2]-tl[,3]+tl[,4]+tl[,5]-tl[,6]
yatesppl.lmat[rows.nitr,]
yatesppl.lmat <- sweep(yatesppl.lmat, 2,
                       apply(abs(yatesppl.lmat[rows.nitr,]), 2, sum)/2, "/")
zapsmall(yatesppl.lmat)[rows.nitr,]
## End Pedagogical construction

if.R(r={
  focus.lmat <- rbind(nitrogen1=-colSums(yatesppl.lmat[rows.nitr,]),
                      yatesppl.lmat[rows.nitr,])
  row.names(focus.lmat) <- levels(yatesppl$nitrogen)
  ywa.mmc <- mmc(yatesppl.wrong.anova, focus.lmat=focus.lmat, focus="nitrogen")
  old.omd <- par(omd=c(0,.95,0,1))
  plot(ywa.mmc, ry=c(73,127))
  par(old.omd)
}, s={
  old.mar <- par(mar=c(5,4,4,4)+.1)
  ywa.mmc <- multicomp.mmc(yatesppl.wrong.anova, focus="nitrogen",
                           lmat=yatesppl.lmat, lmat.rows=rows.nitr)
  par(old.mar)
})
  print(ywa.mmc)
## export.eps(hh("dsgntwo/figure/yatesppl-mmc.eps"))


## yatesppl.ex.s
### \begin{enumerate}
### \item The whole plot column space is defined by the  \verb+plots %in% blocks+
### dummy variables generated by the
### \begin{verbatim}
## alternate residuals formula
yatesppl.resida.aov <- aov(y ~ blocks/plots,
                           data=yatesppl, x=TRUE)
summary(yatesppl.resida.aov)
t(yatesppl.resida.aov$x)
### \end{verbatim}%$
###
### This is the same column space defined by the \verb;variety + blocks:variety;
### dummy variables generated by the
### \begin{verbatim}
## computational shortcut
yatesppl.short.aov <-
  aov(terms(y ~ blocks + variety + blocks*variety +
            nitrogen + variety*nitrogen,
            keep.order=TRUE),  ## try it without keep.order=TRUE
      data=yatesppl, x=TRUE)
summary(yatesppl.short.aov)
t(yatesppl.short.aov$x)
### \end{verbatim}%$
###
### We illustrate this by regressing the response variable {\tt y} on
### the  \verb;variety + blocks:variety;
### dummy variables
### \begin{verbatim}
## project y onto blocks/plots dummy variables
plots.aov <- lm(y ~ yatesppl.resida.aov$x[,7:18], data=yatesppl)
summary.aov(plots.aov)
y.bp <- predict(plots.aov)
variety.aov <- aov(y.bp ~ blocks*variety, data=yatesppl)
summary(variety.aov)
### \end{verbatim}%$
### and seeing that we reproduce the \verb+plots %in% blocks+
### stratum of the ANOVA table
### \begin{verbatim}
### Error: plots %in% blocks
###           Df Sum of Sq  Mean Sq F Value     Pr(F)
###   variety  2  1786.361 893.1806 1.48534 0.2723869
### Residuals 10  6013.306 601.3306
### \begin{verbatim}
### obtained from the complete five-factor specification.
###
### \begin{verbatim}
## split plot analysis
yatesppl.anova <- aov(y ~ variety*nitrogen +
                          Error(blocks/plots/subplots),
                      data=yatesppl)
summary(yatesppl.anova)
### \end{verbatim}
###
### \end{enumerate}


## 2.8-2-full.s
data(R282)
R282.aov <- aov(y ~ blocks + (AA+BB+CC+DD+EE+FF+GG+HH)^2, data=R282, x=TRUE)
anova(R282.aov)
R282.aov$x
cor(R282.aov$x[,-1])  ## suppress the constant column of 1
summary.lm(R282.aov)


## 2.8-2.s
y ~ blocks + (a+b+c+d+e+f+g+h)^2



## circuit.s
## Analysis of Integrated Circuit Data
data(circuit)

circuit.aov <- aov( yield ~ A + B + C + A:B, data=circuit)

summary(circuit.aov)

tapply(circuit[,"yield"], circuit[,"A"], mean)
tapply(circuit[,"yield"], circuit[,"B"], mean)
tapply(circuit[,"yield"], circuit[,"C"], mean)
tapply(circuit[,"yield"], circuit[,c("A","B")], mean)
tapply(circuit[,"yield"], circuit[,"D"], mean)
tapply(circuit[,"yield"], circuit[,c("A","C")], mean)

interaction2wt(yield ~ A*B*C, data=circuit,
               par.strip.text=list(cex=1.4),
               scale=list(x=list(cex=1), y=list(cex=1, alternating=FALSE)),
               main.cex=1.6)
## export.eps(hh("dsgntwo/figure/circuit.eps"))



## cc135.s
#Cochran and Cox, p135.  Heiberger p 601.  HH pp449-451.
# Residual effects design

data(cc135)
cc135 <- cc135
if.R(s=
     a1c <-  aov(terms(yield ~ cow + square/period + treat + res.treat,
                       keep.order=TRUE), data=cc135)
     ,r= ## bug in R 1.9.1 model.tables
     a1c <-  aov(terms(yield ~ cow + square:period + treat + res.treat,
                       keep.order=TRUE), data=cc135)
     )
##                  terms()   keeps the order as specified
summary(a1c)

model.tables(a1c, type="means")

if.R(s=
     a1cr <- aov(terms(yield ~ cow + square/period + res.treat + treat,
                       keep.order=TRUE), data=cc135)
     ,r= ## bug in R 1.9.1 model.tables
     a1cr <- aov(terms(yield ~ cow + square:period + res.treat + treat,
                       keep.order=TRUE), data=cc135)
     )
##                   terms()   keeps the order as specified
summary(a1cr)
model.tables(a1cr, type="means")

if.R(s=
     apply(summary(a1cr)[,1:2], 2, sum)
     ,r=
     apply(summary(a1cr)[[1]][,1:2], 2, sum)
     )


## cc135-bwplot.s

## construct the yield adjusted for the blocking factors
cc135.block.aov <- aov(yield ~ cow + square/period, data=cc135)
summary(cc135.block.aov)
cc135$y.adj <- mean(cc135$yield) - resid(cc135.block.aov)

a1cr.adj <- aov(terms(y.adj ~ cow + square/period + res.treat + treat,
                  keep.order=TRUE), data=cc135)
summary(a1cr.adj)


## construct the yield adjusted for the blocking factors and res.treat
cc135.block.res.aov <- aov(yield ~ cow + square/period + res.treat, data=cc135)
summary(cc135.block.res.aov)
cc135$y.adj.res <- mean(cc135$yield) - resid(cc135.block.res.aov)

a1cr.adj.res <- aov(terms(y.adj.res ~ cow + square/period + res.treat + treat,
                  keep.order=TRUE), data=cc135)
summary(a1cr.adj.res)

print(split = c(1,1,2,1), more = TRUE,  # left
if.R(s=
t(bwplot(treat ~ y.adj.res, data=cc135, xlim=c(35,80),
         scales=list(y=list(cex=.7), x=list(cex=1.4)), xlab=list(cex=1.4),
         main="treatment means adjusted\n for blocks and residual treatments"))
,r=
bwplot(treat ~ y.adj.res, data=cc135, xlim=c(35,80),
         scales=list(y=list(cex=.7), x=list(cex=1.4)), xlab=list(cex=1.4),
         main="treatment means adjusted\n for blocks and residual treatments"))
)



## construct the yield adjusted for the blocking factors and treat
cc135.block.treat.aov <- aov(yield ~ cow + square/period + treat, data=cc135)
summary(cc135.block.treat.aov)
cc135$y.adj.treat <- mean(cc135$yield) - resid(cc135.block.treat.aov)

a1cr.adj.treat <- aov(terms(y.adj.treat ~ cow + square/period + treat + treat,
                  keep.order=TRUE), data=cc135)
summary(a1cr.adj.treat)

print(split = c(2,1,2,1), more = FALSE,  # right
if.R(s=
t(bwplot(res.treat ~ y.adj.treat, data=cc135, xlim=c(35,80),
         scales=list(y=list(cex=.7), x=list(cex=1.4)), xlab=list(cex=1.4),
         main="residual treatment means\n adjusted for blocks and treatments"))
,r=
bwplot(res.treat ~ y.adj.treat, data=cc135, xlim=c(35,80),
         scales=list(y=list(cex=.7), x=list(cex=1.4)), xlab=list(cex=1.4),
         main="residual treatment means\n adjusted for blocks and treatments"))
      )


## this size looks better on paper
print(position=c(-.1, .18, .45, .88), more = TRUE,  # left
if.R(s=
t(bwplot(treat ~ y.adj.res, data=cc135, xlim=c(35,80),
         scales=list(y=list(cex=.9), x=list(cex=1.4)), xlab=list(cex=1.4)))
,r=
bwplot(treat ~ y.adj.res, data=cc135, xlim=c(35,80),
         scales=list(y=list(cex=.9), x=list(cex=1.4)), xlab=list(cex=1.4)))
      )
mtext("treatment means adjusted\n for blocks and residual treatments",
      side=3, at=.03)
par(new=TRUE)
print(position=c(.4, .18, 1.05, .88), more = FALSE,  # right
if.R(s=
t(bwplot(res.treat ~ y.adj.treat, data=cc135, xlim=c(35,80),
         scales=list(y=list(cex=.9), x=list(cex=1.4)), xlab=list(cex=1.4)))
,r=
bwplot(res.treat ~ y.adj.treat, data=cc135, xlim=c(35,80),
         scales=list(y=list(cex=.9), x=list(cex=1.4)), xlab=list(cex=1.4)))
      )
mtext("residual treatment means\n adjusted for blocks and treatments",
      side=3, at=.82)

## export.eps(hh("dsgntwo/figure/cc135.f.bwplot.eps"))



## apple3.r
## apple3.s
cat("See file Ch14-apple3.r\n")



## testscore.s
## testscore data:
## P. O. Johnson and F. Tsao, 1945
## R. L. Anderson and T. A. Bancroft, 1952

data(testscore)
testscore <- testscore

splom( ~ testscore, cex=.5, pch=16,
      main="Original ordering of factor values")
## export.eps(hh("dsgntwo/figure/testscore.f.splom1.eps"))


## reorder the levels of the factors, and the order of the variables
## to improve the simplicity of the graph.
testscore$standing  <- ordered(testscore$standing,
                              levels=rev(c("good", "average", "poor")))
testscore$order     <- ordered(testscore$order,
                              levels=rev(c("high", "medium", "low")))
splom( ~ testscore[,c(1,6,7,4,5,3,2)], cex=.5, pch=16,
      main="Revised ordering of factor values")
## export.eps(hh("dsgntwo/figure/testscore.f.splom2.eps"))



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



## after looking at all of above
## Total Sum of Squares
if.R(s=
     var(testscore$final, S=TRUE)
     ,r=
     var(testscore$final) * (length(testscore$final)-1)
)

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

round(summary.lm(testscore7.aov)$coef, digits=4)
proj(testscore7.aov, onedf=TRUE)
testscore7.aov <- aov(final ~ initial + mental.age +
                      standing + order + sex,
                      data=testscore, x=TRUE)
testscore$final.adj <-
  testscore$final - as.vector(apply(proj(testscore7.aov, onedf=TRUE)[,1:3],1,sum))

if.R(s=
print(position=c(-.05,0,1.05,1),
      t(dotplot(standing ~ final.adj | order + sex, data=testscore,
                ylab=list("standing",cex=1.7),
                xlab=list("adjusted final scores",cex=1.7),
                scales=list(cex=1.4, alternating=FALSE),
                par.strip.text=list(cex=1.7),
                strip=strip.ladder, cex=1)
        ))
,r=
dotplot(final.adj ~ standing | order + sex, data=testscore,
                xlab=list("standing",cex=1.7),
                ylab=list("adjusted final scores",cex=1.7),
                scales=list(cex=1.4, alternating=FALSE),
                par.strip.text=list(cex=1.7),
                strip=strip.ladder, cex=1)
        )
## export.eps(hh("dsgntwo/figure/final.adjs.eps"))

if.R(s=
print(position=c(-.05,0,1.05,1),
      t(dotplot(order ~ final.adj | standing + sex, data=testscore,
                ylab=list("order",cex=1.7),
                xlab=list("adjusted final scores",cex=1.7),
                scales=list(cex=1.4, alternating=FALSE),
                par.strip.text=list(cex=1.7),
                strip=strip.ladder, cex=1)
        ))
,r=
dotplot(final.adj ~ order | standing + sex, data=testscore,
                xlab=list("order",cex=1.7),
                ylab=list("adjusted final scores",cex=1.7),
                scales=list(cex=1.4, alternating=FALSE),
                par.strip.text=list(cex=1.7),
                strip=strip.ladder, cex=1)
        )
## export.eps(hh("dsgntwo/figure/final.adjo.eps"))



interaction2wt(final.adj ~ order + standing + sex, data=testscore,
               par.strip.text=list(cex=1.),
               main.cex=1.6,
               scales=list(cex=1, y=list(alternating=1)))

interaction2wt(final.adj ~ order + standing,
               data=testscore[testscore$sex=="male",],
               ## above is necessary, below controls formatting.
               strip=FALSE,
               main="final.adj for male", xlab="",
               main.cex=1.6,
               ylim=c(-3.5,3),
               scales=list(
                 cex=1.2,
                 x=list(cex=1),
                 y=list(cex=1, alternating=1)),
               responselab="",
               key.print=FALSE,
               key.cex.title=1.4,
               key.cex.text=1)
## export.eps(hh("dsgntwo/figure/final.2wtm.eps"))


interaction2wt(final.adj ~ order + standing,
               data=testscore[testscore$sex=="female",],
               ## above is necessary, below controls formatting.
               strip=FALSE,
               main="final.adj for female", xlab="",
               main.cex=1.6,
               ylim=c(-3.5,3),
               scales=list(
                 cex=1.2,
                 x=list(cex=1),
                 y=list(cex=1, alternating=1)),
               responselab="",
               key.cex.title=1.4,
               key.cex.text=1)
## export.eps(hh("dsgntwo/figure/final.2wtf.eps"))




## prediction
newdata <- cbind(initial=mean(testscore$initial),
                 mental.age=mean(testscore$mental.age),
                 testscore[c(1:9,28:36), c("standing","order","sex")])
newdata
final.pred <- predict(testscore7.aov, newdata=newdata)
final.pred

final.pred <- tapply(final.pred, newdata[,3:5], c)
class(final.pred) <- "table"
final.pred

## now summarize this over each factor to get predicted values

apply(final.pred, 1, mean) ## For each scholastic standing
apply(final.pred, 2, mean) ## For each individual order




## crash.s
data(crash)
crash <- crash

## crash-interaction.eps.gz
par(mfrow=c(2,2))
interaction.plot(crash$agerange, crash$passengers, crash$crashrate,
                 ylim=c(0,7), yaxt="n")
axis(2, at=0:7)
interaction.plot(crash$passengers, crash$agerange, crash$crashrate,
                 ylim=c(0,7), yaxt="n")
axis(2, at=0:7)
mtext("Interactions for Crash Rates by Driver Age and Passenger Presence",
      outer=TRUE, side=3, line=-3, cex=1.2)
par(mfrow=c(1,1))
## export.eps(hh("dsgntwo/figure/crash-interaction.eps"))


## For black and white interaction2wt plots
## we prefer line types c(1,2,3,5,4,6,7) rather than the default 1:7
##   tpg.sl <- trellis.par.get("superpose.line")
##   tpg.sl.previous <- tpg.sl
##   tpg.sl$lty <- c(1,2,3,5,4,6,7)
##   trellis.par.set("superpose.line", tpg.sl)

## crash-original.eps.gz
interaction2wt(crashrate ~ agerange + passengers,
               data=crash,
               main="                                  a. original scale,  k = 1",
               ## above is necessary, below controls formatting.
               strip=FALSE,
               xlab="",
               main.cex=1.6,
               scales=list(
                 cex=1.2,
                 x=list(cex=.9),
                 y=list(cex=1, alternating=1)),
               responselab="",
               key.in=list(cex.title=1.4, cex=1, plot=FALSE))
## export.eps(hh("dsgntwo/figure/crash-original.eps"))

## original author's presentation
## crash-bar.eps.gz
old.par <- par(mfrow=c(2,3), oma=c(2,0,0,0))
for (i in levels(crash$agerange)) {
  barplot(crash$crashrate[crash$agerange==i],
          names=levels(crash$passengers),
          ylim=c(0,7), yaxt="n",
          xlab=paste("Ages", i),
          col=55)
  axis(2, at=0:7)
}
mtext("Crash Rates by Driver Age and Passenger Presence per 10,000 Trips",
      outer=TRUE, side=3, line=-3, cex=1.2)
mtext("Number of Passengers", outer=TRUE, side=1, line=-27, cex=1)
par(old.par)
## export.eps(hh("dsgntwo/figure/crash-bar.eps"))



## Based on odoffna paper, hh/articles/odoffna.pdf
## with detail based on UREDA

crash.aov <- aov(crashrate ~ agerange + passengers, data=crash, qr=TRUE)
summary(crash.aov)

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


## as.vector() is needed to respond to the change in definition of
## model.tables in S-Plus 6.1 from earlier releases.
## The class="table" was not precisely defined in earlier releases.
##
## cv: comparison value
cv <- outer(as.vector(mte$tables$agerange),
            as.vector(mte$tables$passengers)) / mtm$tables$"Grand mean"[1]
dimnames(cv) <- list(levels(crash$agerange), levels(crash$passengers))
cv

## diagnostic plot. UREDA page 200--204
## crash-diag.eps.gz
crashr.lm <- lm(resid(crash.aov) ~ as.vector(cv))
anova(crashr.lm)
coef(crashr.lm)
plot(resid(crash.aov) ~ as.vector(cv), pch=16,
     main=paste("resid(crash.aov) ~ as.vector(cv)\nslope =",
       round(coef(crashr.lm)[2],4)))
abline(crashr.lm)
## export.eps(hh("dsgntwo/figure/crash-diag.eps"))

crash2.aov <- aov(crashrate ~ agerange + passengers + as.vector(cv),
                  data=crash, qr=TRUE)
summary(crash2.aov)
coef(crash2.aov)


## The regression coefficient of the cv term is 1.551647, from either
## crashr.lm or crash2.lm.  The recommended power for a transformation
## is 1-1.551647 = -0.551647
## We illustrate power = 0 -.5 -1


interaction2wt(log(crashrate) ~ agerange + passengers,
               data=crash,
               main="                                 b. log scale,  k = 0",
               ## above is necessary, below controls formatting.
               strip=FALSE,
               xlab="",
               main.cex=1.6,
               scales=list(
                 cex=1.2,
                 x=list(cex=.9),
                 y=list(cex=1, alternating=1)),
               responselab="",
               key.in=list(cex.title=1.4, cex=1))
## export.eps(hh("dsgntwo/figure/crash-log.eps"))

interaction2wt( I(-1/sqrt(crashrate)) ~ agerange + passengers,
               data=crash,
               main="                                  c. negative reciprocal square root scale,  k = -.5",
               ## above is necessary, below controls formatting.
               strip=FALSE,
               xlab="",
               main.cex=1.6,
               scales=list(
                 cex=1.2,
                 x=list(cex=.9),
                 y=list(cex=1, alternating=1)),
               responselab="",
               key.in=list(cex.title=1.4, cex=1, plot=FALSE))
## export.eps(hh("dsgntwo/figure/crash-neg-rec-sqrt.eps"))

interaction2wt( I(-1/(crashrate)) ~ agerange + passengers,
               data=crash,
               main="                                  d. negative reciprocal scale,  k = -1",
               ## above is necessary, below controls formatting.
               strip=FALSE,
               xlab="",
               main.cex=1.6,
               scales=list(
                 cex=1.2,
                 x=list(cex=.9),
                 y=list(cex=1, alternating=1)),
               responselab="",
               key.in=list(cex.title=1.4, cex=1))
## export.eps(hh("dsgntwo/figure/crash-neg-rec.eps"))

## restore the original line types
##   trellis.par.set("superpose.line", tpg.sl.previous)


## The appearance of the -.5 and -1 transformations are similar.  We
## choose the reciprocal (power = -1) because it is easy to explain.
## The units are "crashes per mile".  I dropped the negative in the
## anova table and table of means.  I left the negative in the graph
## so it would go in the same order as the other graphs.

crashi.aov <- aov(1/crashrate ~ agerange + passengers, data=crash, qr=TRUE)
summary(crashi.aov)
model.tables(crashi.aov, type="means")


## Now that we can see the different behavior for the passengers
## conditional on the agerange, let us make the anova table show it.

## exploration of dummy variables
crashin.aov <- aov(1/crashrate ~ agerange/passengers, data=crash, qr=TRUE)
summary(crashin.aov)
coef(summary.lm(crashin.aov))
if.R(s=
summary(crashin.aov,
        split=list("passengers %in% agerange"=
          list(teens=1:2, adults=3, rest=4:9)))
,r=
summary(crashin.aov,
        split=list("agerange:passengers"=
          list(teens=1:2, adults=3, rest=4:9)))
)
model.tables(crashin.aov, type="means")

## selection of just the linear contrasts
pass <- as.numeric(crash$passengers)
crashinlin.aov <- aov(1/crashrate ~ agerange/pass,
                      data=crash, qr=TRUE)
summary(crashinlin.aov)
coef(summary.lm(crashinlin.aov))
coef(summary.lm(crashinlin.aov))[4:6,]

if.R(s=
summary(crashinlin.aov,
        split=list("pass %in% agerange"=
          list(teens=1:2, adults=3)))
,r=
summary(crashinlin.aov,
        split=list("agerange:pass"=
          list(teens=1:2, adults=3)))
)


## barplot display (original author's format) of the positive reciprocal
old.par <- par(mfrow=c(2,3), oma=c(2,0,0,0))
for (i in levels(crash$agerange)) {
  barplot(1/crash$crashrate[crash$agerange==i],
          names=levels(crash$passengers),
          ylim=c(0,3), yaxt="n",
          xlab=paste("Ages", i),
          col=55)
  axis(2, at=0:7)
}
mtext("Trips per .0001 Crashes by Driver Age and Passenger Presence",
      outer=TRUE, side=3, line=-3, cex=1.2)
mtext("Number of Passengers", outer=TRUE, side=1, line=-27, cex=1)
par(old.par)
## export.eps(hh("dsgntwo/figure/crash-bar-rec.eps"))



## This is the plot of the crash data following Tukey's method in the
## 1949 odoffna paper.  It has been superceded by the presentation above
## based on UREDA (Hoaglin, Mosteller, and Tukey 1983)
## odoffna plot, based on odoffna page 237
scp <- tapply(crash$crashrate, crash[,1:2], c) %*%
       as.vector(mte$tables$passengers)
plot(scp ~ as.vector(mtm$tables$agerange))
crash.lm <- lm(scp ~ as.vector(mtm$tables$agerange))
abline(crash.lm)

## error.limits formula, odoffna page 240
s2 <- sum(crash2.aov$residuals^2) / crash2.aov$df.residual
abline(a=coef(crash.lm)[1] - sum(mte$tables$passengers^2)^.5 * s2^.5,
       b=coef(crash.lm)[2])
abline(a=coef(crash.lm)[1] + sum(mte$tables$passengers^2)^.5 * s2^.5,
       b=coef(crash.lm)[2])



## crash-cover.s
## after crash.s

position(crash$agerange) <- c(1.2, 2.5, 3.7)
cover <-
interaction2wt( -1/(crashrate) ~ agerange + passengers,
               data=crash,
               main="                                  d. negative reciprocal scale,  k = -1",
               ## above is necessary, below controls formatting.
               strip=FALSE,
               xlab="",
               main.cex=1.6,
               scales=list(
                 cex=1.2,
                 x=list(cex=.9),
                 y=list(cex=1, alternating=1)),
               responselab="",
               key.in=list(cex.title=1.4, cex=1))
cover


s.l <- trellis.par.get("superpose.line")
s.l.old <- s.l
b.r <- trellis.par.get("box.rectangle")
b.r.old <- b.r
b.u <- trellis.par.get("box.umbrella")
b.u.old <- b.u

s.l$lwd[] <- 6
s.l$lty[] <- 1
trellis.par.set("superpose.line", s.l)
b.r$lwd[] <- 6
b.r$lty[] <- 1
trellis.par.set("box.rectangle", b.r)
b.u$lwd[] <- 6
b.u$lty[] <- 1
trellis.par.set("box.umbrella", b.u)
cover
## export.eps(hh("dsgntwo/figure/crash-neg-rec-color.eps"))


cover <-
interaction2wt( -1/(crashrate) ~ agerange + passengers,
               data=crash,
               main="",
               ## above is necessary, below controls formatting.
               strip=FALSE,
               xlab.in=FALSE,
               main.cex=1.6,
               scales=list(
                 cex=1.2,
                 x=list(cex=.9, tick=FALSE, labels=FALSE),
                 y=list(cex=1, alternating=1, labels=FALSE, tick=FALSE)),
               responselab="",
               key.in=list(cex.title=1.4, cex=1, plot=FALSE))
cover
if.R(r={
  cover$legend <- NULL
  cover <- update(cover, scales=list(alternating=0, tck=0))
  cover$lattice.options$axis.options <- list(bottom=NULL, right=NULL)
}, s={
  cover$factor.levels$agerange <- c("","","")
  cover$factor.levels$passengers <- c("","","","")
})

print(cover, position = c(-.2,0,1,.8))
## export.eps(hh("dsgntwo/figure/crash-neg-rec-cover.eps"))


## restore original settings
trellis.par.set("superpose.line", s.l.old)
trellis.par.set("box.rectangle", b.r.old)
trellis.par.set("box.umbrella", b.u.old)
