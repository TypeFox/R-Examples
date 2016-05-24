### R code from vignette source 'LSmeans.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: LSmeans.Rnw:54-57
###################################################
require( doBy )
prettyVersion <- packageDescription("doBy")$Version
prettyDate <- format(Sys.Date())


###################################################
### code chunk number 2: LSmeans.Rnw:95-96
###################################################
options("width"=90, "digits"=3)


###################################################
### code chunk number 3: LSmeans.Rnw:107-111
###################################################
dir.create("figures")
oopt <- options()
options("digits"=4, "width"=80, "prompt"="> ", "continue"="  ")
##options(useFancyQuotes="UTF-8")


###################################################
### code chunk number 4: LSmeans.Rnw:151-152 (eval = FALSE)
###################################################
## lm( y ~ treat + block + year)


###################################################
### code chunk number 5: LSmeans.Rnw:171-173 (eval = FALSE)
###################################################
## library(lme4)
## lmer( y ~ treat + (1|block) + (1|year))


###################################################
### code chunk number 6: LSmeans.Rnw:208-213
###################################################
simdat<-structure(list(treat = structure(c(1L, 1L, 1L, 2L, 1L, 2L, 2L, 2L
), .Label = c("t1", "t2"), class = "factor"), year = structure(c(1L,
1L, 1L, 1L, 2L, 2L, 2L, 2L), .Label = c("1", "2"), class = "factor"),
    y = c(0.5, 1, 1.5, 3, 3, 4.5, 5, 5.5)), .Names = c("treat", "year",
"y"), row.names = c(NA, -8L), class = "data.frame")


###################################################
### code chunk number 7: LSmeans.Rnw:218-219
###################################################
simdat


###################################################
### code chunk number 8: simdat-fig
###################################################
library(ggplot2)
qplot(treat, y, data=simdat, color=year, size=I(3))


###################################################
### code chunk number 9: LSmeans.Rnw:231-233
###################################################
msim <- lm(y ~ treat + year, data=simdat)
LSmeans( msim, effect="treat")


###################################################
### code chunk number 10: LSmeans.Rnw:237-238
###################################################
summaryBy(y~treat, data=simdat)


###################################################
### code chunk number 11: LSmeans.Rnw:250-253
###################################################
summary( warpbreaks )
head( warpbreaks, 4 )
ftable(xtabs( ~ wool + tension, data=warpbreaks))


###################################################
### code chunk number 12: LSmeans.Rnw:257-264
###################################################
 opar <- par(mfrow = c(1, 2), oma = c(0, 0, 1.1, 0))
     plot(breaks ~ tension, data = warpbreaks, col = "lightgray",
          varwidth = TRUE, subset = wool == "A", main = "Wool A")
     plot(breaks ~ tension, data = warpbreaks, col = "lightgray",
          varwidth = TRUE, subset = wool == "B", main = "Wool B")
     mtext("warpbreaks data", side = 3, outer = TRUE)
     par(opar)


###################################################
### code chunk number 13: LSmeans.Rnw:268-269
###################################################
(warp.lm <- lm(breaks ~ wool + tension, data=warpbreaks))


###################################################
### code chunk number 14: LSmeans.Rnw:274-276
###################################################
uni <- unique(warpbreaks[,2:3])
prd <- cbind(breaks=predict(warp.lm, newdata=uni), uni); prd


###################################################
### code chunk number 15: LSmeans.Rnw:291-292
###################################################
LSmeans(warp.lm, effect="tension")


###################################################
### code chunk number 16: LSmeans.Rnw:298-299
###################################################
doBy::summaryBy(breaks ~ tension, data=warpbreaks)


###################################################
### code chunk number 17: LSmeans.Rnw:311-313
###################################################
warp.lm2 <- update(warp.lm, .~.+wool:tension)
coef( summary( warp.lm2 ))


###################################################
### code chunk number 18: LSmeans.Rnw:318-320
###################################################
K2 <- LSmatrix(warp.lm2, effect="tension"); K2
linest(warp.lm2, K=K2)


###################################################
### code chunk number 19: chick-fig
###################################################
library(ggplot2)
ChickWeight$Diet <- factor(ChickWeight$Diet)
qplot(Time, weight, data=ChickWeight, colour=Chick, facets=~Diet,
      geom=c("point","line"))


###################################################
### code chunk number 20: LSmeans.Rnw:350-353
###################################################
library(lme4)
rr <- lmer(weight~Time*Diet + (0+Time|Chick), data=ChickWeight)
coef(summary(rr))


###################################################
### code chunk number 21: LSmeans.Rnw:359-360
###################################################
LSmatrix(rr, effect="Diet")


###################################################
### code chunk number 22: LSmeans.Rnw:367-368
###################################################
K1 <- LSmatrix(rr, effect="Diet", at=list(Time=1)); K1


###################################################
### code chunk number 23: LSmeans.Rnw:374-377
###################################################
K0 <- LSmatrix(rr, effect="Diet", at=list(Time=0))
K1-K0
LSmeans(rr, K=K1-K0)


###################################################
### code chunk number 24: LSmeans.Rnw:382-389
###################################################
LSmeans_trend <- function(object, effect, trend){

    K<-LSmatrix(object, effect=effect, at=as.list(setNames(1, trend))) -
        LSmatrix(object, effect=effect, at=as.list(setNames(0, trend)))
    LSmeans(object, K=K)
}
LSmeans_trend(rr, effect="Diet", trend="Time")


###################################################
### code chunk number 25: LSmeans.Rnw:397-402
###################################################
data(CO2)
CO2 <- transform(CO2, Treat=Treatment, Treatment=NULL)
levels(CO2$Treat) <- c("nchil","chil")
levels(CO2$Type) <- c("Que","Mis")
ftable(xtabs( ~ Plant + Type + Treat, data=CO2), col.vars=2:3)


###################################################
### code chunk number 26: co2-fig
###################################################
qplot(x=log(conc), y=uptake, data=CO2, color=Treat, facets=~Type, size=I(3))


###################################################
### code chunk number 27: LSmeans.Rnw:413-415
###################################################
co2.lm1 <- lm(uptake ~ conc + Type + Treat, data=CO2)
LSmeans(co2.lm1, effect="Treat")


###################################################
### code chunk number 28: LSmeans.Rnw:421-423 (eval = FALSE)
###################################################
## co2.lm <- lm(uptake ~ log(conc) + Type + Treat, data=CO2)
## LSmeans(co2.lm, effect="Treat")


###################################################
### code chunk number 29: LSmeans.Rnw:428-431
###################################################
co2.lm2 <- lm(uptake ~ log.conc + Type + Treat,
             data=transform(CO2, log.conc=log(conc)))
LSmeans(co2.lm2, effect="Treat")


###################################################
### code chunk number 30: LSmeans.Rnw:438-440
###################################################
co2.lm3 <- lm(uptake ~ conc + I(conc^2) + Type + Treat, data=CO2)
LSmeans(co2.lm3, effect="Treat")


###################################################
### code chunk number 31: LSmeans.Rnw:448-451
###################################################
co2.lm4 <- lm(uptake ~ conc + conc2 + Type + Treat, data=
              transform(CO2, conc2=conc^2))
LSmeans(co2.lm4, effect="Treat")


###################################################
### code chunk number 32: LSmeans.Rnw:456-457
###################################################
LSmeans(co2.lm4, effect="Treat", at=list(conc=10, conc2=100))


###################################################
### code chunk number 33: LSmeans.Rnw:476-479
###################################################
warp.poi <- glm(breaks ~ wool + tension, family=poisson, data=warpbreaks)
LSmeans(warp.poi, effect="tension", type="link")
LSmeans(warp.poi, effect="tension", type="response")


###################################################
### code chunk number 34: LSmeans.Rnw:490-493
###################################################
warp.qpoi <- glm(breaks ~ wool + tension, family=quasipoisson, data=warpbreaks)
LSmeans(warp.qpoi, effect="tension", type="link")
LSmeans(warp.qpoi, effect="tension", type="response")


###################################################
### code chunk number 35: LSmeans.Rnw:499-502
###################################################
warp.poi2 <- glm(breaks ~ wool + tension, family=poisson(link=identity),
                 data=warpbreaks)
LSmeans(warp.poi2, effect="tension", type="link")


###################################################
### code chunk number 36: LSmeans.Rnw:516-519
###################################################
warp.gam <- glm(breaks ~ wool + tension, family=Gamma(link=identity),
                 data=warpbreaks)
LSmeans(warp.gam, effect="tension", type="link")


###################################################
### code chunk number 37: LSmeans.Rnw:539-542
###################################################
warp.poi3 <- glm(breaks ~ wool + tension, family=quasipoisson(link=identity),
                 data=warpbreaks)
LSmeans(warp.poi3, effect="tension")


###################################################
### code chunk number 38: LSmeans.Rnw:550-553
###################################################
library(lme4)
warp.mm <- lmer(breaks ~ tension + (1|wool), data=warpbreaks)
LSmeans(warp.mm, effect="tension")


###################################################
### code chunk number 39: LSmeans.Rnw:561-562
###################################################
VarCorr(warp.mm)


###################################################
### code chunk number 40: LSmeans.Rnw:569-570
###################################################
LSmeans(warp.mm, effect="tension", adjust.df=FALSE)


###################################################
### code chunk number 41: LSmeans.Rnw:578-582
###################################################
library(geepack)
warp.gee <- geeglm(breaks ~ tension, id=wool, family=poisson, data=warpbreaks)
LSmeans(warp.gee, effect="tension")
LSmeans(warp.gee, effect="tension", type="response")


###################################################
### code chunk number 42: LSmeans.Rnw:593-594
###################################################
K <- LSmatrix(warp.lm, effect="tension"); K


###################################################
### code chunk number 43: LSmeans.Rnw:598-599
###################################################
linest( warp.lm, K=K )


###################################################
### code chunk number 44: LSmeans.Rnw:606-614
###################################################
## Make balanced dataset
dat.bal <- expand.grid(list(AA=factor(1:2), BB=factor(1:3), CC=factor(1:3)))
dat.bal$y <- rnorm(nrow(dat.bal))

## Make unbalanced dataset:  'BB' is nested within 'CC' so BB=1
## is only found when CC=1 and BB=2,3 are found in each CC=2,3,4
dat.nst <- dat.bal
dat.nst$CC <-factor(c(1,1,2,2,2,2,1,1,3,3,3,3,1,1,4,4,4,4))


###################################################
### code chunk number 45: LSmeans.Rnw:620-622
###################################################
head(dat.nst)
ftable(xtabs( ~ AA + BB + CC, data=dat.nst))


###################################################
### code chunk number 46: LSmeans.Rnw:627-629
###################################################
mod.nst  <- lm(y ~ AA + BB : CC, data=dat.nst)
coef( mod.nst )


###################################################
### code chunk number 47: LSmeans.Rnw:636-637
###################################################
LSmeans(mod.nst, effect=c("BB", "CC"))


###################################################
### code chunk number 48: LSmeans.Rnw:649-650
###################################################
X <- model.matrix( mod.nst ); as(X,"Matrix")


###################################################
### code chunk number 49: LSmeans.Rnw:664-666
###################################################
K <- LSmatrix(mod.nst, effect="BB", at=list(CC=2));K
LSmeans(mod.nst, K=K)


###################################################
### code chunk number 50: LSmeans.Rnw:677-680
###################################################
XtXinv <- MASS::ginv(t(X)%*%X)
bhat <- as.numeric(XtXinv %*% t(X) %*% dat.nst$y)
zapsmall(bhat)


###################################################
### code chunk number 51: LSmeans.Rnw:686-687
###################################################
K %*% bhat


###################################################
### code chunk number 52: LSmeans.Rnw:714-716
###################################################
S<-svd(X)
names(S)


###################################################
### code chunk number 53: LSmeans.Rnw:720-722
###################################################
B<-S$v[, S$d<1e-10, drop=FALSE ]; zapsmall(B) ## Basis for N(X)
zapsmall( rowSums(K%*%B) )


###################################################
### code chunk number 54: LSmeans.Rnw:736-739
###################################################
library("multcomp")
g1 <- glht(warp.lm, mcp(tension="Tukey"))
summary( g1 )


###################################################
### code chunk number 55: LSmeans.Rnw:744-745
###################################################
K1 <- g1$linfct; K1


