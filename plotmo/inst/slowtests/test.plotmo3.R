# test.plotmo3.R: extra tests for plotmo version 3 and higher

options(warn=1) # print warnings as they occur

if(!interactive())
    postscript(paper="letter")

printf <- function(format, ...) cat(sprintf(format, ...), sep="") # like c printf

strip.space <- function(s) gsub("[ \t\n]", "", s)

# test that we got an error as expected from a try() call
expect.err <- function(object, expected.msg="")
{
    if(class(object)[1] == "try-error") {
        msg <- attr(object, "condition")$message[1]
        if(length(grep(expected.msg, msg, fixed=TRUE)))
            cat("Got error as expected from ",
                deparse(substitute(object)), "\n", sep="")
        else
            stop(sprintf("Expected: %s\n  Got:      %s",
                         expected.msg, substr(msg, 1, 1000)))
    } else
        stop("did not get expected error ", expected.msg)
}
printf("library(earth)\n")
library(earth)
data(ozone1)
data(etitanic)

check.naken <- function(s, expected)
{
    nude <- plotmo:::naken.formula.string(s)
    printf("%-60.60s %-s\n", s, nude)
    stopifnot(nude == expected)
}
printf("=== checking naken.formula.string\n")
check.naken("y ~ x1 : x2 + x3", "y~x1+x2+x3")
check.naken("y ~ x1 + x2 - x3", "y~x1+x2+x3")
check.naken("cbind(damage, 6-damage)~temp", "cbind(damage,6-damage)~temp")
check.naken("depIndex~q_4+q_2102+q_2104+q_3105+q_3106", "depIndex~q_4+q_2102+q_2104+q_3105+q_3106")
check.naken("doy ~ (vh+wind+humidity)^2", "doy~vh+wind+humidity")
check.naken("doy ~ s(wind) + s(humidity,wind) + s(vh)", "doy~wind+humidity+vh")
check.naken("log(doy) ~ I(vh*wind) + I(humidity*temp)+log(doy)", "log(doy)~vh+wind+humidity+temp+doy")
check.naken("log(doy)~vh+wind+humidity+I(wind*humidity)+temp+log(ibh)", "log(doy)~vh+wind+humidity+temp+ibh")
check.naken("O3 ~ s(humidity)+s(temp)+s(ibt)+s(temp,ibt)", "O3~humidity+temp+ibt")
check.naken("Ozone^(1/3) ~ lo(Solar.R) + lo(Wind, Temp)", "Ozone^(1/3)~Solar.R+Wind+Temp")
check.naken("Volume~(Girth*Height2)-Height", "Volume~Girth+Height2+Height")
check.naken("y ~ s(x) + s(x,z1)", "y~x+z1")
check.naken("y~s(x0,x1,k=12)+s(x2)+s(x3,k=20,fx=20)", "y~x0+x1+x2+x3")
check.naken("y~x[,1]+x[,2]", "y~x[,1]+x[,2]")
check.naken("y~x[,1]+x[,my.list$j]", "y~x[,1]+x[,my.list$j]")
check.naken("y~x[,i]+x[,2]", "y~x[,i]+x[,2]")
check.naken("Salary~Hitters[,1]", "Salary~Hitters[,1]")
check.naken("Salary~Hitters[,-1]", "Salary~Hitters[,-1]")
check.naken("Salary~Hitters[,c(1,2)]", "Salary~Hitters[,c(1,2)]")
check.naken("Salary~Hitters[,1:2]", "Salary~Hitters[,1:2]")
check.naken("Salary~Hitters[,c(1,2)]", "Salary~Hitters[,c(1,2)]")
check.naken("x[,c(1,2)] + x[,3]", "x[,c(1,2)]+x[,3]")
check.naken("x[,1] + x[,2] + x[,3] + x[,29] + x[,-14]", "x[,1]+x[,2]+x[,3]+x[,29]+x[,-14]")
check.naken("x[,c(1,2)] + x[,3] + x[,5:6] + x[,-1]", "x[,c(1,2)]+x[,3]+x[,5:6]+x[,-1]")
check.naken("log(y) ~ x9 + ns(x2,4) + s(x3,x4,df=4) + x5:sqrt(x6)", "log(y)~x9+x2+x3+x4+x5+x6")

plotmo1 <- function(object, ..., trace=0, SHOWCALL=TRUE, caption=NULL) {
    if(is.null(caption))
        caption <- paste(deparse(substitute(object)), collapse=" ")
    call <- match.call(expand.dots=TRUE)
    call <- strip.space(paste(deparse(substitute(call)), collapse=" "))
    printf("%s\n", call)
    plotmo(object, trace=trace, SHOWCALL=SHOWCALL, caption=caption, ...)
}
plotres1 <- function(object, ..., trace=0, SHOWCALL=TRUE, caption=NULL) {
    if(is.null(caption))
        caption <- paste(deparse(substitute(object)), collapse=" ")
    call <- match.call(expand.dots=TRUE)
    call <- strip.space(paste(deparse(substitute(call)), collapse=" "))
    printf("%s\n", call)
    plotres(object, trace=trace, SHOWCALL=SHOWCALL, caption=caption, ...)
}
# basic tests of plotmo on abbreviated titanic data

get.tit <- function()
{
    tit <- etitanic
    pclass <- as.character(tit$pclass)
    # change the order of the factors so not alphabetical
    pclass[pclass == "1st"] <- "first"
    pclass[pclass == "2nd"] <- "class2"
    pclass[pclass == "3rd"] <- "classthird"
    tit$pclass <- factor(pclass, levels=c("class2", "classthird", "first"))
    # log age is so we have a continuous predictor even when model is age~.
    set.seed(2015)
    tit$logage <- log(tit$age) + rnorm(nrow(tit))
    tit$parch <- NULL
    # by=12 gives us a small fast model with an additive and a interaction term
    tit <- tit[seq(1, nrow(etitanic), by=12), ]
}
tit <- get.tit()

mod.lm.age <- lm(age~., data=tit)
plotmo1(mod.lm.age)
plotmo1(mod.lm.age, level=.95)
plotmo1(mod.lm.age, level=.95, col.resp=3)

sexn <- as.numeric(tit$sex)
mod.lm.sexn <- lm(sexn~.-sex, data=tit)
plotmo1(mod.lm.sexn)
plotmo1(mod.lm.sexn, level=.95)

mod.earth.age <- earth(age~., data=tit, degree=2, nfold=3, ncross=3, varmod.method="lm")
plotmo1(mod.earth.age)
plotmo1(mod.earth.age, level=.9, degree2=0)

# tit[,4] is age
mod.earth.tit <- earth(tit[,-4], tit[,4], degree=2, nfold=3, ncross=3, varmod.method="lm")
plotmo1(mod.earth.tit)
plotmo1(mod.earth.tit, level=.9, degree2=0)

a.earth.sex <- earth(sex~., data=tit, degree=2, nfold=3, ncross=3, varmod.method="lm")
plotmo1(a.earth.sex)
plotmo1(a.earth.sex, level=.9)
plotmo1(a.earth.sex, type="class")
expect.err(try(plotmo1(a.earth.sex, level=.9, degree2=0, type="class")), "predicted values are strings")

# tit[,3] is sex
mod.earth.tit <- earth(tit[,-3], tit[,3], degree=2, nfold=3, ncross=3, varmod.method="lm")
plotmo1(mod.earth.tit)
plotmo1(mod.earth.tit, level=.9, degree2=0)
plotmo1(mod.earth.tit, type="class")
expect.err(try(plotmo1(mod.earth.tit, level=.9, degree2=0, type="class")), "predicted values are strings")

mod.earth.sex <- earth(sex~., data=tit, degree=2, nfold=3, ncross=3, varmod.method="earth", glm=list(family=binomial))
plotmo1(mod.earth.sex)
plotmo1(mod.earth.sex, type="link")
plotmo1(mod.earth.sex, type="class")
plotmo1(mod.earth.sex, level=.9, type="earth")

# tit[,3] is sex
mod.earth.tit <- earth(tit[,-3], tit[,3], degree=2, nfold=3, ncross=3, varmod.method="earth", glm=list(family=binomial))
plotmo1(mod.earth.tit)
plotmo1(mod.earth.tit, type="link")
plotmo1(mod.earth.tit, type="class")
plotmo1(mod.earth.tit, level=.9, type="earth")

# check factor handling when factors are not ordered alphabetically
tit.orgpclass <- etitanic[seq(1, nrow(etitanic), by=12), ]
tit  <- get.tit()
tit$logage <- NULL
tit.orgpclass$parch <- NULL
stopifnot(names(tit.orgpclass) == names(tit))
a.tit.orgpclass <- earth(pclass~., degree=2, data=tit.orgpclass)
a.tit           <- earth(pclass~., degree=2, data=tit)
expect.err(try(plotmo(a.tit)), "predict.earth returned multiple columns")
# following two graphs should be identical
plotmo1(a.tit.orgpclass, nresponse="1st",   all1=T, col.resp=3, type2="im")
plotmo1(a.tit,           nresponse="first", all1=T, col.resp=3, type2="im")
# following two graphs should be identical
plotmo1(a.tit.orgpclass, nresponse="2nd",    all1=T)
plotmo1(a.tit,           nresponse="class2", all1=T)

tit  <- get.tit()
mod.earth.pclass <- earth(pclass~., data=tit, degree=2)
expect.err(try(plotmo1(mod.earth.pclass)), "nresponse is not specified")
plotmo1(mod.earth.pclass, nresponse="fi")
plotmo1(mod.earth.pclass, nresponse="first")
plotmo1(mod.earth.pclass, nresponse=3)
plotmo1(mod.earth.pclass, type="class")
plotmo1(mod.earth.pclass, nresponse=1,
       type="class", grid.levels=list(sex="fem"),
       smooth.col="indianred", smooth.lwd=2,
       pt.col=as.numeric(tit$pclass)+1,
       pt.pch=1)

# tit[,1] is pclass
mod.earth.tit <- earth(tit[,-1], tit[,1], degree=2)
expect.err(try(plotmo1(mod.earth.tit)), "nresponse is not specified")
plotmo1(mod.earth.tit, nresponse="first")
plotmo1(mod.earth.tit, type="class")

mod.earth.pclass2 <- earth(pclass~., data=tit, degree=2, glm=list(family=binomial))
# expect.err(try(plotmo1(mod.earth.pclass2)), "nresponse is not specified")
plotmo1(mod.earth.pclass2, nresponse=3)
plotmo1(mod.earth.pclass2, type="link", nresponse=3)
plotmo1(mod.earth.pclass2, type="class")

# tit[,1] is pclass
mod.earth.tit <- earth(tit[,-1], tit[,1], degree=2, glm=list(family=binomial))
plotmo1(mod.earth.tit, nresponse=3)
plotmo1(mod.earth.tit, type="link", nresponse=3)
plotmo1(mod.earth.tit, type="class")

# plotmo vignette examples

# use a small set of variables for illustration
printf("library(earth)\n")
library(earth) # for ozone1 data
data(ozone1)
oz <- ozone1[, c("O3", "humidity", "temp", "ibt")]

lm.model.vignette <- lm(O3 ~ humidity + temp*ibt, data=oz) # linear model
plotmo1(lm.model.vignette, pt.col="gray", nrug=-1)
plotmo1(lm.model.vignette, level=.9)

printf("library(mda)\n")
library(mda)
mars.model.vignette1 <- mars(oz[,-1], oz[,1], degree=2)
plotmo1(mars.model.vignette1)
plotres1(mars.model.vignette1)
mars.model.vignette2 <- mars(oz[,-1,drop=FALSE], oz[,1,drop=FALSE], degree=2)
plotmo1(mars.model.vignette2)
# TODO causes Error in lm.fit(object$x, y, singular.ok = FALSE) : (list) object cannot be coerced to type 'double'
#      although still works
#      the error is mars.to.earth try(hatvalues.lm.fit(lm.fit(object$x, y, singular.ok=FALSE)))
plotres1(mars.model.vignette2, trace=1)

printf("library(rpart)\n")
library(rpart)                                          # rpart
rpart.model.vignette <- rpart(O3 ~ ., data=oz)
plotmo1(rpart.model.vignette, all2=TRUE)
expect.err(try(plotmo1(rpart.model.vignette, level=.9)), "not supported for rpart objects")

# commented out because is slow and already tested in test.non.earth.R
# printf("library(randomForest)\n")
# library(randomForest)                                   # randomForest
# rf.model.vignette <- randomForest(O3~., data=oz)
# plotmo1(rf.model.vignette)
# partialPlot(rf.model.vignette, oz, temp) # compare to partial-dependence plot

printf("library(gbm)\n")
library(gbm)                                            # gbm
gbm.model.vignette <- gbm(O3~., data=oz, dist="gaussian", inter=2, n.trees=100)
# commented out following because they always take the whole page
# plot(gbm.model.vignette, i.var=2) # compare to partial-dependence plots
# plot(gbm.model.vignette, i.var=c(2,3))
plotmo1(gbm.model.vignette, caption="gbm.model.vignette")

# commented out because is slow and already tested elsewhere
# printf("library(mgcv)\n")
# library(mgcv)                                           # gam
# gam.model.vignette <- gam(O3 ~ s(humidity)+s(temp)+s(ibt)+s(temp,ibt), data=oz)
# plotmo1(gam.model.vignette, level=.95, all2=TRUE)

printf("library(nnet)\n")
library(nnet)                                           # nnet
set.seed(4)
nnet.model.vignette <- nnet(O3~., data=scale(oz), size=2, decay=0.01, trace=FALSE)
plotmo1(nnet.model.vignette, type="raw", all2=T)

printf("library(MASS)\n")
library(MASS)                                           # qda
lcush <- data.frame(Type=as.numeric(Cushings$Type),log(Cushings[,1:2]))
lcush <- lcush[1:21,]
qda.model.vignette <- qda(Type~., data=lcush)
plotmo1(qda.model.vignette, type="class", all2=TRUE,
       type2="contour", ngrid2=100, nlevels=2, drawlabels=FALSE,
       pt.col=as.numeric(lcush$Type)+1,
       pt.pch=as.character(lcush$Type))

# miscellaneous other examples

tit <- get.tit()

mod.glm.sex <- glm(sex~., data=tit, family=binomial)
plotmo1(mod.glm.sex, pt.col=as.numeric(tit$pclass)+1)

# tit[,4] is age, tit[,1] is pclass
printf("library(lars)\n")
library(lars)
set.seed(2015)
xmat <- as.matrix(tit[,c(2,5,6)])
mod.lars.xmat <- lars(xmat, tit[,4])
par(mfrow=c(2,2))
plot(mod.lars.xmat)
plotmo1(mod.lars.xmat, nresponse=4, do.par=F)
plotres(mod.lars.xmat, trace=0, nresponse=4)

printf("library(cosso)\n")
library(cosso)
cosso <- cosso(xmat,tit[,4],family="Gaussian")
# TODO tell maintainer of cosso that you have to do this
class(cosso) <- "cosso"
plotmo1(cosso)
plotres(cosso)

# examples from James, Witten, et al. ISLR book
# I tested all models in their scripts manually.
# All worked except for exceptions below.

printf("library(pls)\n")
library(pls)
printf("library(ISLR)\n")
library(ISLR)
Hitters=na.omit(Hitters)

set.seed(1)
x <- model.matrix(Salary~.,Hitters)[,-1]
y <- Hitters$Salary
train=sample(1:nrow(x), nrow(x)/2)
pcr.fit1=pcr(Salary~., data=Hitters,subset=train,scale=TRUE, validation="CV")
plotmo1(pcr.fit1, nresponse=10)

# set.seed(1)
# x <- model.matrix(Salary~.,Hitters)[,-1]
# y <- Hitters$Salary
# train=sample(1:nrow(x), nrow(x)/2)
# pcr.fit2=pcr(y~x,scale=TRUE,ncomp=7)
# # TODO following gives Error: predictions returned the wrong length (got 263 but expected 50)
# plotmo1(pcr.fit2, nresponse=5)

library(splines)
fit.lm2=lm(wage~bs(age,knots=c(25,40,60)),data=Wage)
par(mfrow=c(1,2),mar=c(4.5,4.5,1,1),oma=c(0,0,4,0))
agelims=range(Wage$age)
age.grid=seq(from=agelims[1],to=agelims[2])
pred=predict(fit.lm2,newdata=list(age=age.grid),se=T)
plot(Wage$age,Wage$wage,col="gray", ylim=c(0,320))
lines(age.grid,pred$fit,lwd=2)
lines(age.grid,pred$fit+2*pred$se,lty="dashed")
lines(age.grid,pred$fit-2*pred$se,lty="dashed")
fit.lm2=lm(wage~bs(age,knots=c(25,40,60)),data=Wage,model=F) # TODO delete
plotmo1(fit.lm2, col.resp=2, do.par=F, level=.95, ylim=c(0,320),
        nrug=TRUE, caption="fit.lm2", ylab="wage")

fit.glm2 <- glm(I(wage>250)~poly(age,4),data=Wage,family=binomial)
par(mfrow=c(1,2),mar=c(4.5,4.5,1,1),oma=c(0,0,4,0))
agelims=range(Wage$age)
age.grid=seq(from=agelims[1],to=agelims[2])
# their plot
preds=predict(fit.glm2,newdata=list(age=age.grid),se=T)
pfit=exp(preds$fit)/(1+exp(preds$fit))
se.bands.logit = cbind(preds$fit+2*preds$se.fit, preds$fit-2*preds$se.fit)
se.bands = exp(se.bands.logit)/(1+exp(se.bands.logit))
preds=predict(fit.glm2,newdata=list(age=age.grid),type="response",se=T)
plot(Wage$age,I(Wage$wage>250),xlim=agelims,type="n",ylim=c(0,.2))
points(jitter(Wage$age), I((Wage$wage>250)/5),cex=.5,pch="|",col="darkgrey")
lines(age.grid,pfit,lwd=2, col="blue")
matlines(age.grid,se.bands,lwd=1,col="blue",lty=3)
# plotmo plot, side by side
# TODO Warning: the level argument may not be properly supported on glm objects built with weights
plotmo1(fit.glm2, level=.95, degree1.col="blue", ylim=c(0,.2), do.par=FALSE, nrug=-1, caption="fit.glm2", ylab="I(wage > 250)")

# Test deparsing of the formula in plotmo.pairs.default
# TODO Height is included in the plots even though formula says -Height
Height2 <- trees$Height^2
a <- lm(Volume~(Girth*Height2)-Height, data=trees, x=TRUE, model=FALSE)
plotmo(a)

# test "the variable on the right side of the formula is a matrix or data.frame"
# TODO would like to solve this problem

old.warn <- options("warn")
options(warn=2)

data(gasoline, package="pls")
earth.octane <- earth(octane ~ NIR, data=gasoline)
expect.err(try(plotmo(earth.octane)), "the variable on the right side of the formula is a matrix or data.frame")

library(ElemStatLearn)
x <- mixture.example$x
g <- mixture.example$y
lm.mixture.example <- lm(g ~ x)
expect.err(try(plotmo(lm.mixture.example)), "the variable on the right side of the formula is a matrix or data.frame")

options(warn=old.warn$warn)

# test variable names with $ are not supported

a <- earth(O3~ozone1$doy, data=ozone1)
expect.err(try(plotmo(a)), "cannot get the original model predictors")

a <- earth(O3~ozone1$doy + temp, data=ozone1)
expect.err(try(plotmo(a)), "cannot get the original model predictors")

a <- lm(O3~ozone1$doy, data=ozone1)
expect.err(try(plotmo(a)), "cannot get the original model predictors")

a <- lm(O3~ozone1$doy + temp, data=ozone1)
expect.err(try(plotmo(a)), "cannot get the original model predictors")

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}
