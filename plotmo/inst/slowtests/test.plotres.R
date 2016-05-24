# test.plotres.R

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

plotlm1 <- function(object)
{
    old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
    par(mfrow=c(2,2), oma=c(0,0,3,0), mar=c(4, 3, 3, 1.5),
        mgp=c(1.5, 0.4, 0), tcl=-.3, font.main=1, cex.main=1)
    plot(object, sub.caption="standard call to plot.lm")
}
plotlm.using.plotres <- function(object)
{
    old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
    par(mfrow=c(2,2), oma=c(0,0,3,0), mar=c(4, 3, 3, 1.5),
        mgp=c(1.5, 0.4, 0), tcl=-.3, font.main=1, cex.main=1)
                                                            # residuals vs fitted
    plotres(object, pch=1, which=3,
            caption=paste(deparse(object$call), collapse=" "))
                                                            # QQ plot
    plotres(object, pch=1, which=4, standardize=TRUE)
                                                            # scale-location plot
    plotres(object, pch=1, which=6, standardize=TRUE)
                                                            # leverage plot
    plotres(object, pch=1, which=3, versus=4, standardize=TRUE)
}
lm.mod <- lm(Volume~., data=trees)
plotlm1(lm.mod)
plotlm.using.plotres(lm.mod)

# various arguments

plotres(lm.mod, SHOWCALL=TRUE)
plotres(lm.mod, level=.95, id.n=-3, SHOWCALL=TRUE)
lm.tit <- lm(survived~., data=tit)
col <- ifelse(tit$survived, "green", "red")
pch <- ifelse(tit$sex == "male", 20, 6)
plotres(lm.tit, level=.95, col=col, pch=pch,
         level.shade="gray", level.shade2="lightgray", SHOWCALL=TRUE)
plotres(lm.tit, col.resp=3, cum.col=2, cum.cex=1.2, grid.col=5, qq.col=1, qq.cex=.3, SHOWCALL=TRUE)
plotres(lm.tit, pt.col="pink", smooth.col=0, SHOWCALL=TRUE)
plotres(lm.tit, smooth.col=3, smooth.lwd=1.2, smooth.lty=2, smooth.f=.2,
         label.col=4, label.cex=.9, label.font=2, SHOWCALL=TRUE)
foo <- function()
{
    afoo <- earth(O3~., data=ozone1, deg=2)
    old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
    par(mfrow=c(2,2), oma=c(0,0,3,0), mar=c(4, 3, 3, 1.5),
        mgp=c(1.5, 0.4, 0), tcl=-.3, font.main=1, cex.main=1)
    # test xlim ylim etc. on qq and cum plots
    plotres(afoo, which=2, trace=0, xlim=c(0,20), ylim=c(-.2,1.1), grid.col="pink", info=TRUE)
    plotres(afoo, which=2, trace=0,
            grid.col="pink", info=TRUE, cum.col=2, cum.cex=1.4)
    plotres(afoo, which=4)
    plotres(afoo, which=4, trace=0, xlim=c(-7,7), ylim=c(-20, 20),
            qq.col=2, qq.cex=.5, label.col=1, qqline.col="orange")
    # check xlim and ylim apply only to resids plots if multiple plots
    plotres(afoo, which=c(2:5), trace=0, xlim=c(-1,5), ylim=c(-8, 8),
            qq.col=2, qq.cex=.5, label.col=1, qqline.col="orange", smooth.col=3, smooth.lwd=2)
}
foo()

# test id.n and npoints
set.seed(1066)
a20 <- earth(Volume~., data=trees, ncr=3, nfo=3, varmod.method="lm", keepxy=TRUE)
par(mfrow=c(2,2), oma=c(0,0,3,0), mar=c(4, 3, 3, 1.5), mgp=c(1.5, 0.4, 0), tcl=-.3, cex=1)
plot(a20, which=3, standardize=TRUE, smooth.col=0, id.n=-1, main="a20-00, smooth.col=0, id.n=-1",
     caption="test id.n and npoints")
plot(a20, which=3, standardize=TRUE, smooth.col=0, id.n=10, main="a20-01, smooth.col=0, id.n=10")
# this tests cex with do.par=FALSE
plot(a20, which=3, standardize=TRUE, smooth.col=0, npoints=10, cex=.8, main="a20-02, smooth.col=0, npoints=10, cex=.8")
# TODO labels are hosed in the following
plot(a20, which=3, standardize=TRUE, smooth.col=0, npoints=5, id.n=10, main="a20-03, labels hosed\nsmooth.col=0, npoints=10, id.n=10")

# test leverages and handling of unity leverages
lm.mod <- lm(Volume~., data=trees)
par(mfrow=c(2,2), oma=c(0,0,3,0), mar=c(4, 3, 3, 1.5), mgp=c(1.5, 0.4, 0), tcl=-.3, cex=1)
a20$leverages[31] <- 1 # fake a unity leverage
plot(a20, which=3, versus=4, standardize=TRUE, main="resids vs leverage\nunity leverage",
     caption="leverage plots")
plotres(a20, which=3, standardize=TRUE, main="resids vs fitted\nunity leverage")
plotres(lm.mod, which=3, versus=4, standardize=TRUE, main="lever plot for lm.mod")
plotres(lm.mod, which=3, versus=4, standardize=TRUE, main="cook args",
        cook.levels=c(.5, .8, 1), cook.col="blue", cook.lty=2)

plot(a20, which=3, versus=4, standardize=TRUE, info=TRUE, main="resids vs leverage\nunity leverage",
     caption="leverage plots with info=TRUE")
plotres(a20, which=3, standardize=TRUE, info=TRUE, main="resids vs fitted\nunity leverage")
plotres(lm.mod, which=3, versus=4, standardize=TRUE, info=TRUE, main="lever plot for lm.mod")
plotres(lm.mod, which=3, versus=4, standardize=TRUE, info=TRUE, main="cook args",
        cook.levels=c(.5, .8, 1), cook.col="blue", cook.lty=2)

# back compat tests
par(mfrow=c(2,2), oma=c(0,0,3,0), mar=c(4, 3, 3, 1.5), mgp=c(1.5, 0.4, 0), tcl=-.3)
plotres(a20, which=3, col.smooth=4, smooth.lwd=2, smooth.lty=2,
        main="a20-04 col.smooth=4, smooth.lwd=2, smooth.lty=2",
        caption="back compat tests with plot.earth")
plotres(a20, which=4, qq.col=3,
             qqline.col="lightblue", qqline.lty=2, main="a20-05 qq.col=3")
plotres(a20, which=4, qqline.col=0, main="a20-06 qqline.col=0")
# set.seed(1066)
# mod.earth.tit <- earth(tit[,-3], tit[,3], degree=2, nfold=3, ncross=3, varmod.method="earth", keepxy=TRUE)
plot(0,0)
plot(a20, which=1, col.grid="pink", col.rsq=3, lty.rsq=1, main="a20-07 col.grid=\"pink\", col.rsq=3, lty.rsq=1")

# TODO following not working?
plot(a20, which=3, col.cv=4, col.grid="pink", main="a20-08 col.cv=4, col.grid=\"pink\"")

plot(a20, which=3, col.points="orange", cex.points=1.5, main="a20-09 col.points=\"orange\", cex.points=1.5")
plot(a20, which=3, col.residuals="orange", smooth.f=.2, col.line=3, main="a20-10 col.residuals=\"orange\", smooth.f=.2, col.line=3")

# test graphics args outside do.par
par(col.main="#456789")
cat("before par: cex=", par("cex"), " col.main=", par("col.main"), " col.axis=", par("col.axis"), "\n", sep="")
plot(a20, which=c(2,3),          caption="a20 which=c(2,3) (i.e. do.par=TRUE)  no cex")
plot(a20, which=c(2,3), cex=1,   caption="a20 which=c(2,3) (i.e. do.par=TRUE)  cex=1, plot should be identical to previous page")
plot(a20, which=c(2,3), cex=1.2, caption="a20 which=c(2,3) (i.e. do.par=TRUE)  cex=1.2")
plot(a20, which=3,         main="no cex", caption="a20 test graphics args with do.par=FALSE")
plot(a20, which=3, cex=1,  main="cex=1")
plot(a20, which=3, cex=.8, main="cex=.8")
plot(a20, which=3, cex=1.1, col.main=2, col.axis="blue", col.lab=3, font.lab=2,
          main="cex=1.1, col.main=2, col.axis=\"blue\", col.lab=3, font.lab=2")
# all of these should have been restored
cat("after par: cex=", par("cex"), " col.main=", par("col.main"), " col.axis=", par("col.axis"), "\n", sep="")
stopifnot(par("col.main") == "#456789")
par(col.main=1)

survived <- as.numeric(tit$survived) # 0 or 1
sex      <- as.numeric(tit$sex)      # 1 or 2
pclass   <- as.numeric(tit$pclass)   # 1,2, or 3
age      <- tit$age                  # .2 to 80

printf("======== basic operation, compare to plot.lm etc.\n")
par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(1.5,0.5,0), oma=c(0,0,2.5,0))
lm <- lm(survived~sex+pclass+age)
plot(lm, which=5, pch=20)
plot(0, 0)
plot(lm, which=1, pch=20)
plot(lm, which=2, pch=20)
plotres(lm, standardize=1, cook.levels=c(.1,.2,.3), SHOWCALL=TRUE)
elm <- earth(survived~sex+pclass+age, linpreds=TRUE, thresh=0, penalty=-1)
plotres(elm, col=survived+2, SHOWCALL=TRUE)
plotres(elm, col=survived+2, SHOWCALL=TRUE)
plotres(elm, col=survived+2, col.rsq="darkorange", lty.rsq=1, SHOWCALL=TRUE)
set.seed(2015)
elm.glm <- earth(survived~sex+pclass+age, linpreds=TRUE, thresh=0, penalty=-1,
                 glm=list(family=binomial),
                 ncr=3, nfold=3, varmod.method="lm")
plotres(elm.glm, col=survived+2, SHOWCALL=TRUE)

printf("======== check type arg with earth\n")
par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(1.5,0.5,0), oma=c(0,0,2.5,0))
# following two are equivalent
plotres(elm.glm, col=survived+2, standardize=TRUE,
        which=3, do.par=FALSE, main="standardize=TRUE")
mtext("elm.glm with various type options", outer=TRUE, font=2, line=1, cex=1)
plotres(elm.glm, col=survived+2, type="standardize",
        which=3, do.par=FALSE, main="type=\"standardize\"\nequivalent to standardize=TRUE")
# TODO double standardization, should not be allowed
plotres(elm.glm, col=survived+2, standardize=TRUE, type="standardize",
        which=3, do.par=FALSE,
        main="standard=TRUE, type=\"deviance\"\ndouble standardization")
plotres(elm.glm, col=survived+2, type="deviance",
        which=3, do.par=FALSE, main="type=\"deviance\"")

printf("======== multiple response earth models\n")
par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(1.5,0.5,0), oma=c(0,0,2.5,0))
set.seed(2015)
emulti <- earth(cbind(Volume, Volume + 100 + 5 * rnorm(nrow(trees)))~., data=trees)
plot(emulti, nresponse=2,
     which=3, do.par=FALSE, main="emulti nresponse=2")
mtext("multiple response earth models", outer=TRUE, font=2, line=1, cex=1)
plot(emulti, nresponse=2, FORCEPREDICT=TRUE,
     which=3, do.par=FALSE, main="emulti, nresponse=2\nFORCEPREDICT=TRUE")

printf("======== earth model with a factor response\n")
epclass <- earth(pclass~., data=tit)
par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(1.5,0.5,0), oma=c(0,0,2.5,0))
set.seed(2015)
plot(epclass, nresponse="first", trace=1,
     which=3, do.par=FALSE, main="pclass response, nresponse=\"first\"")
mtext("earth model with a factor response", outer=TRUE, font=2, line=1, cex=1)
plot(epclass, nresponse="first", trace=1, FORCEPREDICT=TRUE,
     which=3, do.par=FALSE,
     main="pclass response, nresponse=\"first\"\nFORCEPREDICT=TRUE")

printf("======== glm\n")
glm <- glm(survived~sex+pclass+age, family=binomial)
par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(1.5,0.5,0), oma=c(0,0,2.5,0))
plot(glm, which=1, pch=20, main="plot.lm")
mtext("glm model with plot.lm and plotres", outer=TRUE, font=2, line=1, cex=1)
plotres(glm, which=3, main="plotres glm survived")
# with plotres we can also plot pearson etc. residuals
plotres(glm, which=3, type="pearson", main="plotres glm survived\ntype=\"pearson\"")

printf("======== rpart\n")
library(rpart)
par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(1.5,0.5,0), oma=c(0,0,2.5,0))
rpart <- rpart(survived~sex+pclass+age)
plotres(rpart, SHOWCALL=TRUE)
plotres(rpart, SHOWCALL=TRUE, FORCEPREDICT=TRUE) # identical
# TODO following fails in plotmo.predict.rpart (which is called to get the fitted values)
# plotres(rpart, type="pearson")
plotres(rpart, jitter=3, w1.extra=100, w1.under=TRUE, w1.branch.type=5,
        col=survived+2, smooth.col=NA, label.col=1, SHOWCALL=TRUE)

fit <- rpart(Kyphosis ~ Age + Number + Start, data = kyphosis)
plotres(fit, nresponse=1, SHOWCALL=TRUE, jitter=5)
plotres(fit, nresponse=2, SHOWCALL=TRUE, jitter=TRUE)

printf("======== versus=\"b:\"\n")

library(gam)
gam.package.loaded  <- "package:gam"  %in% search()
mgcv.package.loaded <- "package:mgcv" %in% search()
if(mgcv.package.loaded && gam.package.loaded) {
    # prevent downstream confusing error messages
    stop0("both 'gam' and 'mgcv' are loaded")
}
library(earth)
data(ozone1)
data(ozone1)
oz <- ozone1[, c("O3", "humidity", "temp", "ibt")]
gam.mod <- gam(O3^(1/3) ~ lo(humidity)+lo(ibt,temp), data=oz)
plotmo(gam.mod, SHOWCALL=TRUE)
plotres(gam.mod, SHOWCALL=TRUE)
plotres(gam.mod, versus="b:", SHOWCALL=TRUE)
plotres(gam.mod, versus="b:ib", info=TRUE, SHOWCALL=TRUE)

gam.linear.humidity.only <- gam(O3^(1/3) ~ humidity, data=oz)
plotres(gam.linear.humidity.only, versus="b:", SHOWCALL=TRUE)

library(mda)
mars <- mars(ozone1[,2:3], ozone1[,1], degree=2)
mars.to.earth <- mars.to.earth(mars)
plotres(mars,  versus="b:", caption="mars model, versus=\"b:\"", SHOWCALL=TRUE)
plotres(mars.to.earth, versus="b:", caption="earth model, versus=\"b:\", should be same as previous page", SHOWCALL=TRUE)
plotres(mars, versus="b:1", caption="mars model, versus=\"b:1\"", SHOWCALL=TRUE)

# lars is tested in plotmo3.R
# gbm  is tested in plotmo3.R
# TODO fda is not tested

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}
