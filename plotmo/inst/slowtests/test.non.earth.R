# test.non.earth.R: test plotmo on non-earth models
# Stephen Milborrow, Basley KwaZulu-Natal Mar 2011

library(plotmo)
library(earth)
data(ozone1)
data(etitanic)
options(warn=1) # print warnings as they occur
if(!interactive())
    postscript(paper="letter")
dopar <- function(nrows, ncols, caption = "")
{
    cat("                             ", caption, "\n")
    par(mfrow=c(nrows, ncols))
    par(oma = c(0, 0, 3, 0))
    par(mar = c(3, 3, 1.7, 0.5))
    par(mgp = c(1.6, 0.6, 0))
    par(cex = 0.7)
}
expect.err <- function(obj) # test that we got an error as expected from a try() call
{
    if(class(obj)[1] == "try-error")
        cat("Got error as expected\n")
    else
        stop("did not get expected error")
}
caption <- "test lm(log(doy) ~ vh+wind+humidity+temp+log(ibh), data=ozone1)"
dopar(4,5,caption)
a <- lm(log(doy) ~ vh + wind + humidity + temp + log(ibh), data=ozone1)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA, col.response=3, pt.pch=20, smooth.col="indianred")
termplot(a)

caption <- "test lm(log(doy) ~ vh+wind+humidity+I(wind*humidity)+temp+log(ibh), data=ozone1)"
dopar(4,5,caption)
a <- lm(log(doy) ~ vh + wind + humidity + temp + log(ibh), data=ozone1)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA, col.resp=3, pt.pch=20, clip=FALSE, smooth.col="indianred")
termplot(a)

caption <- "test lm(doy ~ (vh+wind+humidity)^2, data=ozone1)"
dopar(4,3,caption)
a <- lm(doy ~ (vh+wind+humidity)^2, data=ozone1)
plotmo(a, do.par=FALSE, caption=caption, ylim=NULL)
# termplot(a) # termplot fails with Error in `[.data.frame`(mf, , i): undefined columns selected

caption <- "test lm(doy^2 ~ vh+wind+humidity+I(wind*humidity)+temp+log(ibh), data=ozone1)"
dopar(4,3,caption)
a <- lm(doy^2 ~ vh+wind+humidity+I(wind*humidity)+temp+log(ibh), data=ozone1)
plotmo(a, do.par=FALSE, caption=caption, ylim=NULL)
termplot(a) # termplot draws a funky second wind plot

caption <- "test lm with data=ozone versus attach(ozone)"
dopar(4,3,caption)
a <- lm(log(doy) ~ I(vh*wind) + wind + I(humidity*temp) + log(ibh), data=ozone1)
plotmo(a, do.par=FALSE, caption=caption, degree1=c(1,2,4,5))
attach(ozone1)
a <- lm(log(doy) ~ I(vh*wind) + wind + I(humidity*temp) + log(ibh))
plotmo(a, do.par=FALSE, degree1=c(1,2,4,5))
detach(ozone1)

# commented out because "$" in names is not yet supported
# a <- lm(log(ozone1$doy) ~ I(ozone1$vh*ozone1$wind) + log(ozone1$ibh))
# plotmo(a)

set.seed(1)
caption <- "test lm and glm a900..a902: damage~temp family=binomial data=orings"
dopar(2,3,caption)
library(faraway)
data(orings)
a900 <- lm(I(damage/6) ~ temp, data=orings)
plotmo(a900, do.par=FALSE, caption=caption, col.response=2, nrug=-1,
    main="lm(damage/6~temp)", smooth.col="indianred", trace=0)
response <- cbind(orings$damage, 6-orings$damage)
a901 <- glm(response ~ temp, family="binomial", data=orings)
plotmo(a901, do.par=FALSE, col.response=2, nrug=-1,
       main="glm(response~temp)", smooth.col="indianred", trace=2)
a902 <- glm(cbind(damage, 6-damage)~temp, family="binomial", data=orings)
plotmo(a902, do.par=FALSE, col.response=2, nrug=TRUE,
       main="glm(cbind(damage,6-damage)~temp)", trace=0)
termplot(a902, main="termplot")
plotmo(a902, type="link", main="type=\"link\"", do.par=F)
plotmo(a902, type="response", main="type=\"response\"", col.response=2, do.par=F)
par(mfrow=c(1,1))

set.seed(1)
caption <- "test glm(lot2~log(u),data=clotting,family=Gamma)"
dopar(2,2,caption)
u = c(5,10,15,20,30,40,60,80,100)
lota = c(118,58,42,35,27,25,21,19,18)
clotting <- data.frame(u = u, lota = lota)
a <- glm(lota ~ log(u), data=clotting, family=Gamma)
plotmo(a, do.par=FALSE, caption=caption, col.response=3, clip=FALSE, nrug=-1)
termplot(a)
plotmo(a, type="link", caption=paste("type=\"link\"", caption))

if(length(grep("package:gam", search())))
    detach("package:gam")
library(mgcv)
set.seed(1)
caption <- "test plot.gam, with mgcv::gam(y ~ s(x) + s(x,z)) with response and func (and extra image plot)"
dopar(3,2,caption)
par(mar = c(3, 5, 1.7, 0.5))    # more space for left and bottom axis
test1 <- function(x,sx=0.3,sz=0.4)
    (pi**sx*sz)*(1.2*exp(-(x[,1]-0.2)^2/sx^2-(x[,2]-0.3)^2/sz^2)+
    0.8*exp(-(x[,1]-0.7)^2/sx^2-(x[,2]-0.8)^2/sz^2))
n <- 100
set.seed(1)
x <- runif(n);
z1 <- runif(n);
y <- test1(cbind(x,z1)) + rnorm(n) * 0.1
a <- gam(y ~ s(x) + s(x,z1))
plotmo(a, do.par=FALSE, type2="contour", caption=caption, col.response=3, smooth.col="indianred",
      func=test1, func.col="indianred", func.lwd=5, func.lty=2, smooth.lwd=3)

plotmo(a, do.par=FALSE, degree1=F, degree2=1, type2="image", ylim=NA)
plot(a, select=1)
plot(a, select=2)
plot(a, select=3)
n<-400
sig<-2
set.seed(1)
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
f0 <- function(x) 2 * sin(pi * x)
f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
f <- f0(x0) + f1(x1) + f2(x2)
e <- rnorm(n, 0, sig)
y <- f + e
test.func <- function(x) f0(x[,1]) + f1(x[,2]) + f2(x[,3])
library(mgcv)
caption <- "test mgcv::gam(y~s(x0,x1,k=12)+s(x2)+s(x3,k=20,fx=20)) (and extra persp plot)"
dopar(3,3,caption)
a <- gam(y~s(x0,x1,k=12)+s(x2)+s(x3,k=20,fx=20))
plot(a, select=2)
plot(a, select=3)
plot(a, select=1)
plotmo(a, do.par=FALSE, type2="contour", caption=caption, xlab=NULL, main="", func=test.func, ngrid2=10, drawlabels=FALSE)
plotmo(a, do.par=FALSE, degree1=F, degree2=1, persp.the=-35)

set.seed(1)
caption <- "test plot.gam, with mgcv::gam(doy~s(wind)+s(humidity,wind)+s(vh)+temp,data=ozone1)"
dopar(3,3,caption)
a <- gam(doy ~ s(wind) + s(humidity,wind) + s(vh) + temp, data=ozone1)
plotmo(a, do.par=FALSE, caption=caption, type2="contour", degree1=c("wind","vh"), swapxy=T, xlab=NULL, main="", clip=FALSE)
plot(a, select=1)
plot(a, select=3)
plot(a, select=2)
plot(a, select=4)

detach("package:mgcv")
library(gam)
caption <- "test gam:gam(Ozone^(1/3)~lo(Solar.R)+lo(Wind, Temp),data=airquality)"
set.seed(1)
dopar(3,2,caption)
data(airquality)
airquality <- na.omit(airquality)   # plotmo doesn't know how to deal with NAs yet
a <- gam(Ozone^(1/3) ~ lo(Solar.R) + lo(Wind, Temp), data = airquality)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA, col.response=3)
# termplot gives fishy looking wind plot, plotmo looks ok
# termplot(a) #TODO this fails with R2.5: dim(data) <- dim: attempt to set an attribute on NULL
detach("package:gam")

library(mda)
caption <- "test mars and earth (expect not a close match)"
dopar(6,3,caption)
a <- mars( ozone1[, -1], ozone1[,1], degree=2)
b <- earth(ozone1[, -1], ozone1[,1], degree=2)
# this also tests trace=2 on a non formula model
plotmo(a, do.par=FALSE, caption=caption, trace=2)
plotmo(b, do.par=FALSE)

caption <- "test mars and mars.to.earth(mars) (expect no degree2 for mars)"
dopar(6,3,caption)
a <- mars(ozone1[, -1], ozone1[,1], degree=2)
b <- mars.to.earth(a)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA)
plotmo(b, do.par=FALSE, ylim=NA)

# check fix for bug reported by Martin Maechler:
# form <- Volume ~ .; a <- earth(form, data = trees); plotmo(a) fails

dopar(4,4, "test f <- O3 ~ .; a <- earth(f, data=ozone1)")
fa <- log(O3) ~ .
a <- earth(fa, data=ozone1, degree=2)
print(summary(a))
plot(a, do.par=FALSE)
plotmo(a, do.par=FALSE, degree1=2:3, degree2=c(1,2), col.response = "pink", smooth.col="indianred")
a <- lm(log(doy) ~ I(vh*wind) + I(humidity*temp) + log(ibh), data=ozone1)
plotmo(a, do.par=FALSE, degree1=1:2)
fa <- log(doy) ~ I(vh*wind) + I(humidity*temp) + log(ibh)
a <- lm(fa, data=ozone1)
plotmo(a, do.par=FALSE, degree1=1:2)

# test inverse.func and func

caption <- "test inverse.func=exp"
a <- lm(log(Volume) ~ Girth + Height + I(Girth*Height), data=trees)
my.func <- function(x) -60 + 5 * x[,1] + x[,2] / 3
plotmo(a, caption=caption, inverse.func = exp, col.response = "pink", func=my.func, func.col="gray", ngrid1=1000, type2="p", smooth.col="indianred")

# se testing

caption = "level=.95, lm(doy~., data=ozone1) versus termplot"
dopar(6,3,caption)
a <- lm(doy~., data=ozone1)
plotmo(a, level=.95, do.par=FALSE, caption=caption)
termplot(a, se=2)

caption <- "test different se options, level=.95, lm(log(doy)~vh+wind+log(humidity),data=ozone1)"
dopar(4,3,caption)
a <- lm(log(doy) ~ vh + wind + log(humidity), data=ozone1)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA, level=.95)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA, level=.95, level.shade="pink", level.shade2=3)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA, level=.95, level.shade=3)
plotmo(a, do.par=FALSE, caption=caption, ylim=NULL, level=.95, level.shade=3)

caption <- "test level=.95, lm(log(doy)~vh+wind+log(humidity),data=ozone1)"
dopar(2,3,caption)
a <- lm(log(doy) ~ vh + wind + log(humidity), data=ozone1)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA, level=.95)
termplot(a, se=2)

caption <- "test level=.95 and inverse.func, lm(log(doy)~vh+wind+log(humidity),data=ozone1)"
dopar(3,3,caption)
a <- lm(log(doy) ~ vh + wind + log(humidity), data=ozone1)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA, level=.95)
plotmo(a, do.par=FALSE, caption=caption, ylim=NULL, level=.95, inverse.func=exp)
termplot(a, se=2)

caption <- "test level=.95, glm(lot2~log(u),data=clotting,family=Gamma)"
set.seed(1)
dopar(2,2,caption)
u = c(5,10,15,20,30,40,60,80,100)
lota = c(118,58,42,35,27,25,21,19,18)
clotting <- data.frame(u = u, lota = lota)
a <- glm(lota ~ log(u), data=clotting, family=Gamma)
plotmo(a, do.par=FALSE, caption=caption, col.response=4, pt.pch=7, clip=FALSE, nrug=-1, level=.95, smooth.col="indianred")
termplot(a, se=2)

if(length(grep("package:gam", search())))
    detach("package:gam")
library(mgcv)
set.seed(1)
caption <- "test level=.95, plot.gam, with mgcv::gam(y ~ s(x) + s(x,z1)) with response and func (and extra image plot)"
dopar(3,2,caption)
par(mar = c(3, 5, 1.7, 0.5))    # more space for left and bottom axis
test1 <- function(x,sx=0.3,sz=0.4)
    (pi**sx*sz)*(1.2*exp(-(x[,1]-0.2)^2/sx^2-(x[,2]-0.3)^2/sz^2)+
    0.8*exp(-(x[,1]-0.7)^2/sx^2-(x[,2]-0.8)^2/sz^2))
n <- 100
set.seed(1)
x <- runif(n);
z1 <- runif(n);
y <- test1(cbind(x,z1)) + rnorm(n) * 0.1
a <- gam(y ~ s(x) + s(x,z1))
plotmo(a, do.par=FALSE, type2="contour", caption=caption, col.response=3, func=test1, func.col="magenta", level=.95)
plotmo(a, do.par=FALSE, degree1=F, degree2=1, type2="image", image.col=topo.colors(10),
        ylim=NA, level=.95, main="topo.colors")
plot(a, select=1)
plot(a, select=2)
plot(a, select=3)

# TODO Following commented out because it causes:
#      Error: gam objects in the "gam" package do not support confidence intervals on new data
# detach("package:mgcv")
# library(gam)
# set.seed(1)
# caption <- "test level=.95, gam:gam(Ozone^(1/3)~lo(Solar.R)+lo(Wind, Temp),data=airquality)"
# dopar(3,2,caption)
# data(airquality)
# airquality <- na.omit(airquality)   # plotmo doesn't know how to deal with NAs yet
# a <- gam(Ozone^(1/3) ~ lo(Solar.R) + lo(Wind, Temp), data = airquality)
# plotmo(a, do.par=FALSE, caption=caption, ylim=NA, col.response=3, level=.95)
# # termplot(a)  #TODO this fails with R2.5: dim(data) <- dim: attempt to set an attribute on NULL
# detach("package:gam")

# test factors by changing wind to a factor

ozone2 <- ozone1
ozone2[,"wind"] <- factor(ozone2[,"wind"], labels=c(
    "wind0", "wind2", "wind3", "wind4", "wind5", "wind6",
    "wind7", "wind8", "wind9", "wind10", "wind11"))

# commented out because factors are not yet supported by plotmo.earth
# caption <- "test wind=factor, earth(O3 ~ ., data=ozone2)"
# a <- earth(doy ~ ., data=ozone2)
# set.seed(1)
# dopar(4,3,caption)
# plotmo(a, col.response="gray", level=.95, nrug=-1, do.par=FALSE, caption=caption)
# termplot(a)

caption <- "test wind=factor, lm(doy ~ vh + wind + I(humidity*temp) + log(ibh), data=ozone2)"
a <- lm(doy ~ vh + wind + I(humidity*temp) + log(ibh), data=ozone2)
set.seed(1)
dopar(4,3,caption)
plotmo(a, col.response="gray", level=.95, nrug=-1, do.par=FALSE, caption=caption, smooth.col="indianred")
termplot(a, se=2)

caption <- "test level options"
dopar(2,2,caption)
plotmo(a, do.par=FALSE, degree1=2, degree2=FALSE, level=.95, level.shade=0, caption=caption)
plotmo(a, do.par=FALSE, degree1=2, degree2=FALSE, level=.95, level.shade="orange")
plotmo(a, do.par=FALSE, degree1=2, degree2=FALSE, level=.95, level.shade2=0)

caption <- "test wind=factor, glm(y ~ i + j, family=poisson())"
y <- c(18,17,15,20,10,20,25,13,12)
i <- gl(3,1,9)
j <- gl(3,3)
a <- glm(y ~ i + j, family=poisson())
set.seed(1)
dopar(2,2,caption)
plotmo(a, do.par=F, level=.95, nrug=-1, caption=caption)
termplot(a, se=1, rug=T)

if(length(grep("package:gam", search())))
   detach("package:gam")
caption <- "test wind=factor, gam(doy ~ vh + wind + s(humidity) + s(vh) + temp, data=ozone2)"
library(mgcv)
a <- gam(doy ~ vh + wind + s(humidity) + s(vh) + temp, data=ozone2)
plotmo(a, level=.95, caption=caption)
caption <- "test wind=factor, clip=TRUE, gam(doy ~ vh + wind + s(humidity) + s(vh) + temp, data=ozone2)"
plotmo(a, level=.95, caption=caption, clip=FALSE)
# termplot doesn't work here so code commented out
# dopar(3,3,caption)
# plotmo(a, do.par=FALSE)
# termplot(a)

# test lda and qda, and also col.response, pt.pch, and jitter
library(MASS)
etitanic2 <- etitanic
etitanic2$pclass <- as.numeric(etitanic$pclass)
etitanic2$sex <- as.numeric(etitanic$sex)
etitanic2$sibsp <- NULL
etitanic2$parch <- NULL
lda.model <- lda(survived ~ ., data=etitanic2)
set.seed(7)
plotmo(lda.model, caption="lda", clip=F,
       col.response=as.numeric(etitanic2$survived)+2, type="posterior", nresponse=1, smooth.col="indianred",
       all2=TRUE, type2="image")
set.seed(8)
plotmo(lda.model, caption="lda with no jitter", clip=F,
       col.response=as.numeric(etitanic2$survived)+2, type="posterior", nresponse=1,
       all2=TRUE, type2="image", jitter=0)
qda.model <- qda(survived ~ ., data=etitanic2)
set.seed(9)
plotmo(qda.model, caption="qda", clip=F,
       col.response=as.numeric(etitanic2$survived)+2, type="post", nresponse=2, smooth.col="indianred",
       all2=TRUE, type2="image", jitter.resp=.6, pch.resp=20)

# test plotmo.y from the 2nd argument of the model function (non-formula interface)
lcush <- data.frame(Type=as.numeric(Cushings$Type), log(Cushings[,1:2]))[1:21,]
a <- qda(lcush[,2:3], lcush[,1])
plotmo(a, type="class", all2=TRUE,
       caption= "plotmo.y from 2nd argument of call (qda)",
       type2="contour", ngrid2=100, nlevels=2, drawlabels=FALSE,
       col.response=as.numeric(lcush$Type)+1,
       pt.pch=as.character(lcush$Type))

# # example from MASS (works, but removed because unnecessary test)
# predplot <- function(object, main="", len = 100, ...)
# {
#     plot(Cushings[,1], Cushings[,2], log="xy", type="n",
#          xlab = "Tetrahydrocortisone", ylab = "Pregnanetriol", main = main)
#     for(il in 1:4) {
#         set <- Cushings$Type==levels(Cushings$Type)[il]
#         text(Cushings[set, 1], Cushings[set, 2],
#              labels=as.character(Cushings$Type[set]), col = 2 + il) }
#     xp <- seq(0.6, 4.0, length=len)
#     yp <- seq(-3.25, 2.45, length=len)
#     cushT <- expand.grid(Tetrahydrocortisone = xp,
#                          Pregnanetriol = yp)
#     Z <- predict(object, cushT, ...); zp <- as.numeric(Z$class)
#     zp <- Z$post[,3] - pmax(Z$post[,2], Z$post[,1])
#     contour(exp(xp), exp(yp), matrix(zp, len),
#             add = TRUE, levels = 0, labex = 0)
#     zp <- Z$post[,1] - pmax(Z$post[,2], Z$post[,3])
#     contour(exp(xp), exp(yp), matrix(zp, len),
#             add = TRUE, levels = 0, labex = 0)
#     invisible()
# }
# par(mfrow=c(2,2))
# cush <- log(as.matrix(Cushings[, -3]))
# tp <- Cushings$Type[1:21, drop = TRUE]
# set.seed(203)
# cush.data <- data.frame(tp, cush[1:21,])
# a <- qda(tp~., data=cush.data)
# predplot(a, "QDA example from MASS")
# plotmo(a, type="class", all2=TRUE, type2="contour", degree1=NA, do.par=FALSE,
#        col.response=as.numeric(cush.data$tp)+1)
# plotmo(a, type="class", all2=TRUE, type2="contour", degree1=NA, do.par=FALSE,
#        col.response=as.numeric(cush.data$tp)+1, drawlabels=F, nlevels=2)
# plotmo(a, type="class", all2=TRUE, type2="contour", degree1=NA, do.par=FALSE,
#        col.response=as.numeric(cush.data$tp)+1, drawlabels=F, nlevels=2, ngrid2=100)
# par(mfrow=c(1,1))

library(rpart.plot)
data(kyphosis)
# kyphosis data, earth model
a <- earth(Kyphosis ~ ., data=kyphosis, degree=2, glm=list(family=binomial))
par(mfrow=c(3, 3))
old.mar <- par(mar=c(3, 3, 2, .5))  # small margins to pack figs in
set.seed(9) # for jitter
plotmo(a, do.par=F, type2="image",
       col.response=ifelse(kyphosis$Kyphosis=="present", "red", "lightblue"),
       clip=F)
plotmo(a, do.par=F, clip=F, degree1=0)
par(mar=old.mar)

# kyphosis data, rpart models (also test ngrid2)
fit1 <- rpart(Kyphosis ~ ., data=kyphosis)
par(mfrow=c(3, 3))
old.par <- par(mar=c(.5, 0.5, 2, .5), mgp = c(1.6, 0.6, 0))  # b l t r small margins to pack figs in
prp(fit1, main="rpart kyphosis\nno prior")
plotmo(fit1, degree1=NA, do.par=F, main="", persp.theta=220, nresponse=2)
par(mar=c(4, 4, 2, .5))
plotmo(fit1, nresp=2, degree1=FALSE, do.par=F, main="", type2="image", # test default type="prob"
       col.response=ifelse(kyphosis$Kyphosis=="present", "red", "lightblue"),
       pt.pch=ifelse(kyphosis$Kyphosis=="present", "p", "a"),
       image.col=gray(10:4/10), ngrid2=30)
par(mar=c(.5, 0.5, 2, .5))  # b l t r small margins to pack figs in
plotmo(fit1, type="class", degree1=NA, do.par=F, main="type=\"class\"")
# with type="prob" and response has two columns,
# nresponse should automatically default to column 2
plotmo(fit1, type="prob", degree1=0, do.par=F, main="type=\"prob\"",
       clip=F, ngrid2=50, persp.border=NA, trace=1)
plotmo(fit1, type="prob", nresp=2, degree1=NA, do.par=F, main="", type2="image",
       col.response=ifelse(kyphosis$Kyphosis=="present", "red", "lightblue"),
       pt.pch=20, image.col=gray(10:4/10), ngrid2=5)
# better rpart model with prior
fit2 <- rpart(Kyphosis ~ ., data=kyphosis, parms=list(prior=c(.65,.35)))
prp(fit2, main="rpart kyphosis\nwith prior, better model")
plotmo(fit2, type="v", degree1=NA, do.par=F, main="", persp.theta=220, ngrid2=10)
par(mar=c(4, 4, 2, .5))
plotmo(fit2, type="v", degree1=NA, do.par=F, main="", type2="image",
       col.response=ifelse(kyphosis$Kyphosis=="present", "red", "lightblue"),
       pt.pch=20, image.col=gray(10:4/10), ngrid2=100)
par(mar=old.par$mar, mgp=old.par$mgp)

plotmo(fit1, type="prob", nresponse=1, persp.border=NA, persp.col="pink", all1=TRUE, all2=TRUE,
       caption="plotmo rpart fit1, all1=TRUE, all2=TRUE")
expect.err(try(plotmo(fit1, type="none.such1")))

# rpart model with ozone data
data(ozone1)
par(mfrow=c(4,4))
old.par <- par(mar=c(.5, 0.5, 2, .5), cex=.6, mgp = c(1.6, 0.6, 0))  # b l t r small margins to pack figs in
a1 <- rpart(O3~temp+humidity, data=ozone1)
prp(a1, main="rpart model with ozone data\n(temp and humidity only)\n")
plotmo(a1, do.par=F, degree1=0, main="rpart", persp.ticktype="detail", persp.nticks=2)
expect.err(try(plotmo(a1, type="class")))
# compare to a linear and earth model
a3 <- lm(O3~temp+humidity, data=ozone1)
plotmo(a3, do.par=F, clip=F, main="lm", degree1=0, all2=TRUE, persp.ticktype="detail", persp.nticks=2)
expect.err(try(plotmo(a3, type="none.such2")))
a <- earth(O3~temp+humidity, data=ozone1, degree=2)
plotmo(a, do.par=F, clip=F, main="earth", degree1=NA, persp.ticktype="detail", persp.nticks=2)
expect.err(try(plotmo(a, type="none.such3")))
expect.err(try(plotmo(a, type=c("abc", "def"))))

# detailed rpart model
par(mfrow=c(3,3))
a1 <- rpart(O3~., data=ozone1)
prp(a1, cex=.9, main="rpart model with full ozone data")
plotmo(a1, type="vector", do.par=F, degree1=NA, persp.ticktype="detail",
       persp.nticks=3, degree2=2:3)
par(mar=old.par$mar, cex=old.par$cex, mgp=old.par$mgp)

plotmo(a1, persp.border=NA, all1=TRUE, all2=TRUE,
       caption="plotmo rpart a1, all1=TRUE, all2=TRUE")

library(tree)
tree1 <- tree(O3~., data=ozone1)
plotmo(tree1)
plotres(tree1)

# test xflip and yflip

par(mfrow=c(4, 4))
par(mgp = c(1.6, 0.6, 0))
par(mar=c(4, 4, 2, .5))

flip.test1 <- rpart(Kyphosis ~ ., data=kyphosis)
plotmo(flip.test1, type="prob", nresp=2, degree1=NA, do.par=F, main="", type2="image",
       col.response=ifelse(kyphosis$Kyphosis=="present", "red", "lightblue"),
       pt.pch=20, image.col=gray(10:4/10))
plotmo(flip.test1, type="prob", nresp=2, degree1=NA, do.par=F, main="xflip", type2="image",
       col.response=ifelse(kyphosis$Kyphosis=="present", "red", "lightblue"),
       pt.pch=20, image.col=gray(10:4/10),
       xflip=T)
plotmo(flip.test1, type="prob", nresp=2, degree1=NA, do.par=F, main="yflip", type2="image",
       col.response=ifelse(kyphosis$Kyphosis=="present", "red", "lightblue"),
       pt.pch=20, image.col=gray(10:4/10),
       yflip=T)
plotmo(flip.test1, type="prob", nresp=2, degree1=NA, do.par=F, main="xflip and yflip", type2="image",
       col.response=ifelse(kyphosis$Kyphosis=="present", "red", "lightblue"),
       pt.pch=20, image.col=gray(10:4/10),
       xflip=T, yflip=T)

flip.test2 <- earth(O3~., data=ozone1, degree=2)
plotmo(flip.test2, degree1=NA, degree2=2, do.par=F, main="", type2="cont")
plotmo(flip.test2, degree1=NA, degree2=2, do.par=F, main="xflip", type2="cont",
       xflip=T)
plotmo(flip.test2, degree1=NA, degree2=2, do.par=F, main="yflip", type2="cont",
       yflip=T)
plotmo(flip.test2, degree1=NA, degree2=2, do.par=F, main="xflip and yflip", type2="cont",
       xflip=T, yflip=T)

cat("Expect warnings: ignoring xflip=TRUE for persp plot\n")
plotmo(flip.test2, degree1=NA, degree2=2, do.par=F, main="xflip and yflip", type2="persp",
       xflip=T, yflip=T)

library(randomForest)
data(ozone1)
set.seed(2015)
rf <- randomForest(age~., data=etitanic, ntree=100)
plotmo(rf, trace=1)
plotres(rf, trace=1)
set.seed(3)
a <- randomForest(O3~., data=ozone1, ntree=20)
plotmo(a, caption="randomForest ozone1")
plotres(a)
set.seed(4)
a <- randomForest(Kyphosis ~ ., data=kyphosis, ntree=5, mtry=2)
plotmo(a, type="prob", nresponse="pre", caption="randomForest kyphosis", ndiscrete=10)
# TODO residuals are in range 0 to 1
plotres(a, type="prob", nresponse="pre", caption="plotres randomForest kyphosis")

# gbm
library(gbm)
library(rpart.plot) # for ptitanic, want data with NAs for testing
data(ptitanic)
ptit <- ptitanic[c(1:10,400:410,600:610),] # small data for fast test
ptit <- ptitanic
ptit$survived <- ptit$survived == "survived"
temp <- ptit$pclass # put pclass at the end so can check ordering of importances
ptit$pclass <- NULL
ptit$pclass <- factor(as.numeric(temp), labels=c("first", "second", "third"))
set.seed(1010)
gbm.model <- gbm(survived~., data=ptit, train.frac=.95, verbose=TRUE,
                 n.trees=30, shrinkage=.1) # small number of trees for fast test

par(mfrow=c(4,4))
par(mar=c(3.5, 3, 2, 0.5))  # small margins and text to pack figs in
par(mgp=c(1.5, .5, 0))      # flatten axis elements
plotmo(gbm.model, persp.ticktype="d", persp.nticks=2, do.par=F,
       degree1=0, degree2=3, main="gbm model", caption="gbm models")
plotmo(gbm.model, type2="im", do.par=F,
       col.response=ptit$survived+2, pt.pch=20, pt.cex=.5)
print(summary(gbm.model))   # will also plot
plotres(gbm.model, col=ptit$survived+2, do.par=FALSE)
par(mfrow=c(1,1))

ozplus <- ozone1
# add more variables so we can test all1 and all2
ozplus$ltemp <- log(ozplus$temp)
ozplus$lhum <- log(ozplus$humidity)
ozplus$ltemphum <- log(ozplus$temp) + log(ozplus$humidity)
ozplus$ibttemp <- ozplus$ibt * ozplus$temp
ozplus <- data.frame(scale(ozplus))
set.seed(2015)
gbm.ozplus <- gbm(O3~., data=ozplus, train.frac=.9,
               n.trees=100, cv.fold=2, shrinkage=.1, interact=3)
# note that we get different RSquareds printed in the trace here
plotres(gbm.ozplus, trace=1, SHOWCALL=TRUE)
plotres(gbm.ozplus, predict.n.trees=42, trace=1, SHOWCALL=TRUE)
plotmo(gbm.ozplus, trace=-1, SHOWCALL=TRUE)
plotmo(gbm.ozplus, trace=-1, all1=TRUE, SHOWCALL=TRUE)
plotmo(gbm.ozplus, trace=-1, all2=TRUE, SHOWCALL=TRUE)

library(caret)
set.seed(2015)
caret.earth.mod <- train(O3~., data=ozone1, method="earth",
                         tuneGrid=data.frame(degree=2, nprune=10))
# TODO pairs are not plotted
plotmo(caret.earth.mod, type="raw", trace=1, SHOWCALL=TRUE)
# but the pairs are plotted here
plotmo(caret.earth.mod$finalModel, trace=1, SHOWCALL=TRUE)
plotres(caret.earth.mod, type="raw", trace=1, SHOWCALL=TRUE)

set.seed(2015)
bag <- bagEarth(O3~., data=ozone1, degree=2, B=3)
print(bag$fit)
# pairs are plotted correctly (I think)
plotmo(bag, type="response", trace=1, SHOWCALL=TRUE)

set.seed(2015)
a.bag1 <- bagEarth(trees[,-3], trees[,3], degree=2, B = 3)
plotres(a.bag1, trace=1, SHOWCALL=TRUE)
plotmo(a.bag1, trace=1, SHOWCALL=TRUE, all2=TRUE, caption="bagEarth, trees")

# TODO following doesn't work properly, factors are plotted as continuous?
a.bag3 <- bagEarth(survived~., data=etitanic, degree=2, B=3)
plotmo(a.bag3, clip=F, caption="bagEarth, etitanic", trace=1, SHOWCALL=TRUE)
plotres(a.bag3, clip=F, trace=1, SHOWCALL=TRUE)

# example by Max Kuhn on stackoverflow
set.seed(2015)
etit <- etitanic
etit$survived <- factor(ifelse(etit$survived == 1, "yes", "no"),
                       levels = c("yes", "no"))
# TODO pairs are not plotted
caret.earth.mod2 <- train(survived ~ .,
            data = etit,
            method = "earth",
            tuneGrid = data.frame(degree = 2, nprune = 9),
            trControl = trainControl(method = "none",
                                     classProbs = TRUE))
plotmo(caret.earth.mod2, type="raw", trace=1, SHOWCALL=TRUE)
plotres(caret.earth.mod2, type="raw", trace=1, SHOWCALL=TRUE)

data(ozone1)
a <- train(O3 ~ ., data = ozone1,  method = "earth",
            tuneGrid = data.frame(degree = 2, nprune = 14))
plotmo(a, type="raw", trace=1, SHOWCALL=TRUE)
plotres(a, type="raw", trace=1, SHOWCALL=TRUE)

library(nnet)
data(iris3)
set.seed(301)
samp <- c(sample(1:50,25), sample(51:100,25), sample(101:150,25))
ird <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]),
                  species=factor(c(rep("seto",50), rep("vers", 50), rep("virg", 50))))
ir.nn2 <- nnet(species ~ ., data = ird, subset = samp, size = 2, rang = 0.1,
               decay = 5e-4, maxit = 20)
plotmo(ir.nn2, nresponse=1, type="class", all2=T, degree2=2:6)
plotmo(ir.nn2, nresponse=2, clip=F, all2=T, degree2=1:5)
plotres(ir.nn2, nresponse=2)

#--- fda ------------------------------------------------------------------------------

par(mfrow=c(1,1))

par(mfrow=c(4,5))
par(mar = c(3, 2, 3, .1)) # b, l, t, r
par(mgp = c(1.5, .5, 0))
fda.earth <- fda(Species~., data=iris, keep.fitted=TRUE, method=earth, keepxy=TRUE)
fda.polyreg <- fda(Species~., data=iris, keep.fitted=TRUE, keepxy=TRUE)
fda.bruto <- fda(Species~., data=iris, keep.fitted=TRUE, method=bruto)

# 'fda.polyreg$fit' does not have a 'call' field or 'x' and 'y' fields
expect.err(try(plotmo(fda.polyreg$fit, type="variates", nresponse=1, clip=F, do.par=F)))

plot(1, main="plotmo with fda", xaxt="n", yaxt="n", xlab="", ylab="",
     type="n", bty="n", cex.main=1.2, xpd=NA)

plotmo(fda.earth, type="variates", nresponse=1, clip=F, do.par=F)

plot(1, main="plotmo with fda.earth$fit", xaxt="n", yaxt="n", xlab="", ylab="",
     type="n", bty="n", cex.main=1.2, xpd=NA)

plotmo(fda.earth$fit, nresponse=1, clip=F, do.par=F)

plot(1, main="", xaxt="n", yaxt="n", xlab="", ylab="",
     type="n", bty="n", cex.main=1.5, xpd=NA)

plot(fda.earth)
plotmo(fda.earth, clip=F, do.par=F) # default type is class

plot(fda.polyreg)
plotmo(fda.polyreg, type="variates", nresponse=1, clip=F, do.par=F, degree1=c(1,3,4))
plot(1, main="", xaxt="n", yaxt="n", xlab="", ylab="",
     type="n", bty="n", cex.main=1.5, xpd=NA)

par(mfrow=c(3,3))
par(mar = c(3, 2, 3, .1)) # b, l, t, r
par(mgp = c(1.5, .5, 0))
plot(fda.bruto)
plotmo(fda.bruto, type="variates", nresponse=1, do.par=F)

# neural net package
# for speed we use artificial data because neuralnet is very slow on say trees
library(neuralnet)
n <- 20
set.seed(3)
x1 <- runif(n, min=-1, max=1)
x2 <- runif(n, min=-1, max=1) # x2 is noise
y <- x1^2
data <- data.frame(y=y, x1=x1, x2=x2)
colnames(data) <- c("y","x1", "x2")
set.seed(3)
nn <- neuralnet(y~x1+x2, data=data, hidden=3, rep=3)
print(head(plotmo:::predict.nn(nn, rep="best", trace=TRUE)))
plotmo(nn, trace=1, col.response=2, all2=TRUE, SHOWCALL=TRUE)
plotmo(nn, trace=1, col.response=2, predict.rep="best", SHOWCALL=TRUE)
plotres(nn, trace=1, info=TRUE, SHOWCALL=TRUE)
plotres(nn, trace=1, info=TRUE, predict.rep="best", SHOWCALL=TRUE)

library(biglm)
data(trees)
ff <- log(Volume)~log(Girth)+log(Height)
chunk1 <- trees[1:20,]
chunk2 <- trees[20:31,]
biglm <- biglm(ff,chunk1)
biglm <- update(biglm, chunk2)
plotmo(biglm, pt.col=2, SHOWCALL=TRUE)
plotres(biglm, SHOWCALL=TRUE)

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}
