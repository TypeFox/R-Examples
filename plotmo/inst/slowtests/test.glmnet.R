# test.glmnet.R: glmnet tests for plotmo and plotres

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
printf("library(glmnet)\n")
library(glmnet)

data(ozone1)
data(etitanic)

get.tit <- function() # abbreviated titanic data
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
tit <- get.tit()

set.seed(2015)
xmat <- as.matrix(tit[,c(2,5,6)])
set.seed(2015)
mod.glmnet.xmat <- glmnet(xmat, tit[,4])
# plotmo on glmnet mods is boring but we test it anyway
plotmo1(mod.glmnet.xmat)
plotres1(mod.glmnet.xmat)

set.seed(2015)
mod.cv.glmnet.xmat <- cv.glmnet(xmat, tit[,4], nfolds=3)

# following was needed before plotmo 3.1.3 (before adding plotmo.prolog.cv.glmnet)
# mod.cv.glmnet.xmat$x <- as.data.frame(xmat)
# mod.cv.glmnet.xmat$y <- tit[,4]

cat("==Test plotmo trace=1 and lambda.min\n")
plotmo1(mod.cv.glmnet.xmat, predict.s="lambda.min", trace=1)
cat("==Test plotmo trace=2 and lambda.min\n")
plotmo1(mod.cv.glmnet.xmat, predict.s="lambda.min", trace=2)
cat("==Test plotres trace=1 and lambda.1se\n")
plotres1(mod.cv.glmnet.xmat, predict.s="lambda.1se", trace=1)
cat("==Test plotres trace=2 and lambda.1se\n")
plotres1(mod.cv.glmnet.xmat, predict.s="lambda.1se", trace=2)

set.seed(2015)
x <- matrix(rnorm(100*20),100,20)
y <- rnorm(100)
mod.glmnet.x <- glmnet(x,y)
plotmo1(mod.glmnet.x)

# glmnet with sparse matrices
set.seed(2015)
n <- 100
p <- 20
nzc <- trunc(p/10)
x <- matrix(rnorm(n*p),n,p)
iz <- sample(1:(n*p),size=n*p*.85,replace=FALSE)
x[iz] <- 0
sx <- Matrix(x,sparse=TRUE)
# colnames(sx) <- paste("x", 1:ncol(sx), sep="") # need column names for plotmo
inherits(sx,"sparseMatrix") # confirm that it is sparse
beta <- rnorm(nzc)
fx <- x[,seq(nzc)]%*%beta
eps <- rnorm(n)
y <- fx+eps
px <- exp(fx)
px <- px/(1+px)
ly <- rbinom(n=length(px),prob=px,size=1)
mod.glmnet.sx <- glmnet(sx,y)
plotmo1(mod.glmnet.sx, all2=TRUE)

y <- trees$Volume
x <- as.matrix(data.frame(Girth=trees$Girth, Height=trees$Height))
glmnet <- glmnet(x, y)
par(mfrow=c(2,4), mar=c(3,3,3,1), mgp=c(1.5,0.5,0), oma=c(0,0,2.5,0))
plotres(glmnet, do.par=FALSE, caption="glmnet and lm: top and bottom should be the same")
lm <- lm(Volume~., data=trees)
plotres(lm, do.par=FALSE, SHOWCALL=TRUE)

par(mfrow=c(3,2), mar=c(3,3,3,1), mgp=c(1.5,0.5,0), oma=c(0,0,2.5,0))
plotres(glmnet, do.par=FALSE, which=c(1,3), w1.xvar="norm",
        caption="glmnet with various options", SHOWCALL=TRUE)
plotres(glmnet, trace=1, do.par=FALSE, which=c(1,3), SHOWCALL=TRUE)
plotres(glmnet, trace=1, do.par=FALSE, which=c(1,3), predict.s=5, SHOWCALL=TRUE)

printf("======== glmnet multinomial (multnet)\n")
set.seed(2015)
n <- 20
p <- 3
nresp <- 2
x <- matrix(rnorm(n*p),n,p)
y <- rep(1:nresp, each=n/nresp)
fit3 <- glmnet(x, y, family="multinomial")
# TODO why are all the residuals positive?
plotres(fit3, nresponse=1, caption="why are all the residuals positive?")
# TODO seems to work but prints plotmo_y, why?
plotres(fit3, nresponse=2, SHOWCALL=TRUE)

printf("======== glmnet additional tests\n")
set.seed(2015)
p <- 10
n <- 30
x <- cbind(matrix(rnorm(n*p),n,p))
y <- rowSums(x[,1:3]^3)
glmnet <- glmnet(x,y)
plotres(glmnet, SHOWCALL=TRUE, caption="glmnet: y <- rowSums(x[,1:3]^3)")
plotres(glmnet, SHOWCALL=TRUE, w1.xvar="norm")
par(mfrow=c(1,1))
omar <- par("mar")
ocex.axis <- par("cex.axis")
ocex.lab <- par("cex.lab")
plotres(glmnet, SHOWCALL=TRUE, which=1)
stopifnot(par("mar") == omar)
stopifnot(par("cex.axis") == ocex.axis)
stopifnot(par("cex.lab") == ocex.lab)

# test some args for plot.glmnetx
plotres(glmnet, predict.s=.05, SHOWCALL=TRUE, trace=0, col.main=2,
        w1.xlab="my xlab", w1.ylab="my ylab", w1.main="my main",
        w1.col=terrain.colors(n=5))

old.par <- par(no.readonly=TRUE)
plotres(glmnet, predict.s=.05, SHOWCALL=TRUE, which=c(1,3), grid.col="gray", do.par=2)
plotres(glmnet, predict.s=.05, SHOWCALL=TRUE, which=c(1,3), w1.s.col=0, do.par=0)
par(old.par)

# TODO the following issues a stream of warnings: restarting interrupted promise evaluation
expect.err(try(plotres(glmnet, w1.col=nonesuch)), "cannot evaluate 'col'")

printf("======== glmnet additional tests, multiple response models\n")
set.seed(2015)
x <- cbind((1:n)/n, matrix(rnorm(n*(p-1)),n,p-1))
colnames(x) <- paste0("x", 1:p)
# ymultresp <- cbind(rowSums(x[,1:5]^3), rowSums(x[,5:p]^3), 1:n)
set.seed(1)
ymultresp <- cbind(x[,1]+.001*rnorm(n), rowSums(x[,2:5]^3), rnorm(n))
glmnet.multresp <- glmnet(x, ymultresp, family="mgaussian")

plotres(glmnet.multresp, nresponse=1, SHOWCALL=TRUE, which=c(1:3), do.par=2)
# manually calculate the residuals
plot(x=predict(glmnet.multresp, newx=x, s=0)[,1,1],
     y=ymultresp[,1] - predict(glmnet.multresp, newx=x, s=0)[,1,1],
     pch=20, xlab="Fitted", ylab="Residuals",
     main="Manually calculated residuals, nresponse=1, s=0")
abline(h=0, col="gray")

plotres(glmnet.multresp, nresponse=2, w1.label=5, trace=1, SHOWCALL=TRUE)

graphics::par(mfrow=c(2,2), mgp=c(1.5,0.4,0), tcl=-0.3, cex.main=1,
              font.main=1, mar=c(4,3,1.2,0.8), oma=c(0,0,4,0), cex=0.83)
plotres(glmnet.multresp, nresponse=2, SHOWCALL=TRUE, which=3, do.par=FALSE,
        caption="glmnet.multresp compare to manually calculated residuals")
plot(x=predict(glmnet.multresp, newx=x, s=0)[,2,1],
     y=ymultresp[,2] - predict(glmnet.multresp, newx=x, s=0)[,2,1],
     pch=20, xlab="Fitted", ylab="Residuals",
     main="Manual residuals, nresponse=2, s=0")
abline(h=0, col="gray")

plotres(glmnet.multresp, nresponse=2, predict.s=.5, SHOWCALL=TRUE, which=3, do.par=FALSE)
plot(x=predict(glmnet.multresp, newx=x, s=.5)[,2,1],
     y=ymultresp[,2] - predict(glmnet.multresp, newx=x, s=.5)[,2,1],
     pch=20, xlab="Fitted", ylab="Residuals",
     main="Manual residuals, nresponse=2, s=.5")
abline(h=0, col="gray")

plotres(glmnet.multresp, predict.s=.05, nresponse=3, info=TRUE, SHOWCALL=TRUE) # essentially random

data(trees)
set.seed(2015)
par(mfrow=c(2,3), mar=c(3,3,3,.5), oma=c(0,0,3,0), mgp=c(1.5,0.4,0), tcl=-0.3)
# variable with a long name
x50 <- cbind(trees[,1:2], Girth12345678901234567890=rnorm(nrow(trees)))
mod.with.long.name <- glmnet(data.matrix(x50),data.matrix(trees$Volume))
plotres(mod.with.long.name, which=1, caption="test plot.glmnetx with x50 and x60")

# one inactive variable (all coefs are zero for variable "rand")
# (this requires set.seed(2015) above and formation of x50 above)
x60 <- cbind(trees[,1], rand=rnorm(nrow(trees)), trees[,2])
# complicate the issue: use an unnamed column (column 3)
colnames(x60) <- c("Girth", "rand", "")
mod.with.inactive.var <- glmnet(data.matrix(x60),data.matrix(trees$Volume))
stopifnot(all(mod.with.inactive.var$beta["rand",] == 0))
plotres(mod.with.inactive.var, which=1)
plotres(mod.with.inactive.var, which=1, w1.xvar="norm")
# compare to plot.glmnet (but note that labels aren't always plotted unless par=c(1,1)?)
plot(mod.with.inactive.var, xvar="norm", label=TRUE)
# plotmo calls the unnamed column "x3", fair enough
plotmo(mod.with.inactive.var, do.par=FALSE, pt.col=2)

# single active variable
x70 <- cbind(trees[,1,drop=F], 0)
a <- glmnet(data.matrix(x70), data.matrix(trees$Volume))
par(mfrow=c(2,2), mar=c(3,3,2,4))
plotres(a, which=1, predict.s=1, caption="single active variable")
plotres(a, which=1, w1.xvar="norm")
plotres(a, which=1, w1.xvar="lambda")
plotres(a, which=1, w1.xvar="dev")

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}
