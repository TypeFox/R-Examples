# test.fac.R: test factor plotting in plotmo. This also tests swapxy, xflip, and yflip
# Stephen Milborrow, Berea Mar 2011

library(plotmo)
library(earth)
library(rpart)
data(ozone1)
data(etitanic)
options(warn=1) # print warnings as they occur
if(!interactive())
    postscript(paper="letter")

cat("==test plotmo with factors==\n")
test.fac.with.rpart <- function(ngrid2=20)
{
    et <- etitanic

    col.response <- as.numeric(et$sex)+2
    et$pclass.fac <- et$pclass
    et$parch.num <- et$parch
    parch.fac <- et$parch
    parch.fac[parch.fac >= 3] <- 3
    # use non alphabetically sorted factor levels
    et$parch.fac <- factor(parch.fac, labels=c( "levz", "lev1", "lev2", "levf"))
    et$pclass.num <- as.numeric(et$pclass)
    et$pclass <- et$sex <- et$age <- et$sibsp <- et$parch <- NULL
    cat("names(et):", names(et), "\n") # survived pclass.fac parch.num parch.fac pclass.num

    old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
    par(mfrow=c(4,5))
    par(mar = c(2, 2, 3, 0.5), cex=.6)

    # numeric x numeric
    a2 <- rpart(survived ~ pclass.num+parch.num, data=et)
    set.seed(145)
    plotmo(a2, do.par=F, type2="im", degree1=2,
           col.response=col.response, pt.cex=.3)
    set.seed(145)
    plotmo(a2, do.par=F, type2="con", degree1=NA,
           col.response=col.response, pt.cex=.3)
    set.seed(145)
    plotmo(a2, do.par=F, type2="persp", degree1=NA,
           ngrid2=40, persp.theta=NA, persp.ticktype="d", cex.lab=.8, persp.ntick=2)

    # factor x numeric
    a3 <- rpart(survived ~ pclass.fac+parch.num, data=et)
    set.seed(145)
    plotmo(a3, do.par=F, type2="im",
           col.response=col.response, pt.cex=.3)
    set.seed(145)
    plotmo(a3, do.par=F, type2="con", degree1=NA, ylim=c(0,1),
           col.response=col.response, pt.cex=.3)
    set.seed(145)
    plotmo(a3, do.par=F, type2="persp", degree1=NA,
           ngrid2=40, persp.theta=NA, persp.ticktype="d", persp.border=NA, cex.lab=.8, persp.ntick=2)

    # numeric x factor
    a4 <- rpart(survived ~ pclass.num+parch.fac, data=et)
    set.seed(145)
    plotmo(a4, do.par=F, type2="im", tra=1,
           col.response=col.response, pt.cex=.3)
    set.seed(145)
    plotmo(a4, do.par=F, type2="con", degree1=NA,
           col.response=col.response, pt.cex=.3)
    set.seed(145)
    plotmo(a4, do.par=F, type2="persp", degree1=NA,
           ngrid2=40, persp.theta=NA, persp.ticktype="d", persp.border=NA, cex.lab=.8, persp.ntick=2)

    # factor x factor
    a5 <- rpart(survived ~ pclass.fac+parch.fac, data=et)
    set.seed(145)
    plotmo(a5, do.par=F, type2="im", nrug=100,
           col.response=col.response, pt.cex=.3)
    set.seed(145)
    plotmo(a5, do.par=F, type2="con", degree1=NA,
           col.response=col.response, pt.cex=.3)
    set.seed(145)
    plotmo(a5, do.par=F, type2="persp", degree1=NA,
           ngrid2=40, persp.theta=NA, persp.ticktype="d", persp.border=NA, cex.lab=.8, persp.ntick=2)

    # test ndiscrete
    par(mfrow=c(3,5))
    par(mar = c(2, 2, 3, 0.5), cex=.6)

    plotmo(a2, do.par=F, type2="persp", degree1=2, ndiscrete=0, main="ndiscrete=0",
           persp.theta=NA, persp.ticktype="d", persp.ntick=2,
           col.response=col.response, pt.cex=.3)
    plotmo(a2, do.par=F, type2="im", degree1=NA, ndiscrete=0)
    plotmo(a2, do.par=F, type2="con", degree1=NA, ndiscrete=0)
    plotmo(a2, do.par=F, type2="persp", degree1=2, degree2=NA, ndiscrete=0, main="center", center=TRUE,
           col.response=col.response, pt.cex=.3)

    plotmo(a2, do.par=F, type2="persp", degree1=2, ndiscrete=3, main="ndiscrete=3",
           persp.theta=NA, persp.ticktype="d", persp.ntick=2,
           col.response=col.response, pt.cex=.3)
    plotmo(a2, do.par=F, type2="im", degree1=NA, ndiscrete=3)
    plotmo(a2, do.par=F, type2="con", degree1=NA, ndiscrete=3)
    plotmo(a2, do.par=F, type2="persp", degree1=2, degree2=NA, ndiscrete=3, main="center", center=TRUE,
           col.response=col.response, pt.cex=.3)

    plotmo(a2, do.par=F, type2="persp", degree1=2, ndiscrete=10, main="ndiscrete=10",
           persp.theta=NA, persp.ticktype="d", persp.ntick=2,
           col.response=col.response, pt.cex=.3)
    plotmo(a2, do.par=F, type2="im", degree1=NA, ndiscrete=10)
    plotmo(a2, do.par=F, type2="con", degree1=NA, ndiscrete=10)
    plotmo(a2, do.par=F, type2="persp", degree1=2, degree2=NA, ndiscrete=10, main="center", center=TRUE,
           col.response=col.response, pt.cex=.3)
}
test.fac.with.rpart()
cat("==test plotmo swapxy with factors==\n")
test.swapxy.with.rpart <- function(ngrid2=20)
{
    et <- etitanic[c(1:50,300:350,600:650),]

    col.response <- as.numeric(et$sex)+2
    et$pclass.fac <- et$pclass
    et$parch.num <- et$parch
    parch.fac <- et$parch
    parch.fac[parch.fac > 2] <- 2
    # use non alphabetically sorted factor levels
    et$parch.fac <- factor(parch.fac, labels=c("lev.zero", "lev.one", "lev.two.or.more"))
    print(et$parch.fac)
    et$pclass.num <- as.numeric(et$pclass)
    et$pclass <- et$sex <- et$age <- et$sibsp <- et$parch <- NULL
    cat("names(et):", names(et), "\n") # survived pclass.fac parch.num parch.fac pclass.num

    old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
    par(mfrow=c(4,4))
    par(mar = c(2, 3, 5, 0.5), cex=.6)

    # factor x factor
    a5 <- rpart(survived ~ pclass.fac+parch.fac, data=et)
    for(swapxy in c(F,T)) {
        for(xflip in c(F,T))
            for(yflip in c(F,T)) {
                set.seed(145)
                plotmo(a5, do.par=F, type2="im", degree1=NA,
                       swapxy=swapxy, xflip=xflip, yflip=yflip,
                       main=paste("swapxy=", swapxy, "\nxflip=", xflip, "\nyflip=", yflip),
                       col.response=col.response, pt.cex=3,
                       pt.pch=".")
                set.seed(145)
                plotmo(a5, do.par=F, type2="con", degree1=NA,
                       swapxy=swapxy, xflip=xflip, yflip=yflip,
                       main=paste("swapxy=", swapxy, "\nxflip=", xflip, "\nyflip=", yflip),
                       col.response=col.response, pt.cex=.3)
            }
    }
    par(mfrow=c(2,2))
    set.seed(146)
    plotmo(a5, do.par=F, type2="persp", degree1=NA,
           swapxy=FALSE, main=paste("swapxy=", FALSE),
           ngrid2=40, persp.theta=145, persp.ticktype="d", cex.lab=.8, persp.ntick=5)
    set.seed(146)
    plotmo(a5, do.par=F, type2="persp", degree1=NA,
           swapxy=TRUE, main=paste("swapxy=", TRUE),
           ngrid2=40, persp.theta=145, persp.ticktype="d", cex.lab=.8, persp.ntick=5)
    set.seed(146)
    plotmo(a5, do.par=F, type2="im", degree1=2,
           swapxy=FALSE, main=paste("swapxy=", FALSE))
}
test.swapxy.with.rpart()

aflip <- earth(O3~vh + wind + humidity + temp, data=ozone1, degree=2)
col.response<- ifelse(ozone1$O3 == 38, "red", "pink")

# test xflip arg, degree1 plots
par(mfrow=c(2,2))
set.seed(102)
plotmo(aflip, degree1=1:2, degree2=0, do.par=F, col.response=col.response, nrug=-1, ylab="O3", smooth.col="gray")
plotmo(aflip, degree1=1:2, degree2=F, do.par=F, col.response=col.response, nrug=-1, ylab="O3", xflip=T, main="xflip=TRUE, degree1 plots", , smooth.col="gray")

col.response<- ifelse(ozone1$O3 == 1, "green", "pink")

# test flip args, type2=persp
par(mfrow=c(2,2))
plotmo(aflip, degree1=0, degree2=2, do.par=F, persp.ticktype="d")
plotmo(aflip, degree1=0, degree2=2, do.par=F, persp.tickt="d", swapxy=T, main="swapxy=TRUE")
plot(0, 0, type="n", axes=FALSE, xlab="", ylab="")
plot(0, 0, type="n", axes=FALSE, xlab="", ylab="")

# test swapxy args, type2=image
par(mfrow=c(3,3))

plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="im", col.response=col.response, main="test swapxy on image plots\nreference plot")
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="im", col.response=col.response, swapxy=T, main="swapxy=T")
plot(0, 0, type="n", axes=FALSE, xlab="", ylab="")

plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="im", col.response=col.response, xflip=T, main="xflip=T")
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="im", col.response=col.response, yflip=T, main="yflip=T")
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="im", col.response=col.response, xflip=T, yflip=T, main="xflip=T, yflip=T")

plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="im", col.response=col.response, swapxy=T, xflip=T, main="swapxy=T, xflip=T")
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="im", col.response=col.response, swapxy=T, yflip=T, main="swapxy=T, yflip=T")
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="im", col.response=col.response, swapxy=T, xflip=T, yflip=T, main="swapxy=T, xflip=T, yflip=T")

# test flip args, type2=contour
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="con", col.response=col.response, main="test flip on contour plots\nreference plot")
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="con", col.response=col.response, swapxy=T)
plot(0, 0, type="n", axes=FALSE, xlab="", ylab="")

plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="con", col.response=col.response, xflip=T)
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="con", col.response=col.response, yflip=T)
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="con", col.response=col.response, xflip=T, yflip=T)

plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="con", col.response=col.response, swapxy=T, xflip=T)
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="con", col.response=col.response, swapxy=T, yflip=T)
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="con", col.response=col.response, swapxy=T, xflip=T, yflip=T)

# ordered factor

cat("==test plotmo with ordered factor==\n")
par(mfcol=c(2,2))
par(mar=c(3, 3, 3, 1))
par(mgp=c(1.5, .5, 0))
a <- lm(height~., data=Loblolly)
termplot(a, partial.resid=T, rug=T, terms=2, main="Seed is an ordered factor") # compare to termplot
plotmo(a, do.par=F, col.resp="gray", nrug=T, all2=T)

#---------------------------------------------------------------------------
# test ndiscrete with integer and non integer predictors, with missing values

par(mfcol=c(2,4))
par(mar=c(3, 3, 3, 1))
par(mgp=c(1.5, .5, 0))
et <- etitanic
et$var <- et$parch
et$var[et$var==1] <- 0 # want a "hole" in var's value, for testing
et$var[1:3] <- 6
cat("table(et$var):")
print(table(et$var))
cat("\n")
a <- earth(survived~var+age, data=et, degree=2, pm="none")

plotmo(a, trace=FALSE, ndiscrete=0,
       main="integral var\n(var levels are 0 2 3 4 5 6)\nndiscrete=0", cex.lab=.8,
       do.par=F, smooth.col="indianred", persp.ticktype="d", clip=F, degree1=0, persp.theta=40)

plotmo(a, ndiscrete=0,
       do.par=F, smooth.col="indianred", ylim=c(-.5,1), degree2=0, degree1=1)

#------------
plotmo(a, ndiscrete=10,  main="integral var\nndiscrete=10", cex.lab=.8,
       do.par=F, smooth.col="indianred", persp.ticktype="d", clip=F, degree1=0, persp.theta=40)

plotmo(a, trace=0, ndiscrete=10,
       do.par=F, smooth.col="indianred", ylim=c(-.5,1), degree2=0, degree1=1)

#------------
et$var <- et$var / 2
cat("table(et$var):")
print(table(et$var))
cat("\n")
a <- earth(survived~var+age, data=et, degree=2, pm="none")

plotmo(a, ndiscrete=0,
       main="integral var\n(var levels are 0 1 1.5 2 2.5 3)\nndiscrete=0", cex.lab=.8,
       do.par=F, smooth.col="indianred", persp.ticktype="d", clip=F, degree1=0, persp.theta=40)

plotmo(a, ndiscrete=0,
       do.par=F, smooth.col="indianred", ylim=c(-.5,1), degree2=0, degree1=1)

#------------
plotmo(a, ndiscrete=10, main="non integral var\nndiscrete=10", cex.lab=.8,
       do.par=F, smooth.col="indianred", persp.ticktype="d", clip=F, degree1=0, persp.theta=40)

plotmo(a, ndiscrete=10,
       do.par=F, smooth.col="indianred", ylim=c(-.5,1), degree2=0, degree1=1)

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}
