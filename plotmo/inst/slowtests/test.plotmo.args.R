# test.plotmo.args..R: test dot and other argument handling in plotmo

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

old.warn <- options("warn")
options(warn=2)
lm.mod <- lm(O3~wind, data=ozone1)

expect.err(try(plotmo(lm.mod, se=2, leve=.95)), "plotmo's 'se' argument is deprecated, please use 'level' instead")

expect.err(try(plotmo(lm.mod, se=T)), "plotmo's 'se' argument is deprecated, please use 'level=.95' instead")

expect.err(try(plotmo(lm.mod, se=.8)), "plotmo's 'se' argument is deprecated, please use 'level=.95' instead")

expect.err(try(plotmo(lm.mod, level=2)), "level=2 is out of range, try level=.95")

oz2 <- ozone1[1:40,]
set.seed(2015)
a <- earth(O3~temp+wind, dat=oz2, deg=2, nk=21, ncr=3, nfo=3, varmod.me="lm")

expect.err(try(plotmo(a, lw=2, trace=1, thresh=.9, SHOWCALL=TRUE)), "predict.earth ignored argument 'lw'")

options(warn=old.warn$warn)

# test col.response and friends
plotmo(a, col.response=2, pch.response=c(1, 2, 20), type2="co", SHOWCALL=TRUE) # pch.response tests back compat
plotmo(a, pt.col=c(1,2,3), pt.pch=c(1, 2, 20), type2="im", SHOWCALL=TRUE)
plotmo(a, pt.col=c(1,2,3), pt.pch=paste(1:nrow(oz2)), pt.cex=.8, type2="im", do.par=2, SHOWCALL=TRUE)
plotmo(a, pt.col=c(1,2,3), pt.pch=paste(1:nrow(oz2)), pt.cex=.8, type2="co", degree1=0, do.par=F)
plotmo(a, col=2, SHOWCALL=TRUE) # will cause red response points
plotmo(a, pt.col=4, col=3, persp.col="pink", SHOWCALL=TRUE) # col now goes to lines

# test cex and nrug and smooth
plotmo(a, cex=.8, SHOWCALL=TRUE, nrug=-1, rug.col=2, rug.lwd=1, smooth.col=3,
       bty="n", col.lab="darkorange", xlab="an x label", cex.lab=1.2) # esoteric, but they work
plotmo(a, SHOWCALL=TRUE, density.col=2, density.lty=2, smooth.col=3, smooth.f=.3, col="lightblue")
plotmo(a, cex=1.2, SHOWCALL=TRUE, nrug="density")

# test caption, grid, interval options
plotmo(a, caption.col=3, caption.font=2, grid.col="pink",
       level=.8, SHOWCALL=TRUE)
plotmo(a, caption.col=2, caption.font=2, caption.cex=.8, grid.col=TRUE, bty="n",
       level=.8, level.shade="lightblue", level.shade2="red",
       grid.lty=3, grid.lwd=4, grid.nx=NA, SHOWCALL=TRUE)

# test overall plot args handled by par() and graphics args outside do.par
par(mfrow=c(2,2), mar = c(3,3,3,1), mgp = c(1.5,.5,0), oma=c(0,0,4,0))
par(col.main="#456789")
old.mar <- par("mar")
old.mfcol <- par("mfcol")
cat("before par: cex=", par("cex"), " col.main=", par("col.main"),
    " col.axis=", par("col.axis"), " mar=", par("mar"), " mfcol=", par("mfcol"),
    "\n", sep="")
plotmo(a, mfcol=c(2,3), cex.main=1.4, oma=c(5,5,5,5), SHOWCALL=TRUE)
plotmo(a,          caption="no cex")
plotmo(a, cex=1,   caption="cex=1, plot should be identical to previous page")
plotmo(a, cex=1.2, caption="cex=1.2")
plotmo(a, do.par=FALSE, degree2=0, degree1=1,         main="do.par=FALSE no cex", caption="a test graphics args with do.par=FALSE")
plotmo(a, do.par=FALSE, degree2=0, degree1=1, cex=1,  main="do.par=FALSE cex=1")
plotmo(a, do.par=FALSE, degree2=0, degree1=1, cex=.8, main="do.par=FALSE cex=.8")
plotmo(a, do.par=FALSE, degree2=0, degree1=1, cex=1.1, xlab="xlab", col.main=2, col.axis="blue", col.lab=3, font.lab=2,
          main="do.par=FALSE cex=1.1, col.main=2\ncol.axis=\"blue\", col.lab=3, font.lab=2")
plotmo(a, do.par=FALSE, degree1=1, degree2=1, persp.ticktype="d",
          main="do.par=FALSE persp.ticktype=\"d\"")
# all of these should have been restored
cat("after par: cex=", par("cex"), " col.main=", par("col.main"),
    " col.axis=", par("col.axis"), " mar=", par("mar"), " mfcol=", par("mfcol"),
    "\n", sep="")
stopifnot(par("col.main") == "#456789")
stopifnot(par("mar") == old.mar)
stopifnot(par("mfcol") == old.mfcol)
par(col.main=1)

# test aliasing of col with other args, and back compat of col.degree1 vs degree1.col
data(etitanic)
a20 <- earth(pclass ~ ., data=etitanic, degree=2)
plotmo(a20, nresponse=1, col=2, col.degree1=3, persp.col="pink", SHOWCALL=1, degree1=1:2, degree2=1:2)
plotmo(a20, nresponse=1, lty=2, persp.lty=1, SHOWCALL=1, degree1=1:2, degree2=1:2)

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}
