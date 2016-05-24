# test.dots.plotmo.R: test dots functions with the plotmo and earth libraries

options(warn=1) # print warnings as they occur

if(!interactive())
    postscript(paper="letter")

strip.space <- function(s) gsub("[ \t\n]", "", s)

cat0 <- function(...) cat(..., sep="")

# test that we got an error as expected from a try() call
expect.err <- function(object, expected.msg="")
{
    if(class(object)[1] == "try-error") {
        msg <- attr(object, "condition")$message[1]
        if(length(grep(expected.msg, msg, fixed=TRUE)))
            cat0("Got error as expected from ",
                deparse(substitute(object)), "\n")
        else
            stop(sprintf("Expected: %s\n  Got:      %s",
                         expected.msg, substr(msg, 1, 1000)))
    } else
        stop("did not get expected error ", expected.msg)
}
library(plotmo)

library(earth) # will also load plotmo
data(ozone1)

a <- earth(O3~., data=ozone1, degree=2)
expect.err(try(plotmo(a, persp.s=99)), "'s' matches both the 'sub' and 'scale' arguments of persp()")

# Commented out because we now silently drop partial plot args like cex.l
# expect.err(try(plotmo(a, cex.l=.8, cex.la=.9)), "arguments 'cex.l' and 'cex.la' both match 'cex.lab' in draw.plot.degree1")
# expect.err(try(plotmo(a, persp.shad=1, persp.sh=2)), "'persp.shad' and 'persp.sh' both match the 'shade' argument of persp()")

old.warn <- options("warn")
options(warn=2) # treat warnings as errors

# Commented out because we now silently drop partial plot args like cex.l
# expect.err(try(plotmo(a, cex.l=.8)), "\"cex.l\" is not a graphical parameter")
# expect.err(try(plotmo(a, cex.lxx=.8)), "\"cex.lxx\" is not a graphical parameter")
# expect.err(try(plotmo(a, cex.labx=.8)), "\"cex.labx\" is not a graphical parameter")
# expect.err(try(plotmo(a, cex.l=.8, cex.lab=.9)), "\"cex.l\" is not a graphical parameter")

expect.err(try(plotmo(a, nonesuch=.8)), "predict.earth ignored argument 'nonesuch'")

expect.err(try(plotmo(a, lw=2)), "predict.earth ignored argument 'lw'")

# test main, xlab, ylab, etc. arguments with recycling
a <- earth(O3~., data=ozone1, degree=2)
plotmo(a, caption="test main, xlab, ylab, ticktype arguments",
       main=c("main1", "main2", "main3", "main4"), xlab=c("x1", "x2"),
       persp.nticks=2, persp.ticktype="d", ylab=c("y1", "y2", "y3"))

par(mfrow=c(2,2), mar = c(3,3,3,1), mgp = c(1.5,.5,0), oma=c(0,0,4,0))
plotmo(a, trace=1, do.par=FALSE, degree1=1, degree2=1, caption="top: standard\nbottom: lwd=2 thresh=.9") # no errors or warnings
plotmo(a, lwd=2, trace=1, thresh=.9, do.par=FALSE, degree1=1, degree2=1) # no errors or warnings

options(warn=old.warn$warn)

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}
