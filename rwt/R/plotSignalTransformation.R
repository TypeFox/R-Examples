###
### $Id: plotSignalTransformation.R 22 2014-06-20 20:59:33Z plroebuck $
###

##-----------------------------------------------------------------------------
plotSignalTransformation <- function(x, s, title, col.x="blue", col.s="red") {
    plot(x, main=title, ylim=c(1, -1), type="n", axes=FALSE)
    box()
    axis(1, seq(0, length(x), by=100))
    axis(2, seq(1, -1, by=-0.2))
    lines(x, col=col.x)
    lines(s, col=col.s)
    legend(length(s), 0, legend="signal", text.col=col.s, bty="n", xjust=1)
    invisible()
}

