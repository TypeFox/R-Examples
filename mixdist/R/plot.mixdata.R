## last modified May 2008

plot.mixdata <- function(x, mixpar = NULL, dist = "norm", root = FALSE, 
    ytop = NULL, clwd = 1, main, sub, xlab, ylab, bty, ...) 
{
    mixdataobj<-x
    plot.mix(mixdataobj, mixpar, dist, root, ytop, clwd, main, 
        sub, xlab, ylab, bty, ...)
    invisible()
}
