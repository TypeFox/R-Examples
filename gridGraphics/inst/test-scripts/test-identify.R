
# Unfortunately, this cannot be run non-interactively (or at least
# I have not thought of a way to do so)

library(gridGraphics)
                             
notrun <- function() {
    plot(1)
    identify(1)
    dl <- recordPlot()
    dev.off()
    plotdiff(expression(replayPlot(dl)), "identify")
}

