plot.cnv <-
function (x, ...)
{
    if (is.null(attr(x, "meanRatio")))
        plot.cnv.probabilities(x, ...)
    else plot.cnv.intensities(x, ...)
    invisible()
}
