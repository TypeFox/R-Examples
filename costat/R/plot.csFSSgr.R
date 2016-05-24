plot.csFSSgr <-
function (x, plotclust = TRUE, plotscale = TRUE, sollabels=FALSE, ...) 
{
    if (plotscale == TRUE) {
        if (sollabels==FALSE)
		plot(x$epscale[, 1], x$epscale[, 2], ...)
        else	{
		plot(x$epscale[, 1], x$epscale[, 2], type="n", ...)
		text(x$epscale[, 1], x$epscale[, 2], labels=as.character(x$epclust$labels))
		}
    }
    if (plotclust == TRUE) {
        plot(x$epclust, ...)
    }
}
