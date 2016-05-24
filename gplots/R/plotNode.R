here <- function() {}

plot.dendrogram <- stats:::plot.dendrogram
environment(plot.dendrogram) <- environment(here)

plotNodeLimit <- stats:::plotNodeLimit
environment(plotNodeLimit) <- environment(here)

.memberDend <- stats:::.memberDend
environment(.memberDend) <- environment(here)

.midDend <- stats:::.midDend
environment(.midDend) <- environment(here)

unByteCode <- function(fun)
    {
        FUN <- eval(parse(text=deparse(fun)))
        environment(FUN) <- environment(fun)
        FUN
    }

plotNode <- unByteCode(stats:::plotNode)
environment(plotNode) <- environment(here)
