plot.proftools_profData <- function(x,
                                    type = c("call", "tree", "flame", "time"),
                                    ...) {
    type <- match.arg(type)
    switch(type,
           call = plotProfileCallGraph(x, ...),
           tree = calleeTreeMap(x, ...),
           flame = flameGraph(x, order = "hot", ...),
           time = flameGraph(x, order = "time", ...))
}
