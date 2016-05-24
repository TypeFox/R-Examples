
library(playwith)
require(latticeExtra)

## A mini application (a playwith toolbar) for clustering.
## Clustering method and distance method can be chosen.
## The "Cut Tree" button allows clusters to be defined.
## Buttons open linked plots (marginal distibutions,
## parallel plot, or MDS plot) showing clusters.

clusterApp <- function(data, ...)
{
    dataName <- toString(deparse(substitute(data)), 30)

    param <- list()
    param$cut.tree <-
        list(function(playState) {
            foo <- playPointInput(playState)
            if (is.null(foo$coords)) {
                ## click was invalid; remove cut line and clusters
                playState$env$cut.height <- NULL
                playClear(playState, type = "annotations",
                          redraw = TRUE)
                ## (redraw to trigger the updating functions)
            } else {
                ## store cut height and draw the line
                height <- foo$coords$y[1]
                playState$env$cut.height <- height
                annLine <- call("panel.abline", h = height,
                                col = "red", lty = 2)
                playAnnotate(playState, annLine, add = FALSE,
                             redraw = TRUE)
                ## (redraw to trigger the updating functions)
            }
        }, label = "Cut Tree")

    param$clust.method <-
        list(c("complete", "ward", "single", "average",
               "mcquitty", "median", "centroid"),
             label = "Clustering method")
    param$dist.method <-
        list(c("euclidean", "maximum", "manhattan",
               "canberra", "binary"),
             label = "Distance method")

    param$parallel <-
        list(function(playState) {
            refObj <- playState
            playwith(parallel(~ refObj$env$data,
                         groups = refObj$env$clusters,
                         main = "Parallel Coordinates plot"),
                     title = "Parallel plot",
                     new = TRUE,
                     link.to = playState)
        }, label = "Parallel")

    param$marginals <-
        list(function(playState) {
            refObj <- playState
            playwith(marginal.plot(~ refObj$env$data,
                         groups = refObj$env$clusters,
                         main = "Marginal distributions"),
                     main = "Marginal distributions",
                     new = TRUE,
                     link.to = playState)
        }, label = "Marginals")

    param$mds <-
        list(function(playState) {
            refObj <- playState
            playwith({
                mds <- with(refObj$env,
                            cmdscale(dist(data, method = dist.method)))
                plot(mds, main = "Classical Multidimensional Scaling")
                text(mds, rownames(mds), cex = 0.7)},
                     title = "MDS plot",
                     new = TRUE,
                     update.actions = function(playState)
                         drawClusterPoints(playState, refObj),
                     link.to = playState)
        }, label = "MDS Plot")

    ## post-plot (update.action) to store 'clusters' object
    computeClusters <- function(playState) {
        height <- playState$env$cut.height
        oldClust <- playState$env$clusters
        if (is.null(height)) {
            ## no clusters defined
            playState$env$clusters <-
                rep(0, length(xyCoords()$x))
        } else {
            ## store cluster vector
            playState$env$clusters <-
                cutree(as.hclust(callArg(playState, 1)),
                       h = height)
        }
        ## and update any linked plots
        ## (since they refer to this 'clusters' object)
        if (!identical(playState$env$clusters, oldClust))
            updateLinkedSubscribers(playState, redraw = TRUE)
    }

    drawClusterPoints <- function(playState, refObj = playState) {
        ## refObj is the plot environment holding 'clusters'
        clusters <- refObj$env$clusters
        if (is.null(clusters)) return()
        coords <- xyCoords(playState)
        playDo(playState,
               quote(panel.xyplot(coords$x, coords$y,
                  pch = 16, groups = clusters, subscripts = TRUE))
               )
    }

    playwith(plot(hclust(dist(data, method = dist.method),
                         method = clust.method), hang = -1, cex = 0.7),
             title = paste("Cluster App:", dataName),
             parameters = param,
             update.actions = list(computeClusters, drawClusterPoints),
             click.mode = "Brush",
             ...)
}

clusterApp(USArrests)

