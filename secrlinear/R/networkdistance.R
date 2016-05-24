############################################################################################
## package 'secrlinear'
## networkdistance.R
## last changed 2014-09-05, 2014-09-07, 2014-10-26, 2014-10-30,31, 2014-11-01, 2014-11-04
## 2014-12-08 replace NA with zero in tempdist
## 2014-12-09 tentatively increased lwd for Euclidean path in showpath()
############################################################################################

asgraph <- function (mask) {

    ## "mask" is a provisional mask including line termini that will be used for the graph
    ## but dropped from mask points

    ## determine connection by Euclidean distance between points
    tempdist <- as.matrix(dist(mask))
    spacingfactor <- attr(mask, 'spacingfactor')
    tempdist[tempdist > (spacing(mask) * spacingfactor)] <- NA
    tempdist[upper.tri(tempdist, diag = TRUE)] <- NA

    ## value is the number of edges to create between a pair of vertices
    ## use only the lower triangle of distance matrix for undirected graph
    tempdist[!is.na(tempdist)] <- 1
    tempdist[is.na(tempdist)] <- 0
    graph <- graph.adjacency(tempdist, mode = 'lower')

    ## force within-line spacings to actual separation along network,
    ## known because we made it so with along.line
    ed <- get.edges(graph, E(graph))
    x <- mask$x
    y <- mask$y
    LineID <- covariates(mask)$LineID
    Terminal <- covariates(mask)$Terminal
    if (is.null(Terminal))
        terminal <- FALSE
    else
        terminal <- Terminal[ed[,1]] | Terminal[ed[,2]]
    nonconsecutive <- abs(apply(ed,1,diff)) > 1
    same <- (LineID[ed[,1]] == LineID[ed[,2]])
    lengths <- ((x[ed[,1]]- x[ed[,2]])^2  + (y[ed[,1]]- y[ed[,2]])^2)^0.5
    lengths[same & !terminal] <- spacing(mask)

    ## add vertex attributes - these really are needed because termini involved
    V(graph)$LineID <- LineID
    V(graph)$x <- x
    V(graph)$y <- y

    ## add edge attribute (required)
    E(graph)$weight <- lengths

    ## drop edges between nonconsecutive points within a line
    notOK <- same & nonconsecutive & !terminal
    graph[from = ed[notOK,1], to = ed[notOK,2]] <- NULL
    graph
}

networkdistance <- function (xy1, xy2, geometry) {
    ## notify secr.fit that no parameter is required
    if (missing(xy1)) return(character(0))
    ## form matrix of distances from each point in xy1 (rows)
    ## to each point in xy2 (columns), given geometry 'geometry'
    rn <- list(rownames(xy1), rownames(xy2))
    ## in case xy is not matrix or dataframe 2014-09-05
    if (is.null(dim(xy1)))
        xy1 <- matrix(xy1, ncol = 2)
    if (is.null(dim(xy2)))
        xy2 <- matrix(xy2, ncol = 2)
    if (missing(geometry))
        geometry <- xy2
    ## obtain igraph representation
    gr <- attr(geometry, 'graph')
    if (is.null(gr)) gr <- asgraph(geometry)

    ## relate points to geometry
    geometryxy <- cbind(V(gr)$x, V(gr)$y)
    matchedxy1 <- nearesttrap(xy1, geometryxy)
    matchedxy2 <- nearesttrap(xy2, geometryxy)
    uniquematchedxy1 <- unique(matchedxy1)
    ## weights = NULL means "use edge attribute 'weight' if present"
    tmp <- t(igraph::shortest.paths(gr, weights = NULL,
                v =  matchedxy2, to = uniquematchedxy1))
    rematch <- match(matchedxy1, uniquematchedxy1)
    tmp <- tmp[rematch,,drop = FALSE]
    if (!("weight" %in% igraph::list.edge.attributes(gr)))
        tmp <- tmp * spacing(geometry)
    dimnames(tmp) <- rn
    tmp
}
## Call this function interactively to verify the difference between
## Euclidean and network distances. In order to compute network distances
## with networkdistance() we set the `mask' attribute of its second
## argument.

showpath <- function (mask, add = FALSE, ...) {
    plot(mask, cex = 0.5, add = add)
    results <- NULL
    paths <- vector('list')
    npath <- 0
    repeat {
        cat("Click near 'from' and 'to' points\n\n")
        flush.console()
        xy <- as.data.frame(locator(2))
        if (nrow(xy) < 2)
            break
        else {
            matched <- nearesttrap(xy, mask)
            points(xy, pch = 1, cex = 1.3, col = 'black', type = 'p')
            points(mask[matched,], pch = 16, cex = 1.3, col = 'red',
                   type = 'o', lwd = 2, lty = 2)
            path <- get.shortest.paths(asgraph(mask), matched[1], matched[2])
            maskpoints <- data.frame(mask[path$vpath[[1]],])
            lines(maskpoints, col = 'red', type = 'l', ...)
            npath <- npath + 1
            paths[[npath]] <- maskpoints
            xy1 <- xy[1,]; xy2 <- xy[2,]
            attr(xy2, 'mask') <- mask
            Euclid <- as.numeric(dist(mask[matched,]))
            network <- as.numeric(networkdistance(xy1, xy2, geometry = mask))
            temp <- data.frame(from = matched[1],
                               to = matched[2],
                               Euclidean.d = Euclid,
                               network.d = network)
            if (is.null(results)) {
                results <- temp
            }
            else {
                results <- rbind(results,temp)
            }
            cat('Mask points       ', matched, '\n')
            cat('Euclidean distance', Euclid, ' m \n')
            cat('Network distance  ', network, ' m \n\n')
            flush.console()
        }
    }
    results <- list(paths = paths, results = results)
    invisible(results)
}

## buildTopo = function(lines) {
## Another approach - just kept for future reference
## Barry Rowlingson
## http://rstudio-pubs-static.s3.amazonaws.com/1572_7599552b60454033a0d5c5e6d2e34ffb.html
##     g = gIntersection(lines, lines)
##     edges = do.call(rbind, lapply(g@lines[[1]]@Lines, function(ls) {
##         as.vector(t(ls@coords))
##     }))
##     lengths = sqrt((edges[, 1] - edges[, 3])^2 + (edges[, 2] - edges[, 4])^2)
##
##     froms = paste(edges[, 1], edges[, 2])
##     tos = paste(edges[, 3], edges[, 4])
##
##     graph = graph.edgelist(cbind(froms, tos), directed = FALSE)
##     E(graph)$weight = lengths
##
##     xy = do.call(rbind, strsplit(V(graph)$name, " "))
##
##     V(graph)$x = as.numeric(xy[, 1])
##     V(graph)$y = as.numeric(xy[, 2])
##     return(graph)
## }
