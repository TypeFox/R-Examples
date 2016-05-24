


## FIXME: want a convenience function that behaves like heatmap()





## Goal: create a grob that could usefully represent a dendrogram


## can a dendrogram node have more than 2 children?
## long term FIXME: use better graph layout algorithms?




## returns a modified dendrogram object, with an extra attribute
## 'position=c(x, y)' for each node

## FIXME: should have and honor 'center=FALSE' argument

addPositions <-
    function(x, order)
{
    if (!is.null(attr(x, "position"))) return(x)
    else if (is.leaf(x))
    {
        attr(x, "position") <-
            list(x = which(x == order)[1],
                 y = attr(x, "height"))
        return(x)
    }
    else
    {
        for (i in seq_along(x))
        {
            x[[i]] <- addPositions(x[[i]], order)
        }
        attr(x, "position") <-
            list(x = mean(sapply(x, function(x) attr(x, "position")$x )),
                 y = attr(x, "height"))
        return(x)
    }
}

## returns a vector data.frame(x0, y0, x1, y1, ...), to be used in a
## call to segmentsGrob after being combined.  The possibility of
## attaching parameters exists, but is not (or barely) tested

edgeLocation <-
    function(pos.node, pos.child, type, ...)
{
    switch(type,
           rectangle = {
               data.frame(x0 = c(pos.node$x, pos.child$x),
                          y0 = c(pos.node$y, pos.node$y),
                          x1 = c(pos.child$x, pos.child$x),
                          y1 = c(pos.node$y, pos.child$y),
                          ..., stringsAsFactors = FALSE) ## 'col' can be strings
           },
           triangle = {
               data.frame(x0 = pos.node$x,  y0 = pos.node$y,
                          x1 = pos.child$x, y1 = pos.child$y,
                          ..., stringsAsFactors = FALSE) ## 'col' can be strings
           })
}


dendrogramGrob <- 
    function(x, ord = order.dendrogram(x),
             side = c("right", "top"),
             add = list(),
             size = 5,
             size.add = 1,
             type = c("rectangle", "triangle"),
             ...)
{
    ## Note: We use dendrapply() a couple of times.  The return value
    ## is unused (we are only interested in side-effects), but certain
    ## types of return values of FUN can make dendrapply() go into an
    ## infinite loop.  To be safe, we return original node.

    if (size <= 0) return(textGrob(label = NULL))
    type <- match.arg(type)
    native.height <- attr(x, "height")
    native.xscale <- c(1, length(ord)) + c(-1, 1) * lattice.getOption("axis.padding")$factor
    xpos <- addPositions(x, ord) ## version of x with positions

    ## how many non-leaf nodes are there?  For a binary tree, n-1,
    ## where n is the number of leaves (join any 2 ==> nodes++,
    ## leaves--), but we're more tolerant

    nnodes <- 0
    dendrapply(xpos,
               function(x) {
                   if (!is.leaf(x)) nnodes <<- nnodes + 1
                   x
               })
    xseg <- vector(mode = "list", length = nnodes)

    ## FIXME: add something similar to have nodes drawn as points
    i <- 0
    getSegments <- function(x, ...)
    {
        if (!is.leaf(x))
        {
            i <<- i + 1
            pos.node <- attr(x, "position")
            xseg[[i]] <<-
                do.call(rbind,
                        lapply(x,
                               function(child) {
                                   pos.child <- attr(child, "position")
                                   edgeLocation(pos.node, pos.child,
                                                type = type,
                                                ...)
                               }))
        }
        x
    }
    dendrapply(xpos, getSegments)
    all.segs <- do.call(rbind, xseg)
    ## number of additional indicators
    nadd <- length(add)
    ## nleaf <- length(ord)
    native.unit <- 1 / diff(native.xscale) # side of one square
    
    switch(side,
           right = {
               key.layout <-
                   grid.layout(nrow = 1, ncol = 1 + nadd,
                               heights = unit(1, "null"),
                               widths =
                               unit(c(rep(size.add, length = nadd), size),
                                    c(rep("lines", nadd), "lines")),
                               respect = FALSE)
               key.gf <- frameGrob(layout = key.layout)
               ## key.gf <- placeGrob(key.gf, rectGrob(gp = gpar(fill = "pink")))
               for (i in seq_len(nadd))
               {
                   addi <- add[[i]]
                   typei <- names(add)[i]
                   switch(typei,
                          rect = {
                              key.gf <-
                                  placeGrob(key.gf,
                                            rectGrob(y = (order(ord) - native.xscale[1]) * native.unit,
                                                     height = native.unit,
                                                     gp = do.call(gpar, addi)),
                                            row = 1, col = i)
                          })
               }
               key.gf <-
                   placeGrob(key.gf, 
                             with(all.segs,
                                  segmentsGrob((y0 / native.height),
                                               (x0 - native.xscale[1]) * native.unit,
                                               (y1 / native.height),
                                               (x1 - native.xscale[1]) * native.unit)),
                             row = 1, col = 1 + nadd)
               key.gf
           },
           top = {
               key.layout <-
                   grid.layout(nrow = 1 + nadd, ncol = 1,
                               widths = unit(1, "null"),
                               heights =
                               unit(c(size, rep(size.add, length = nadd)),
                                    c("lines", rep("lines", nadd))),
                               respect = FALSE)

               key.gf <- frameGrob(layout = key.layout)
               ## key.gf <- placeGrob(key.gf, rectGrob(gp = gpar(fill = "pink")))

               for (i in seq_len(nadd))
               {
                   addi <- add[[i]]
                   typei <- names(add)[i]
                   switch(typei,
                          rect = {
                              key.gf <-
                                  placeGrob(key.gf,
                                            rectGrob(x = (order(ord) - native.xscale[1]) * native.unit,
                                                     width = native.unit,
                                                     gp = do.call(gpar, addi)),
                                            row = 1 + i, col = 1)
                          })
               }
               key.gf <-
                   placeGrob(key.gf, 
                             with(all.segs,
                                  segmentsGrob((x0 - native.xscale[1]) * native.unit,
                                               (y0 / native.height),
                                               (x1 - native.xscale[1]) * native.unit,
                                               (y1 / native.height))),
                             row = 1, col = 1)
               key.gf
           })
}



