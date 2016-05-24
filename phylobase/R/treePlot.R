##' Phylogeny plotting
##'
##' Plot \code{phylo4} or \code{phylo4d} objects, including associated data.
##'
##'
##' @name treePlot-methods
##' @aliases treePlot plot,ANY,ANY-method plot,pdata,missing-method
##' plot,phylo4,missing-method treePlot-method treePlot,phylo4,phylo4d-method
##' @docType methods
##' @details Currently, \code{treePlot} can only plot numeric values
##' for tree-associated data. The dataset will be subset to only
##' include columns of class \code{numeric}, \code{integer} or
##' \code{double}. If a \code{phylo4d} object is passed to the
##' function and it contains no data, or if the data is in a format
##' that cannot be plotted, the function will produce a warning. You
##' can avoid this by using the argument \code{plot.data=FALSE}.
##' @param phy A \code{phylo4} or \code{phylo4d} object
##' @param x A \code{phylo4} or \code{phylo4d} object
##' @param y (only here for compatibility)
##' @param type A character string indicating the shape of plotted tree
##' @param show.tip.label Logical, indicating whether tip labels should be shown
##' @param show.node.label Logical, indicating whether node labels should be
##' shown
##' @param tip.order If NULL the tree is plotted with tips in preorder, if "rev"
##' this is reversed. Otherwise, it is a character vector of tip labels,
##' indicating their order along the y axis (from top to bottom). Or, a numeric
##' vector of tip node IDs indicating the order.
##' @param plot.data Logical indicating whether \code{phylo4d} data should be
##' plotted
##' @param rot Numeric indicating the rotation of the plot in degrees
##' @param tip.plot.fun A function used to generate plot at the each tip of the
##' phylogenetic trees
##' @param edge.color A vector of colors in the order of \code{edges(phy)}
##' @param node.color A vector of colors indicating the colors of the node
##' labels
##' @param tip.color A vector of colors indicating the colors of the tip labels
##' @param edge.width A vector in the order of \code{edges(phy)} indicating the
##' widths of edge lines
##' @param newpage Logical indicating whether the page should be cleared before
##' plotting
##' @param plot.at.tip should the data plots be at the tip? (logical)
##' @param margins number of lines around the plot (similar to \code{par(mar)}).
##' @param \dots additional arguments
##' @return No return value, function invoked for plotting side effect
##' @section Methods: \describe{ \item{phy = "phylo4"}{plots a tree of class
##' \linkS4class{phylo4}} \item{phy = "phylo4d"}{plots a tree with one or more
##' quantitative traits contained in a \linkS4class{phylo4d} object.} }
##' @author Peter Cowan \email{pdc@@berkeley.edu}, Francois Michonneau
##' @seealso \code{\link{phylobubbles}}
##' @keywords methods
##' @export
##' @examples
##'
##' ## example of plotting two grid plots on the same page
##' library(grid)
##' data(geospiza)
##' geotree <- extractTree(geospiza)
##' grid.newpage()
##' pushViewport(viewport(layout=grid.layout(nrow=1, ncol=2), name="base"))
##'   pushViewport(viewport(layout.pos.col=1, name="plot1"))
##'     treePlot(geotree, newpage=FALSE)
##'   popViewport()
##'
##'   pushViewport(viewport(layout.pos.col=2, name="plot2"))
##'     treePlot(geotree, newpage=FALSE, rot=180)
##' popViewport(2)
`treePlot` <- function(phy,
                     type = c('phylogram', 'cladogram', 'fan'),
                     show.tip.label = TRUE,
                     show.node.label = FALSE,
                     tip.order = NULL,
                     plot.data = is(phy, 'phylo4d'),
                     rot = 0,
                     tip.plot.fun = 'bubbles',
                     plot.at.tip = TRUE,
                     edge.color = 'black',
                     node.color = 'black', # TODO what do with node.color parameter
                     tip.color  = 'black',
                     edge.width = 1, # TODO line-type modification hack
                     newpage = TRUE,
                     margins = c(1.1, 1.1, 1.1, 1.1), # number of lines, same as par(mar)
                     ...
            )
{
    ## TODO three dimensional histogram as example, compute values on full dataset
    ## then generate phylo4d object with summary data and plot

    ## TODO factors not handled in data plots
    ## TODO add symbols at the nodes, allow coloirng and sizing downViewport approach?
    ## TODO cladogram methods incorrect
    ## because we may reoder the tip, we need to update the phy objec

    if (!inherits(phy, 'phylo4')) stop('treePlot requires a phylo4 or phylo4d object')
    if (!isRooted(phy)) stop("treePlot function requires a rooted tree.")
    if (plot.data) {
        if (!hasTipData(phy)) {
            warning("tree has no tip data to plot")
            plot.data <- FALSE
        }
        else {
            ## if new plotting functions are developped that allow users to plot other type of data
            ## this needs to be removed/adjusted
            ##  other checks are being made in phylobubbles()
            if (!any(sapply(tdata(phy, "tip"), function(x) class(x) %in% c("numeric", "double", "integer")))) {
                warning("only numeric data can be plotted at this time")
                plot.data <- FALSE
            }
        }
    }
    if (hasRetic(phy))
        stop("treePlot requires non-reticulated trees.")

    if(newpage) grid::grid.newpage()
    type   <- match.arg(type)
    Nedges <- nEdges(phy)
    Ntips  <- nTips(phy)
    if(!is.null(tip.order) && length(tip.order) > 1) { ## if length of tip.order is more than 1 it can't be "rev"
        if(length(tip.order) != Ntips) {stop('tip.order must be the same length as nTips(phy)')}
        if(is.numeric(tip.order)) {
            tip.order <- tip.order
        }
        else {
            if(is.character(tip.order)) {
                tip.order <- as.numeric(names(tipLabels(phy))[match(tip.order, tipLabels(phy))])
            }
        }
        tip.order <- rev(tip.order)
    }
    ## TODO remove the false cladogram option?
    if(!hasEdgeLength(phy) || type == 'cladogram') {
        edgeLength(phy) <- rep(1, Nedges)
    }
    xxyy   <- phyloXXYY(phy, tip.order)
    if(type == 'cladogram') {
        xxyy$xx[edges(xxyy$phy)[, 2] <= Ntips] <- 1
    }

    ## plotViewport is a convience function that provides margins in lines
    grid::pushViewport(grid::plotViewport(margins=margins))

    if(!plot.data) {
        plotOneTree(xxyy, type, show.tip.label, show.node.label, edge.color,
                                node.color, tip.color, edge.width, rot)
    } else {
        if(!is.function(tip.plot.fun)) {
            if(tip.plot.fun == "bubbles") {
                phylobubbles(
                    type = type,
                    show.node.label = show.node.label,
                    rot = 0,
                    edge.color = edge.color,
                    node.color = node.color, # TODO what do with node.color parameter
                    tip.color  = tip.color,
                    edge.width = edge.width, # TODO line-type modification hack
                    newpage = TRUE,
                    ..., XXYY = xxyy
                )
            } else {
                stop(paste(tip.plot.fun, 'is neither a function or a recognized plot type'))
            }
        } else { ## from -- if(tip.plot.fun == "bubbles")
            ## plot.at.tip <- TRUE
            if (plot.at.tip) {
                tip.data.plot(
                    xxyy = xxyy,
                    type = type,
                    show.tip.label = show.tip.label,
                    show.node.label = show.node.label,
                    rot = 0,
                    tip.plot.fun = tip.plot.fun,
                    edge.color = edge.color,
                    node.color = node.color, # TODO what do with node.color parameter
                    tip.color  = tip.color,
                    edge.width = edge.width, # TODO line-type modification hack
                    newpage = TRUE,
                    ...
                )
                return(invisible())
            } ## if (plot.at.tip)
        } ## else
    } ## else
    grid::upViewport() # margins
}



##' Plot a phylo4 object
##'
##' Plots the phylogenetic tree contained in a \code{phylo4} or \code{phylo4d}
##' object.
##'
##'
##' @param xxyy A list created by the \code{\link{phyloXXYY}} function
##' @param type A character string indicating the shape of plotted tree
##' @param show.tip.label Logical, indicating whether tip labels should be shown
##' @param show.node.label Logical, indicating whether node labels should be
##' shown
##' @param edge.color A vector of colors in the order of \code{edges(phy)}
##' @param node.color A vector of colors indicating the colors of the node
##' labels
##' @param tip.color A vector of colors indicating the colors of the tip labels
##' @param edge.width A vector in the order of \code{edges(phy)} indicating the
##' widths of edge lines
##' @param rot Numeric indicating the rotation of the plot in degrees
##' @return Returns no values, function invoked for the plotting side effect.
##' @author Peter Cowan \email{pdc@@berkeley.edu}
##' @seealso \code{treePlot}, \code{\link{phyloXXYY}}
##' @export
##' @keywords methods
##' @examples
##' library(grid)
##' data(geospiza)
##' grid.newpage()
##' xxyy <- phyloXXYY(geospiza)
##' plotOneTree(xxyy, type = 'phylogram',
##'   show.tip.label = TRUE, show.node.label = TRUE,
##'   edge.color = 'black', node.color = 'orange', tip.color = 'blue',
##'   edge.width = 1, rot = 0
##' )
##'
##' grid.newpage()
##' pushViewport(viewport(w = 0.8, h = 0.8))
##' plotOneTree(xxyy, type = 'phylogram',
##'   show.tip.label = TRUE, show.node.label = TRUE,
##'   edge.color = 'black', node.color = 'orange', tip.color = 'blue',
##'   edge.width = 1, rot = 0
##' )
##' popViewport()
##'

plotOneTree <- function(xxyy, type, show.tip.label, show.node.label, edge.color,
                        node.color, tip.color, edge.width, rot)
{
    # TODO switch to phylobase abstractions
    phy    <- xxyy$phy
    Nedges <- nEdges(phy)
    Nnodes <- nNodes(phy)
    Ntips  <- nTips(phy)
    pedges <- edges(phy)
    tindex <- pedges[pedges[, 2] <= Ntips, 2]
    eindex <- xxyy$eorder
    segs   <- xxyy$segs

    ## TODO check that colors are valid?
    if(length(edge.color) != Nedges) {
        edge.color <- rep(edge.color, length.out = Nedges)
    }
    edge.color <- edge.color[eindex]

    if(length(edge.width) != Nedges) {
        edge.width <- rep(edge.width, length.out = Nedges)
    }
    edge.width <- edge.width[eindex]

    ## TODO check that colors are valid?
    if(length(node.color) != Nnodes) {
        node.color <- rep(node.color, length.out = Nnodes)
    }

    if(show.tip.label) {
        ## calculate several lab dimesisions
        ## labw    -- a vector of string widths
        ## adjlabw -- the max width for adjusting the size of viewports
        ## laboff  -- a vector of half string widths for
        ## offsetting center justified labels, handy for vp rotation
        labw    <- grid::stringWidth(tipLabels(phy))
        adjlabw <- max(labw) + grid::unit(0.1, 'inches')
        laboff  <- labw * 0.5 + grid::unit(0.1, 'inches')
        ## print(foo <<- laboff)
        treelayout <- grid.layout(nrow = 1, ncol = 2,
            widths = grid::unit.c(grid::unit(1, 'null', NULL), grid::convertUnit(adjlabw, 'inches'))
            )
        tindex <- pedges[pedges[, 2] <= Ntips, 2]
        if(length(tip.color) != Ntips) {
            tip.color <- rep(tip.color, length.out = Ntips)
        }
        # keep labels horizontal unless plot is upwards or downwards
        lrot <- ifelse(rot %% 360 %in% c(90, 270), 0, -rot)
    } else {
        treelayout <- grid::grid.layout(nrow = 1, ncol = 1)
    }
    # grid.show.layout(treelayout)
    grid::pushViewport(grid::viewport(
        x = 0.5, y = 0.5,
        width = 1, height = 1,
        layout = treelayout, angle = rot, name = 'treelayout'))
    grid::pushViewport(grid::viewport(
        layout.pos.col = 1,
        name = 'tree'))
    if (type == "fan") {
        dseg <- grid::grid.segments( # draws diag lines
            x0 = segs$v0x, y0 = segs$v0y,
            x1 = segs$h1x, y1 = segs$h1y,
            name = "diag", gp = grid::gpar(col = edge.color, lwd = edge.width))
    } else {
        vseg <- grid::grid.segments( # draws vertical lines
            x0 = segs$v0x, y0 = segs$v0y,
            x1 = segs$v1x, y1 = segs$v1y,
            name = "vert", gp = grid::gpar(col = edge.color, lwd = edge.width))
        hseg <- grid::grid.segments( # draws horizontal lines
            x0 = segs$h0x, y0 = segs$h0y,
            x1 = segs$h1x, y1 = segs$h1y,
            name = "horz", gp = grid::gpar(col = edge.color, lwd = edge.width))
    }
    grid::upViewport() # tree
    if(show.tip.label) {
        grid::pushViewport(grid::viewport(layout.pos.col = 1,
            name = 'tiplabelvp'))
        labtext <- grid::grid.text(
            tipLabels(phy)[tindex],
            x = grid::unit(xxyy$xx[pedges[, 2] %in% tindex], "native") + laboff[tindex],
            y = xxyy$yy[pedges[, 2] %in% tindex], rot = lrot,
            default.units = 'native', name = 'tiplabels',
            just = 'center', gp = grid::gpar(col = tip.color[tindex])
        )
        grid::upViewport() #tiplabelvp
    }
    # TODO probably want to be able to adjust the location of these guys
    if(show.node.label) {
        grid::pushViewport(grid::viewport(layout = treelayout, layout.pos.col = 1, name = 'nodelabelvp'))
        theLabels <- nodeLabels(phy)[match(pedges[pedges[, 2] > Ntips, 2], names(nodeLabels(phy)))]
        ## don't plot NAs
        theLabels[is.na(theLabels)] <- ""
        labtext <- grid::grid.text(
            theLabels,
            x = c(xxyy$xx[pedges[, 2] > Ntips]),
            y = c(xxyy$yy[pedges[, 2] > Ntips]),
            default.units = 'npc', name = 'nodelabels', rot = -rot,
            just = 'center', gp = grid::gpar(col = node.color)
        )
        grid::upViewport() #nodelabelvp
    }
    grid::upViewport() # treelayout
    # grobTree(vseg, hseg, labtext)
}



##' Calculate node x and y coordinates
##'
##' Calculates the node x and y locations for plotting a phylogenetic tree.
##'
##' The y coordinates of the tips are evenly spaced from 0 to 1 in pruningwise
##' order.  Ancestor y nodes are given the mean value of immediate descendants.
##' The root is given the x coordinate 0 and descendant nodes are placed
##' according to the cumulative branch length from the root, with a maximum x
##' value of 1.
##'
##' @param phy A \code{phylo4} or \code{phylo4d} object.
##' @param tip.order A character vector of tip labels, indicating their order
##' along the y axis (from top to bottom). Or, a numeric vector of tip node IDs
##' indicating the order.
##' @return \item{yy}{Internal node and tip y coordinates} \item{xx}{Internal
##' node and tip x coordinates} \item{phy}{A \code{phylo4} or \code{phylo4d}
##' object} \item{segs}{A list of \code{h0x, h1x, v0x, v1x} and \code{h0y, h1y,
##' v0y, v1y} describing the start and end points for the plot line segments}
##' \item{torder}{The tip order provided as \code{tip.order} or if NULL the
##' preoder tip order} \item{eorder}{The an index of the reordered edges
##' compared to the result of \code{edges(phy)}}
##' @author Peter Cowan \email{pdc@@berkeley.edu}
##' @seealso \code{treePlot}, \code{\link{plotOneTree}}
##' @export
##' @keywords methods
##' @examples
##'
##'
##' data(geospiza)
##' coor <- phyloXXYY(geospiza)
##' plot(coor$xx, coor$yy, pch = 20)
##'
##'
phyloXXYY <- function(phy, tip.order=NULL)
{
    phy.orig <- phy
    ## initalize the output
    phy    <- reorder(phy, 'preorder')
    pedges <- edges(phy)
    eindex <- match(pedges[,2], edges(phy.orig)[,2])
    Nedges <- nrow(pedges) ## TODO switch to the accessor once stablized
    Ntips  <- nTips(phy)
    tips <- pedges[, 2] <= Ntips
    xx <- numeric(Nedges)
    yy <- numeric(Nedges)

    treelen <- rep(NA, nEdges(phy))
    segs <- list(v0x = treelen, v0y = treelen, v1x = treelen, v1y = treelen,
                 h0x = treelen, h0y = treelen, h1x = treelen, h1y = treelen)

    ## Set root x value to zero and calculate x positions
    xx[1] <- 0
    segs$v0x[1] <- segs$v1x[1] <- segs$h0x[1] <- 0
    edge1   <- as.integer(pedges[,1])
    edge2   <- as.integer(pedges[,2])
    edgeLen <- edgeLength(phy)
    edgeLen[is.na(edgeLen)] <- 0
    edgeLen <- as.numeric(edgeLen)
    nedges  <- as.integer(nEdges(phy))
    segsv0x <- as.numeric(rep.int(0, Nedges))
    xPos <- .C("phyloxx", edge1, edge2,
            edgeLen, nedges, xx, segsv0x)
    xx <- xPos[[5]]
    segs$v0x <- xPos[[6]]

    ## Set y positions for terminal nodes and calculate remaining y positions
    if(!is.null(tip.order)) {
        if(length(tip.order) == 1 &&  tip.order == "rev") {
            yy[tips] <- seq(1, 0, length = Ntips)
            tip.order <- rev(edge2[edge2 <= Ntips])
        } else {
            yy[tips][match(tip.order, edge2[tips])] <- seq(0, 1, length = Ntips)
        }
    } else {
        yy[tips] <- seq(0, 1, length = Ntips)
        tip.order <- edge2[edge2 <= Ntips]
    }
    segs$h0y[tips] <- segs$h1y[tips] <- yy[tips]
    segs$v1y[tips] <- segs$v0y[tips] <- yy[tips]
    phyloyy <- function() {
        for(i in rev((Ntips + 1):nEdges(phy))) {
            dex <- pedges[, 1] == i
            cur <- pedges[, 2] == i
            yy[cur] <- segs$v0y[dex] <- mean(yy[dex])
        }
        return(list(segs=segs, yy=yy))
    }

    yPos <- phyloyy()
    segs <- yPos$segs
    yy   <- yPos$yy

    ## edgeLen[is.na(edgeLen)] <- 0
    ## edgeLen <- as.numeric(edgeLen)
    ## ntips   <- as.integer(nTips(phy))
    ## yy      <- as.numeric(yy)
    ## segsv0y <- as.numeric(yy)

    ## yPos <- .C("phyloyy", edge1, edge2,
    ##         ntips, nedges, yy, segsv0y)

    segs$h0y <- segs$h1y <- segs$v1y <- yy

    ## scale the x values so they range from 0 to 1
    Xmax <- max(xx)
    segs$v0x <- segs$v0x / Xmax
    xx <- xx / Xmax

    segs$h1x <- xx
    segs$v1x <- segs$h0x <- segs$v0x

    # TODO return an index vector instead of a second phy object
    list(xx = xx, yy = yy, phy = phy, segs = segs, torder=tip.order, eorder=eindex)
}

.bubLegendGrob <- function(tipdata, tipdataS) {
    grid::grob(tipdata=tipdata, tipdataS=tipdataS, cl='bubLegend')
}

drawDetails.bubLegend <- function(x, ...) {
    ## number of bubbles in legend
    leglen  <- 4
    ## the raw tip data
    tipdata <- x$tipdata
    ## the tip data as scaled for bubble plot
    ts <- x$tipdataS
    ## return to the bubble plot viewport to get properly scaled values
    ## this relies on having well named unique viewports
    grid::seekViewport("bubble_plots")
        ## retreive the min and max non-zero bubbles as numerics not units
        bubrange <- grid::convertUnit(
                    grid::unit(c(min(ts[ts != 0], na.rm=TRUE), max(ts[ts != 0], na.rm=TRUE)), "native"),
                    "mm", valueOnly=TRUE)
    grid::seekViewport("bubblelegend")
    ## grid.rect()
    ## Generate the sequence of legend bubble sizes and convert to grid mm units
    legcirS  <- grid::unit(seq(bubrange[1], bubrange[2], length.out=leglen), "mm")
    ## get the corresponding sequence of actual data values
    legcir   <- seq(min(tipdata[tipdata != 0], na.rm=TRUE),
                    max(tipdata[tipdata != 0], na.rm=TRUE), length.out=leglen)
    ccol     <- ifelse(legcir < 0, 'black', 'white')

    leftedge <- abs(grid::convertUnit(legcirS[1], 'npc', valueOnly=TRUE)) + 0.1
    xloc     <- seq(leftedge, 0.5, length.out=leglen)
    textsp   <- grid::convertUnit(max(abs(legcirS)), axisFrom="y", axisTo="y", 'npc', valueOnly=TRUE)
    strsp    <- grid::convertUnit(unit(1, "strheight", "TTT"), axisFrom="y", 'npc', valueOnly=TRUE)
    grid::grid.circle(x=xloc, y=0.9 - textsp - strsp, r=legcirS, gp = grid::gpar(fill=ccol), default.units = 'npc')
    grid::grid.text(as.character(signif(legcir, digits = 2)),
                x=xloc, y=0.75 - 2 * textsp - strsp,
                gp=grid::gpar(cex=0.75),
                default.units='npc'
    )
}



##' Bubble plots for phylo4d objects
##'
##' Plots either circles or squares corresponding to the magnitude of each cell
##' of a \code{phylo4d} object.
##'
##'
##' @param type the type of plot
##' @param place.tip.label A string indicating whether labels should be plotted
##' to the right or to the left of the bubble plot
##' @param show.node.label A logical indicating whether internal node labels
##' should be plotted
##' @param rot The number of degrees that the plot should be rotated
##' @param edge.color A vector of colors for the tree edge segments
##' @param node.color A vector of colors for the coloring the nodes
##' @param tip.color A vector of colors for the coloring the tip labels
##' @param edge.width A vector of line widths for the tree edges
##' @param newpage Logical to control whether the device is cleared before
##' plotting, useful for adding plot inside other plots
##' @param \dots Additional parameters passed to the bubble plotting functions
##' @param XXYY The out put from the phyloXXYY function
##' @param square Logical indicating whether the plot 'bubbles' should be
##' squares
##' @param grid A logical indicating whether a grey grid should be plotted
##' behind the bubbles
##' @author Peter Cowan \email{pdc@@berkeley.edu}
##' @export
##' @seealso \code{\link{phyloXXYY}}, \code{treePlot}
##' @keywords methods
##' @examples
##'
##' ##---- Should be DIRECTLY executable !! ----
##' ##-- ==>  Define data, use random,
##' ##--	or do  help(data=index)  for the standard data sets.
##'
phylobubbles <- function(type = type,
                        place.tip.label = "right",
                        show.node.label = show.node.label,
                        rot = 0,
                        edge.color = edge.color,
                        node.color = node.color, # TODO what do with node.color parameter
                        tip.color  = tip.color,
                        edge.width = edge.width, # TODO line-type modification hack
                        newpage = TRUE,
                        ...,
                        XXYY, square = FALSE, grid = TRUE)
{
    ## TODO add legend command
    ## tys    -- tip y coordinates
    ## nVars  -- number of traits/characters
    ## maxr   -- maximum circle radius, based on nVars or nTips
    ## torder -- the order of tips in the reordered edge matrix
    if(rot != 0) {stop("Rotation of bubble plots not yet implemented")}
    lab.right <- ifelse(place.tip.label %in% c("right", "both"), TRUE, FALSE)
    lab.left  <- ifelse(place.tip.label %in% c("left", "both"), TRUE, FALSE)

    phy       <- XXYY$phy
    tip.order <- XXYY$torder
    tipdata   <- tdata(phy, type = "tip")[tip.order,, drop=FALSE]
    tipClass <- sapply(tipdata, function(x) class(x) %in% c("double", "integer", "numeric"))
    tipdata <- tipdata[, tipClass, drop=FALSE]
    tmin      <- min(tipdata, na.rm = TRUE)
    tmax      <- max(tipdata, na.rm = TRUE)
    pedges    <- edges(phy)

    nVars     <- ncol(tipdata) # number of bubble columns

    dlabwdth <- max(grid::stringWidth(colnames(tipdata))) * 1.2
    if(grid::convertWidth(dlabwdth, 'cm', valueOnly=TRUE) < 2) {dlabwdth <- grid::unit(2, 'cm')}
    phyplotlayout <- grid::grid.layout(nrow = 2, ncol = 2,
        heights = grid::unit.c(grid::unit(1, 'null'), dlabwdth),
        widths = grid::unit(c(1, 1), c('null', 'null'), list(NULL, NULL)))
    grid::pushViewport(viewport(layout = phyplotlayout, name = 'phyplotlayout'))
    grid::pushViewport(viewport(layout.pos.row = 1:2, layout.pos.col = 2,
                height = grid::unit(1, 'npc') +
                                grid::convertUnit(dlabwdth, 'npc'),
                name = 'bubbleplots', default.units = 'native'))

    # tip y coordinates
    tys <- XXYY$yy[pedges[, 2] <= nTips(phy)]
    tys <- tys[match(names(tipLabels(phy))[tip.order], XXYY$torder)]

    maxr <- ifelse(ncol(tipdata) > nTips(phy), 1 / ncol(tipdata), 1 / nTips(phy))
    tipdataS <- apply(tipdata, 2,
                      function(x) (maxr * x) / max(abs(x), na.rm = TRUE))
    if(nVars == 1) {
        xpos <- 0.5
    } else {
        xpos <- seq(0 + maxr + 0.02, 1 - maxr - 0.02, length.out = nVars)
    }

    ## rep coordinates for filling a matrix columnwise
    xrep <- rep(xpos, each = length(tys))
    yrep <- rep(tys, nVars)
    ## color bubbles

    ccol <- ifelse(tipdata < 0, 'black', 'white')

    ## generate matrices of every x and y by filling the repd value columnwise
    ## then subset for datapoints that are NA
    naxs <- matrix(xrep, ncol = nVars)
    nays <- matrix(yrep, ncol = nVars)
    dnas <- is.na(tipdataS)
    naxs <- naxs[dnas]
    nays <- nays[dnas]
    ## set the NA points to zero so that grid.circle doesn't crash
    tipdataS[is.na(tipdataS)] <- 0 + 0.001  # workaround negative circles on PDF

    ## get label widths
    if(lab.right) {
        tiplabwidth  <- max(grid::stringWidth(tipLabels(phy)))
    } else {tiplabwidth <- grid::unit(0, 'null', NULL)}

    ## 2x2 layout -- room at the bottom for data labels, and legend
    bublayout <- grid::grid.layout(nrow = 2, ncol = 2,
        widths  = grid::unit.c(grid::unit(1, 'null', NULL), tiplabwidth),
        heights = grid::unit.c(grid::unit(1, 'null', NULL), dlabwdth))
    grid::pushViewport(viewport(
        x = 0.5, y = 0.5,
        width = 0.95, height = 1,
        layout = bublayout, name = 'bublayout'
    ))
    grid::pushViewport(viewport(
        name = 'bubble_plots',
        layout = bublayout,
        layout.pos.col = 1,
        layout.pos.row = 1
    ))
    if(grid) {
        ## draw light grey grid behind bubbles
        grid::grid.segments(x0 = 0,   x1 = 1,
                      y0 = tys, y1 = tys, gp = grid::gpar(col = 'grey'))
        grid::grid.segments(x0 = xpos, x1 = xpos,
                      y0 = 0,    y1 = 1, gp = grid::gpar(col = 'grey'))
    }
    if (length(naxs) > 0) {
        ## if ther are missing values plot Xs
        grid::grid.points(naxs, nays, pch = 4)
    }

    if(square) {
        ## alternative to circles
        ## to keep the squares square, yet resize nicely use the square npc
        sqedge <- grid::unit(unlist(tipdataS), 'snpc')
        grid::grid.rect(x = xrep, y = yrep,
            width = sqedge,
            height = sqedge,
            gp=grid::gpar(fill = ccol))
    } else {
        ## plot bubbles
        grid::grid.circle(xrep, yrep, r = unlist(tipdataS), gp = grid::gpar(fill = ccol))
    }
    grid::upViewport()

    ## push view ports for tip and data labels fixed locations
    if(lab.right) {
        grid::pushViewport(viewport(
            name = 'bubble_tip_labels',
            layout = bublayout,
            layout.pos.col = 2,
            layout.pos.row = 1
        ))
        tt <- tipLabels(phy)[tip.order] # phy@tip.label
        grid::grid.text(tt, 0.1, tys, just = 'left')
        grid::upViewport()
    }
    grid::pushViewport(viewport(
        name = 'bubble_data_labels',
        layout = bublayout,
        layout.pos.col = 1,
        layout.pos.row = 2
    ))
    ## ideas, for nicer sizing of the data labels
    ## data.label.space <- convertX(unit(1, 'npc'), "points", valueOnly = TRUE)
    ## data.label.fontsize <- data.label.space / ncol(tipdata)
    ## , gp=gpar(fontsize=data.label.fontsize))
    ## offset the data labels from the bottom bubble
    datalaboffset <- grid::convertUnit(grid::unit(15, "mm"), 'npc', valueOnly=TRUE)
    grid::grid.text(colnames(tipdata), xpos, 1-datalaboffset, rot = 90, just = 'right')

    grid::upViewport(3)
    grid::pushViewport(viewport(layout.pos.row=2, layout.pos.col=1,
                name='bubblelegend'))
    yyy <- .bubLegendGrob(tipdata, tipdataS)
    grid::grid.draw(yyy)
    grid::upViewport()

    grid::pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1,
                name = 'tree'))
        plotOneTree(XXYY, type, show.tip.label=lab.left, show.node.label, edge.color,
                                node.color, tip.color, edge.width, rot)
    grid::upViewport(2)

    # to make a nice legend, return the biggest smallest and a scaling factor
    # translate the scale of the current vp to a fixed value
    ## ensure the min is not a zero (or NA) that's replaced by a zero
    ## print(convertUnit(bubscale, 'inches', valueOnly = TRUE))
    ## return(list(max = max(tipdata, na.rm = TRUE),
    ##             min = min(tipdata[tipdata != 0], na.rm = TRUE),
    ##             has.na = length(naxs) > 0,
    ##             bubscale = bubscale))
}



##' Plotting trees and associated data
##'
##' Plotting phylogenetic trees and associated data
##'
##'
##' @param xxyy A list created by the \code{\link{phyloXXYY}} function
##' @param type A character string indicating the shape of plotted tree
##' @param show.tip.label Logical, indicating whether tip labels should be shown
##' @param show.node.label Logical, indicating whether node labels should be
##' shown
##' @param rot Numeric indicating the rotation of the plot in degrees
##' @param tip.plot.fun A function used to plot the data elements of a
##' \code{phylo4d} object
##' @param edge.color A vector of colors in the order of \code{edges(phy)}
##' @param node.color A vector of colors indicating the colors of the node
##' labels
##' @param tip.color A vector of colors indicating the colors of the tip labels
##' @param edge.width A vector in the order of \code{edges(phy)} indicating the
##' widths of edge lines
##' @param \dots Additional parameters passed to \code{tip.plot.fun}
##' @return creates a plot on the current graphics device.
##' @author Peter Cowan
##' @export
##' @keywords methods
tip.data.plot <- function(
                     xxyy,
                     type = c('phylogram', 'cladogram', 'fan'),
                     show.tip.label = TRUE,
                     show.node.label = FALSE,
                     rot = 0,
                     tip.plot.fun = grid.points,
                     edge.color = 'black',
                     node.color = 'black', # TODO what do with node.color parameter
                     tip.color  = 'black',
                     edge.width = 1, # TODO line-type modification hack
                     ...)
{
    phy    <- xxyy$phy
    tip.order <- xxyy$torder
    pedges <- edges(phy)
    Ntips  <- nTips(phy)
    datalayout <- grid::grid.layout(ncol = 2, widths = grid::unit(c(1, 1/Ntips), c('null', 'null')))
    # TODO this is done multiple times,
    grid::pushViewport(viewport(layout = datalayout, angle = rot,
                        name = 'datalayout'))
    grid::pushViewport(viewport(
        yscale = c(-0.5 / Ntips, 1 + 0.5 / Ntips),
        xscale = c(0, 1 + 1 / Ntips),
        layout.pos.col = 1,
        name = 'data_plots'))
    ## TODO should plots float at tips, or only along edge?
    hc <- grid::convertY(grid::unit(1 / Ntips, 'snpc'), 'npc')
    for(i in 1:Ntips) {
        grid::pushViewport(viewport(
            y = xxyy$yy[pedges[, 2] == i],
            x = 1 + 1 / (2 * Ntips), # xxyy$xx[phy@edge[, 2] == i],
            height = hc,
            width = hc,
            # default.units = 'native',
            name = paste('data_plot', i),
            just = "center",
            angle = -rot
            ))
            #grid.rect()
            tvals <- tdata(phy, type = 'tip')[nodeId(phy,'tip'), , drop=FALSE]
            vals = t(tvals[i, ])
            if (!all(is.na(vals))) tip.plot.fun(vals, ...)
        grid::upViewport() # loop viewports
    }
    plotOneTree(xxyy, type, show.tip.label, show.node.label, edge.color,
                            node.color, tip.color, edge.width, rot)
    grid::upViewport(2) ## data_plot & datalayout
}

# phyloStripchart <- function()

##' @rdname treePlot-methods
##' @aliases plot
##' @exportMethod plot
setGeneric('plot')

##' @rdname treePlot-methods
##' @aliases plot,phylo4-method
setMethod('plot', signature(x='phylo4', y='missing'), function(x, y, ...) {
    treePlot(x, ...)
})
