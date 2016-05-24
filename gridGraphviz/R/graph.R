
### Draw an Ragraph laid out graph using grid

## Code for "node" grobs

node <- function(label, x=.5, y=.5,
                 shape="plain",
                 height=NULL, lwidth=NULL, rwidth=NULL,
                 name,
                 color="black",
                 fillcolor="transparent",
                 fontcolor="black",
                 fontsize=10) {    
    # Shrink text if necessary
    cex <- 1
    if (!(is.null(height) || is.null(lwidth) || is.null(rwidth))) {
        tw <- convertWidth(grobWidth(textGrob(label,
                                              gp=gpar(fontsize=fontsize)))*1.4,
                           "inches", valueOnly=TRUE)
        nw <- convertWidth(lwidth + rwidth, "inches",
                           valueOnly=TRUE)
        if (tw > nw)
            cex <- nw/tw
        th <- convertHeight(grobHeight(textGrob(label,
                                                gp=gpar(fontsize=fontsize)))*1.4,
                            "inches", valueOnly=TRUE)
        nh <- convertHeight(height, "inches", valueOnly=TRUE)
        if (th > nh && nh/th < cex)
            cex <- nh/th
    }
    
    # make label grob
    lab <- makeLabelGrob(label, x, y, fontcolor, fontsize, cex, name=name)
    
    if (is.null(height)) {
        height <- grobHeight(lab)
    }
    if (is.null(lwidth) || is.null(rwidth)) {
        lwidth <- 0.5*grobWidth(lab)
        rwidth <- 0.5*grobWidth(lab)
    }

    # make box grob
    box <- makeBoxGrob(shape=shape, name=name, x=x, y=y, height=height,
                       lwidth=lwidth, rwidth=rwidth, color=color,
                       fillcolor=fillcolor)
    
    gTree(children=gList(box, lab),
          name=gsub("\n", "", name), cl="node")
}

# Dimensions of node come from "box" child (which depends on the shape)
grobX.node <- function(x, theta) {
    grobX(getGrob(x, "box"), theta)
}

grobY.node <- function(x, theta) {
    grobY(getGrob(x, "box"), theta)
}

drawNode <- function(node) {
    name <- name(node)
    shape <- shape(node)
    if (is.na(shape)) {
        shape <- "plain"
    }
    height <- getNodeHeight(node)
    lwidth <- getNodeLW(node)
    rwidth <- getNodeRW(node)
    col <- color(node)
    if (is.na(col) || col == "") {
        col <- "black"
    }
    fill <- fillcolor(node)
    if (is.na(fill) || fill == "") {
        fill <- "transparent"
    }
    fontcol <- labelColor(txtLabel(node))
    if (is.na(fontcol)) {
        fontcol <- "black"
    }
    fontsize <- labelFontsize(txtLabel(node))
    label <- labelText(txtLabel(node))
    xy <- getNodeXY(node)
    name <- name(node)
    grid.draw(node(label=label,
                   x=unit(xy$x, "native"),
                   y=unit(xy$y, "native"),
                   shape=shape,
                   height=unit(height, "native"),
                   lwidth=unit(lwidth, "native"),
                   rwidth=unit(rwidth, "native"),
                   name=name, color=col, fillcolor=fill,
                   fontcolor=fontcol, fontsize=fontsize))
}

makeCurve <- function(curve, col, lwd, lty, name) {
    controlPoints <- pointList(curve)
    bezierGrob(unit(sapply(controlPoints, "[" ,1), "native"),
               unit(sapply(controlPoints, "[" ,2), "native"),
	       name=name,
               gp=gpar(col=col, lwd=lwd, lty=lty))
}

drawCurve <- function(curve, col, lwd, lty) {
    grid.draw(curveGrob(curve, col, lwd, lty))
}

makeEdge <- function(edge, edgemode) {
    name <- paste(tail(edge), head(edge), sep="~")
    if (!length(edge@lwd))
        edge@lwd <- 1    
    if (!length(edge@lty))
        edge@lty <- "solid"
    splines <- splines(edge)
    n <- length(splines)
    curveNames <- paste("curve", name, seq(along=splines), sep="-")
    col <- color(edge)
    curves <- mapply(makeCurve, splines, col=col, lwd=edge@lwd, lty=edge@lty,
		 name=curveNames, SIMPLIFY=FALSE)
    
    # Edge label
    if (length(labelText(txtLabel(edge))) != 0) {
        fontcol <- labelColor(txtLabel(edge))
        if (length(fontcol) == 0) {
            fontcol <- "black"
        }
        fontsize <- labelFontsize(txtLabel(edge))
        label <- labelText(txtLabel(edge))
        xy <- labelLoc(txtLabel(edge))
        x <- unit(getX(xy), "native")
        y <- unit(getY(xy), "native")
        cex <- 1
        lab <- list(makeLabelGrob(label, x, y, fontcol, fontsize, cex, name))
    } else {
        lab <- list()
    }
        
    firstCP <- pointList(splines[[n]])[[1]]
    lastCP <- pointList(splines[[n]])[[4]]
    arrowsize <- as.numeric(arrowsize(edge))
    arrowhead <- arrowhead(edge)
    arrowtail <- arrowtail(edge)

    ## APOLOGY and EXPLANATION
    ##
    ## The "back" and "forward" arrows go through some contorted conditionals
    ## before making arrow ojects. This is mostly because of the way edge@sp and
    ## edge@ep are produced. A directed graph will have a sensible edge@ep for
    ## every edge with an arrow at the end (or edge@sp for an arrow at the
    ## start). This is used to plot the arrow between the last control point
    ## of an edge's splines and the border of the node.
    ##
    ## An undirected graph does not require edeg@sp or edge@ep values.
    ##
    ## Unfortunately Rgraphviz v.2.6.0 will produced meaningless edge@sp and
    ## edge@ep values on undirected graphs - values of "0,0", "NA,0" or "0,NA".
    ## For some as-yet-undiscovered reason it will also occasionally
    ## produce these types of values on a directed graph where no arrow is
    ## required.
    ##
    ## As a result, before attempting to draw an arrow the code will
    ## check for the following things:
    ##   1) is the edgemode of this graph undirected?
    ##   2) is the direction of this edge such that no arrow should be drawn
    ##      here? e.g. "back" arrow on "forward" edge
    ##   3) is the X-value in the sp/ep NA?
    ##   4) is the Y-value in the sp/ep NA?
    ##   5) is sp/ep equal to "0,0"?
    ## and if any of these is true an arrow is NOT drawn.
    ##
    ## Only if none of these conditions is true will an arrow be drawn.
    
    # "back" arrow    
    if (edgemode == "undirected" || edge@dir == "forward" ||
        (is.na(getX(sp(edge))) || is.na(getY(sp(edge)))) ||
        (getX(sp(edge)) == 0 && getY(sp(edge)) == 0)) {
      start <- list()
    } else {
      start <- list(makeArrowGrob(arrowtail, arrowsize,
                                  firstCP[1], firstCP[2], 
                                  getX(sp(edge)), getY(sp(edge)), 
                                  col, edge@lwd, edge@lty,
                                  name=paste("arrowtail", name, sep="-")))
    }
    
    # "forward" arrow
    if (edgemode == "undirected" || edge@dir == "back" ||
        (is.na(getX(ep(edge))) || is.na(getY(ep(edge)))) ||
        (getX(ep(edge)) == 0 && getY(ep(edge)) == 0)) {
      end <- list()
    } else {
      end <- list(makeArrowGrob(arrowhead, arrowsize,
                                lastCP[1], lastCP[2], 
                                getX(ep(edge)), getY(ep(edge)), 
                                col, edge@lwd, edge@lty,
                                name=paste("arrowhead", name, sep="-")))
    }
    
    gTree(children=do.call("gList", c(curves, start, end, lab)),
          name=paste("edge", name, sep="-"))
}

drawEdge <- function(edge, edgemode) {
    grid.draw(makeEdge(edge, edgemode))
}

grid.graph <- function(rag, newpage=FALSE, nodesOnTop=TRUE) {
    if (!is(rag, "Ragraph") || !laidout(rag))
        stop("Must have a laid out Ragraph object")
    if (newpage) {
        grid.newpage()
    }
    # (x, y) locations of all the nodes
    # The order is the same as the nodes in
    # the original graphNEL
    bb <- boundBox(rag)
    # Ensure aspect ratio
    pushViewport(viewport(width=unit(getX(upRight(bb)), "points"),
                          height=unit(getY(upRight(bb)), "points"),
                          layout=grid.layout(1, 1,
                              widths=(getX(upRight(bb)) - getX(botLeft(bb))) /
                                     (getY(upRight(bb)) - getY(botLeft(bb))),
                              respect=TRUE),
                          name="graphSizeVP"))
    pushViewport(viewport(layout.pos.col=1,
                          xscale=c(getX(botLeft(bb)), getX(upRight(bb))),
                          yscale=c(getY(botLeft(bb)), getY(upRight(bb))),
                          name="graphScaleVP"))
    if (nodesOnTop) {
        lapply(AgEdge(rag), drawEdge, edgemode(rag))
        lapply(AgNode(rag), drawNode)
    } else {
        lapply(AgNode(rag), drawNode)
        lapply(AgEdge(rag), drawEdge, edgemode(rag))
    }
    upViewport(2)
}
