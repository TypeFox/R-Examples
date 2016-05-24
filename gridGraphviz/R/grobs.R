## create label grobs
makeLabelGrob <- function(label, x, y, col, fontsize, cex, name) {
    textGrob(label, x, y,
             gp=gpar(col=col, fontsize=fontsize, cex=cex),
             name=paste("label", name, sep="-"))
}

## create box grob for node gTree in node()
makeBoxGrob <- function(shape, name, x, y, height, lwidth, rwidth, color,
                        fillcolor) {
    a <- 0.5 * (lwidth + rwidth)
    b <- 0.5 * height
    boxName <- paste("box", name, sep="-")

    if (shape == "circle") {
        box <- circleGrob(x, y, r=a, name=boxName,
                          gp=gpar(col=color, fill=fillcolor))
    } else if (shape == "ellipse") { 
        angle <- seq(0, 2*pi, length=101)
        box <- polygonGrob(x + a*cos(angle),
                           y + b*sin(angle),
                           name=boxName,
                           gp=gpar(col=color, fill=fillcolor))
    } else if (shape == "polygon" || shape == "triangle" ||
               shape == "pentagon" || shape == "hexagon" ||
               shape == "septagon" || shape == "octagon"
               ) {
        vertices <- switch(shape,
                           polygon = 4,
                           triangle = 3,
                           pentagon = 5,
                           hexagon = 6,
                           septagon = 7,
                           octagon = 8)
        angle <- seq(0, 2*pi, length=vertices + 1)[-(vertices + 1)]
        if (vertices %% 2 != 0) {
            angle <- angle + pi/2
        } else {
            angle <- angle + pi/vertices + pi/2
            # expand polygon vertically to fill box
            b <- b*(1/cos(pi/vertices))
        }
        if (vertices %% 4 == 0)
            # expand polygon horizontally to fill box
            a <- a*(1/cos(pi/vertices))
        
        box <- polygonGrob(x + a*cos(angle), y + b*sin(angle),
                           name=boxName,
                           gp=gpar(col=color, fill=fillcolor))
    } else if (shape == "box" || shape == "rect" || shape == "rectangle") {
        box <- rectGrob(x, y,
                        width=lwidth + rwidth,
                        height=height,
                        name=boxName,         
                        gp=gpar(col=color, fill=fillcolor))
    } else if (shape == "square") {
        box <- rectGrob(x, y, width=height, height=height,
                        gp=gpar(col=color, fill=fillcolor))
    } else if (shape == "diamond") {
        xPoints <- unit.c(x - lwidth, x, x + rwidth, x)
        yPoints <- unit.c(y, y + b, y, y - b)
        box <- polygonGrob(xPoints, yPoints, name=boxName,
                           gp=gpar(col=color, fill=fillcolor))
    } else { # plain
        warning("Unsupported node shape; using 'box'")
        box <- rectGrob(x, y,
                        width=lwidth + rwidth,
                        height=height,
                        name=boxName, gp=gpar(col=NA, fill=NA))
    }
    box
}

## create arrow grob
makeArrowGrob <- function(arrowType, arrowsize, startX, startY, endX, endY, 
                          col, lwd, lty, name) {
  if (arrowType == "normal") arrowType <- "closed"
  if (arrowType == "vee") arrowType <- "open"
  if (arrowType == "none" || arrowType == "open" || arrowType == "closed") {
    arrow <- NULL
    if (arrowType != "none") {
      # FIXME: using a calculated 'length' does not scale when adjusted by
      # arrowsize, so graphviz default length 10 used
      arrowlen <- unit(arrowsize*10, "native")
      arrow <- arrow(angle=20, type=arrowType, length=arrowlen)
    }
    z <- segmentsGrob(startX, startY,
                      endX, endY,
                      default.units="native",
                      arrow=arrow, name=name,
                      gp=gpar(col=col, fill=col,
                              lwd=lwd, lty=lty))
  } else if (arrowType == "dot" || arrowType == "odot") {
    # FIXME: does not scale correctly compared to graphviz. documentation shows
    # graphviz sets default radius of 2, but does not work here. using 'length'
    # gets reasonably close results for arrowsize ~= 1
    dx <- endX - startX
    dy <- endY - startY
    length <- sqrt(dx^2 + dy^2) 
    theta <- atan2(dy, dx)
    r <- length/2 * arrowsize
    # FIXME: transparent circle will show edge inside if arrowsize > 1 due to
    # Rgraphviz not passing through arrow info at layout
    if (arrowType == "odot") fill <- "transparent" else fill <- col
    head <- circleGrob(endX - r*cos(theta), endY - r*sin(theta),
                       r=r, default.units="native", name=name,
                       gp=gpar(col=col, fill=fill, lwd=lwd, lty=lty))
    segment <- NULL
    # as different arrowheads cannot currently be passed through Rgraphviz's
    # agopen when laying out a line segment connects the edge to the dot head
    # when the arrowsize is less than the gap 'length'
    if (arrowsize < 1) {
      segment <- segmentsGrob(startX, startY,
                              grobX(name, theta*180/pi),
                              grobY(name, theta*180/pi),
                              default.units="native",
                              gp=gpar(col=col, fill=col,
                                      lwd=lwd, lty=lty))
    } 
    z <- gList(head, segment)
  } else if (arrowType == "box" || arrowType == "obox") {
    #FIXME: calculation of length, r, etc has same problems as with "dot" above
    dx <- endX - startX
    dy <- endY - startY
    length <- sqrt(dx^2 + dy^2) 
    theta <- atan2(dy, dx)
    r <- length/2 * arrowsize
    boxvp <- viewport(x=endX - r*cos(theta), y=endY - r*sin(theta),
                      width=r*2, height=r*2, default.units="native",
                      angle=theta*180/pi)
    if (arrowType == "obox") fill <- "transparent" else fill <- col
    head <- rectGrob(vp=boxvp, gp=gpar(col=col, fill=fill,
                                       lwd=lwd, lty=lty), name=name)
    segment <- NULL
    if (arrowsize < 1) {
      segment <- segmentsGrob(startX, startY,
                              grobX(name, theta*180/pi),
                              grobY(name, theta*180/pi),
                              default.units="native",
                              gp=gpar(col=col, fill=col,
                                      lwd=lwd, lty=lty))
    } 
    z <- gList(head, segment)
  }
}
