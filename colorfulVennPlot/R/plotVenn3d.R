plotVenn3d <-
function(x, labels = c('A','B','C'),
  Colors = c("red","yellow","orange","lightblue","purple","green","grey"),
  Title = NULL, shrink=1, rot=0)
{ # plot a 3-dimensional Venn diagram

#  suppressPackageStartupMessages(library(grid))

  ### Specify necessary functions

    getArcEnds <- function(center1, center2, radius)
    {
      calcdist <- function(x,y) sqrt((x[1] - y[1])^2 + (x[2] - y[2])^2)
      calcangle <- function(x,y) atan((y[2] - x[2]) / (y[1] - x[1]))

      centerDistance <- calcdist(center1, center2)
      connector <- ifelse(center1[1] > center2[1], pi, 0) + calcangle(center1, center2)

      intersection <- acos((centerDistance / 2) / radius)

      c(begin=connector - intersection, end=connector + intersection)
    }

    arcPoints <- function(beginpt, endpt, center, radius)
    {
      angles <- seq(beginpt, endpt, length=nfacets)
      x <- center[1] + radius * cos(angles)
      y <- center[2] + radius * sin(angles)

      list(x=x,y=y)
    }

    getArc <- function(center1, center2, radius)
    {
      ends <- getArcEnds(center1, center2, radius)
      arcPoints(ends["begin"], ends["end"], center1, radius)
    }

    twoWayOverlap <- function(center1, center2, radius, color)
    {
      points1 <- getArc(center1, center2, radius)
      points2 <- getArc(center2, center1, radius)

      points <- list()
      points$x <- c(points1$x,points2$x)
      points$y <- c(points1$y,points2$y)

      grid.polygon(x=points$x, y=points$y, gp=gpar(fill=color))
    }

    centralArcs <- function(centers, i, radius)
    {
      if (i == 1) { j <- 2; k <- 3}
      if (i == 2) { j <- 1; k <- 3}
      if (i == 3) { j <- 2; k <- 1}
      endsone <- getArcEnds(centers[[i]], centers[[j]], radius)
      endstwo <- getArcEnds(centers[[i]], centers[[k]], radius)
      if (endsone[2] < 0) endsone <- endsone + 2 * pi
      if (endstwo[2] < 0) endstwo <- endstwo + 2 * pi
      if (endstwo[1] < endsone[1] & endsone[1] < endstwo[2]) {
        return (arcPoints(endsone[1], endstwo[2], centers[[i]], radius))
      } else {
        return (arcPoints(endstwo[1], endsone[2], centers[[i]], radius))
      }
    }

    fromBase2 <- function(x) 4 * as.numeric(substr(x,1,1)) + 2 * as.numeric(substr(x,2,2)) + as.numeric(substr(x,3,3))

    rotateCoordinates <- function(xy, tiltpi)
    {
      if (is.null(nrow(xy))) adjustCenter <- plotCenter
      else adjustCenter <- matrix(rep(plotCenter,nrow(xy)), byrow=TRUE, nrow=nrow(xy))
      rotationMatrix <- matrix(c(cos(tiltpi), -sin(tiltpi), sin(tiltpi), cos(tiltpi)), byrow=TRUE, ncol=2)
      (xy - adjustCenter) %*% rotationMatrix + adjustCenter
    }

  ### Initialize variables

    rot <- rot %% 360
    tiltpi <- rot * pi / 180
    if (missing(x)) x <- paste(seq(15), Colors)
    values <- rep(0, length(x))
    if (is.null(names(x)))
    names(x) <- c("101","110","010","001","011","100","111")[seq(length(x))]
    valptr <- unlist(lapply(names(x), fromBase2))

    inputColors <- Colors
      for (i in seq_along(x))
        if (valptr[i] %in% 1:7) {
      values[valptr[i]] <- values[valptr[i]] + x[i]
      Colors[valptr[i]] <- inputColors[i]
    }

    nfacets <- 300
    radius <- 0.25
    r0 <- .18
    dy <- r0 * sin(pi/6)
    dx <- r0 * cos(pi/6)
    labx0 <- .29
    labx1 <- 1.5
    labx2 <- .75
    plotCenter <- c(0.5, ifelse(is.null(Title),.5, .43))
    centers <- list(
      c(plotCenter[1], plotCenter[2] + r0),
      c(plotCenter[1] - dx, plotCenter[2] - dy),
      c(plotCenter[1] + dx, plotCenter[2] - dy))

    angle <- seq(0, 2*pi, length=nfacets) [-nfacets]

    points1 <- centralArcs(centers, 1, radius)
    points2 <- centralArcs(centers, 2, radius)
    points3 <- centralArcs(centers, 3, radius)

    points <- list()
    points$x <- c(points1$x, points2$x, points3$x)
    points$y <- c(points1$y, points2$y, points3$y)

   if (tiltpi != 0) {
    pointsxy <- rotateCoordinates(cbind(points$x, points$y), tiltpi)
    points$x <- pointsxy[,1]
    points$y <- pointsxy[,2]
    centers <- sapply(centers, function(x) list(rotateCoordinates(x, tiltpi)))
   }

  ### label and region-center locations

    labelLocations <- rotateCoordinates(matrix(c(
      plotCenter[1], plotCenter[2] + r0 + labx0,
      plotCenter[1] - dx, plotCenter[2] - dy - labx0,
      plotCenter[1] + dx, plotCenter[2] - dy - labx0
    ), byrow=TRUE, ncol=2), tiltpi)

    centerLocations <- rotateCoordinates(matrix(c(
      plotCenter[1] + dx * labx1, plotCenter[2] - dy * labx1,
      plotCenter[1] - dx * labx1, plotCenter[2] - dy * labx1,
      plotCenter[1], plotCenter[2] - r0 * labx2,
      plotCenter[1], plotCenter[2] + r0 * labx1,
      plotCenter[1] + dx * labx2, plotCenter[2] + dy * labx2,
      plotCenter[1] - dx * labx2, plotCenter[2] + dy * labx2,
      plotCenter[1], plotCenter[2]
    ), byrow=TRUE, ncol=2), tiltpi)

  ### Plot

    if (!is.null(Title))
      grid.text(Title, gp = gpar(fontsize=25*shrink, fontface="bold"),x = plotCenter[1], y = 0.97)

    grid.polygon(x=centers[[1]][1] + radius*cos(angle), y=centers[[1]][2] + radius*sin(angle), gp=gpar(fill=Colors[4]))
    grid.polygon(x=centers[[2]][1] + radius*cos(angle), y=centers[[2]][2] + radius*sin(angle), gp=gpar(fill=Colors[2]))
    grid.polygon(x=centers[[3]][1] + radius*cos(angle), y=centers[[3]][2] + radius*sin(angle), gp=gpar(fill=Colors[1]))

    twoWayOverlap(centers[[1]], centers[[2]], radius, Colors[6])
    twoWayOverlap(centers[[2]], centers[[3]], radius, Colors[3])
    twoWayOverlap(centers[[1]], centers[[3]], radius, Colors[5])

    grid.polygon(x=points$x, y=points$y, gp=gpar(fill=Colors[7]))

    for (i in 1:3) grid.text(labels[i], labelLocations[i,1], labelLocations[i,2], gp = gpar(fontsize=18*shrink, fontface="bold"))
    for (i in 1:7) grid.text(values[i], centerLocations[i,1], centerLocations[i,2])

}

