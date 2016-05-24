plotVenn2d  <-
function (x, labels = c('A', 'B'),
  Colors = c("red", "yellow", "green"),
  Title = NULL, shrink = 1, rot=0, radius= c(1,1), resizePlot = 1, reverseLabelOrdering=TRUE)
{ # plot a 2-dimensional Venn diagram

#  suppressPackageStartupMessages(library(grid))

  ### Specify necessary functions

    calcdist <- function(x, y) sqrt((x[1] - y[1])^2 + (x[2] - y[2])^2)
    calcangle <- function(x, y) atan((y[2] - x[2])/(y[1] - x[1]))

    getArcEnds <- function(center1, center2, radius) {
        centerDistance <- calcdist(center1, center2)
        connector <- ifelse(center1[1] > center2[1], pi, 0) +
            calcangle(center1, center2)
        distanceOne <- (radius[1]^2 - radius[2]^2 + centerDistance^2)/(2 * centerDistance)
        intersection <- acos(distanceOne/radius[1])
        return(c(begin = connector - intersection, end = connector + intersection))
    }

    arcPoints <- function(beginpt, endpt, center, radius) {
        angles <- seq(beginpt, endpt, length = nfacets)
        x <- center[1] + radius * cos(angles)
        y <- center[2] + radius * sin(angles)
        list(x = x, y = y)
    }

    getArc <- function(center1, center2, radius) {
        ends <- getArcEnds(center1, center2, radius)
        arcPoints(ends["begin"], ends["end"], center1, radius[1])
    }

    twoWayOverlap <- function(center1, center2, radius, color) {
        points1 <- getArc(center1, center2, radius[1:2])
        points2 <- getArc(center2, center1, radius[2:1])
        points <- list()
        points$x <- c(points1$x, points2$x)
        points$y <- c(points1$y, points2$y)
        grid.polygon(x = points$x, y = points$y, gp = gpar(fill = color))
    }

    centralArcs <- function(centers, i, radius) {
        if (i == 1) {
            j <- 2
            k <- 3
        }
        if (i == 2) {
            j <- 1
            k <- 3
        }
        if (i == 3) {
            j <- 2
            k <- 1
        }
        endsone <- getArcEnds(centers[[i]], centers[[j]], radius)
        endstwo <- getArcEnds(centers[[i]], centers[[k]], radius)
        if (endsone[2] < 0)
            endsone <- endsone + 2 * pi
        if (endstwo[2] < 0)
            endstwo <- endstwo + 2 * pi
        if (endstwo[1] < endsone[1] & endsone[1] < endstwo[2]) {
            return(arcPoints(endsone[1], endstwo[2], centers[[i]], radius))
        }
        else {
            return(arcPoints(endstwo[1], endsone[2], centers[[i]], radius))
        }
    }
    fromBase2 <- function(x) 4 * as.numeric(substr(x, 1, 1)) +
        2 * as.numeric(substr(x, 2, 2)) + as.numeric(substr(x, 3, 3))

    fromBase2_3 <- function(x)
        2 * as.numeric(substr(x, 1, 1)) + as.numeric(substr(x, 2, 2))

    intersectionCenter <- function(center) center[1] + (center[2] - center[1]) * radius[1] / (radius[1] + radius[2])

  ### Initialize variables

    if (is.null(rot)) rot <- 90
    if (is.null(radius)) radius <- c(1,1)
    radius <- radius * .25 * resizePlot

    if (is.null(names(x))) names(x) <- c("01", "10", "11")[seq(length(x))]
	xorder <- order(names(x))   # 4 lines of coded added 11/11/2013
	x <- x[xorder]
	Colors <- Colors[xorder]
	if (!reverseLabelOrdering) labels <- rev(labels)
	valptr <- unlist(lapply(names(x), fromBase2_3))
	
	if (class(x) != "character") {
		values <- rep(0, length(x))
		for (i in seq_along(x)) 
			if (valptr[i] %in% 1:3)
				values[valptr[i]] <- values[valptr[i]] + x[i]
	} else {
		values <- x
		for (i in seq_along(x)) 
			if (valptr[i] %in% 1:3)
				values[valptr[i]] <- x[i]		
	}
    nfacets <- 300

    r0 <- 0.18 * resizePlot
    dy <- r0 * sin(rot * pi / 180)
    dx <- r0 * cos(rot * pi / 180)
    laby0 <- ifelse(dy > 0, 1, -1) * (0.1 + radius)

    centers <- matrix(c(0.5 + dx, 0.5 + dy, 0.5 - dx, 0.5 - dy), byrow=TRUE, ncol=2)
    angle <- seq(0, 2 * pi, length = nfacets)[-nfacets]

  ### Plot

    if (!is.null(Title))
      grid.text(Title, gp = gpar(fontsize=25*shrink, fontface="bold"),x = 0.5, y = 0.97)
    
	centerDistance <- calcdist(centers[1,], centers[2,])

	if (radius[1] > centerDistance + radius[2]) {
	
	    grid.polygon(x = centers[1,1] + radius[1] * cos(angle), y = centers[1,2] + radius[1] * sin(angle),
            gp = gpar(fill = Colors[1]))

        grid.polygon(x = centers[2,1] + radius[2] * cos(angle), y = centers[2,2] + radius[2] * sin(angle),
            gp = gpar(fill = Colors[3]))
	
	} else if (radius[2] > centerDistance + radius[1]) {

        grid.polygon(x = centers[2,1] + radius[2] * cos(angle), y = centers[2,2] + radius[2] * sin(angle),
            gp = gpar(fill = Colors[2]))
	
	    grid.polygon(x = centers[1,1] + radius[1] * cos(angle), y = centers[1,2] + radius[1] * sin(angle),
            gp = gpar(fill = Colors[3]))
	
	} else {
	
	    grid.polygon(x = centers[1,1] + radius[1] * cos(angle), y = centers[1,2] + radius[1] * sin(angle),
            gp = gpar(fill = Colors[1]))

        grid.polygon(x = centers[2,1] + radius[2] * cos(angle), y = centers[2,2] + radius[2] * sin(angle),
            gp = gpar(fill = Colors[2]))

        if (radius[1] + radius[2] > centerDistance) twoWayOverlap(centers[1,], centers[2,], radius, Colors[3])
	
	}
	
    grid.text(labels[1], centers[1,1], centers[1,2] + laby0[1],
        gp = gpar(fontsize = 18 * shrink, fontface = "bold"))
    grid.text(labels[2], centers[2,1], centers[2,2] - laby0[2],
        gp = gpar(fontsize = 18 * shrink, fontface = "bold"))

    if (radius[1] + radius[2] > centerDistance)
		grid.text(values[3], intersectionCenter(centers[,1]), intersectionCenter(centers[,2]))
    grid.text(values[1], centers[1,1], centers[1,2])
    grid.text(values[2], centers[2,1], centers[2,2])

}
