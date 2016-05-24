plotVenn4d <-
function(x, labels = c('A','B','C','D'),
  Colors = c("red", "yellow", "green",'pink','darkgreen','blue','lightblue','tan', 'yellowgreen','orange','purple','white','grey','plum','brown'),
  Title = NULL, shrink = 1, rot=45)
{ # plot a 4-dimensional Venn diagram

#  suppressPackageStartupMessages(library(grid))

  ### Specify necessary functions

    xyangle <- function(x, y, directed = FALSE, deg = TRUE)
    {
        if (missing(y)) {
            y <- x[, 2]
            x <- x[, 1]
        }
        out <- atan2(y, x)
        if (!directed)
            out <- out%%pi
        if (deg)
            out <- out * 180/pi
        out
    }

    draw1ellipse <- function(x, y, a = 1, b = 1, angle = 0, segment = c(0,360),
        arc.only = TRUE, nv = 100, deg = TRUE, border = NULL, col = NA, lty = 1, lwd = 1, ...)
    {
        if (deg) {
            angle <- angle * pi/180
            segment <- segment * pi/180
        }
        z <- seq(segment[1], segment[2], length = nv + 1)
        xx <- a * cos(z)
        yy <- b * sin(z)
        alpha <- xyangle(xx, yy, directed = TRUE, deg = FALSE)
        rad <- sqrt(xx^2 + yy^2)
        xp <- rad * cos(alpha + angle) + x
        yp <- rad * sin(alpha + angle) + y
        if (!arc.only) {
            xp <- c(x, xp, x)
            yp <- c(y, yp, y)
        }
	
	      list(xp, yp)
    }

    getCommonEllipse <- function(ij.1, ij.2)
    {
      corners <- c(ij.1[1,1], ij.1[1,2], ij.2[1,1], ij.2[1,2])
      cornerCounts <- aggregate(corners, list(corners), length)
      cornerCounts[cornerCounts$x == 2, "Group.1"]
    }

    cornerAngle <- function(x, y, e)
    {
      (atan2((y - k[e])/ b[e], (x - h[e])/ a[e]) * 180 / pi + 360) %% 360
    }

    region.polygon <- function(number, region, midpoint, col, rot, ranges, doColor = TRUE)
    {
      tiltpi <- rot * pi / 180
      # order all corners counter clockwise
      ordering <- order(angle <- atan2(region$y - midpoint$y, region$x - midpoint$x), decreasing=FALSE)
      # duplicate the first to make a closed figure
      polygon <- rbind(polygon <- data.frame(region, angle, ordering)[ordering,], head(polygon,1))

      ellipseNumber.previous <- 0
      arcs <- NULL

      # get coordinates for all arcs
      for (i in 1:(nrow(polygon) - 1)) {
          # figure out which ellipse contains this segment
	        ellipseNumber <- getCommonEllipse(polygon[i,c("i","j")], polygon[i+1,c("i","j")])

          # if there is a loop segment, make sure that the ellipse switches to a new ellipse
	        if (length(ellipseNumber) > 1) ellipseNumber <- ellipseNumber[grep(ellipseNumber.previous,ellipseNumber, invert=TRUE)]

	        ellipseNumber.previous <- ellipseNumber
	        # get the starting and ending angles relative to the ellipse center
	        Angle1 <- cornerAngle(polygon[i,"x"], polygon[i,"y"], ellipseNumber)
	        Angle2 <- cornerAngle(polygon[i+1,"x"], polygon[i+1,"y"], ellipseNumber)
          # make sure we take the short way around the circumference
	        if (Angle2 - Angle1 > 200) Angle2 <- Angle2 - 360
	        if (Angle1 - Angle2 > 200) Angle1 <- Angle1 - 360

	        newarc <- draw1ellipse(h[ellipseNumber], k[ellipseNumber], a[ellipseNumber], b[ellipseNumber],
	            angle=0, segment=c(Angle1,Angle2), arc.only=TRUE, border=i, lwd=2)
	        # add new arc segments to the polygon
	        arcs <- rbind(arcs, cbind(newarc[[1]], newarc[[2]]))
      }
      if (tiltpi != 0) arcs <- arcs %*% matrix(c(cos(tiltpi), -sin(tiltpi), sin(tiltpi), cos(tiltpi)), byrow=TRUE, ncol=2)

      # Normalize arc values into range 0 to 1 - to fit into viewport
      arcs1 <- cbind((arcs[,1] + max(0, -ranges[1,1]))/(diff(ranges[,1])), (arcs[,2] + max(0, -ranges[1,2]))/(diff(ranges[,2])))

      if (doColor) grid.polygon(arcs1[,1], arcs1[,2], gp=gpar(fill=col))

      return(arcs1)
    }

  ### Initialize variables

    rot <- rot %% 360
    if (missing(x)) x <- paste(seq(15), Colors)

    if (length(names(x)) == 0) {
      names(x) <- c('1000', '0100', '1100', '0010', '1010', '0110', '1110', '0001', '1001', '0101', '1101', '0011', '1011', '0111', '1111')
    } else {
      ordering <- order(sapply(names(x), function(y) paste(rev(strsplit(y,"")[[1]]), collapse="")))
      x <- x[ordering]
      Colors <- Colors[ordering]
    }

    # Parameters for ellipses
    a <- c(10,10,5,5)
    b <- c(5,5,10,10)
    h <- c(0,0,0,5)
    k <- c(0,5,5,5)

    # Hard-coded crossover regions
    regions <- rbind(
      data.frame(i = 1, j = 2, x = c(rep(-8.66025,3), rep(8.66025, 3)), y = 2.5,
        TF = c('0100', '1000', '1100', '0101', '1001', '1101')),
      data.frame(i = 1, j = 3, x = c(rep(0.01184, 2),rep(4.988876, 4),rep(-4.988876,4)), y = c(rep(-4.99997,2),rep(4.3333,8)),
        TF = c('1000', '1010', '0101', '0111', '1101', '1111', '0100', '0110', '1100', '1110')),
      data.frame(i = 1, j = 4, x = c(3.6852,8,9.648,rep(0,4)), y = c(-4.648,-2.9987,1.3147, rep(5, 4)),
        TF = c('1001', '1001', '1001', '0110', '0111', '1110', '1111')),
      data.frame(i = 2, j = 3, x = c(rep(-4.472,4), rep(4.472,4), rep(-4.472,3),rep(4.472,4)), y = c(rep(0.52786,8), rep(9.472, 7)),
        TF = c('1000', '1010', '1100', '1110', '1001', '1011', '1101', '1111', '0010', '0100', '0110', '0001', '0011', '0101', '0111')),
      data.frame(i = 2, j = 4, x = c(rep(0.667,8), rep(10, 2)), y = c(rep(9.989, 4), rep(0.011, 4), rep(4.969,2)),
        TF = c('0010', '0011', '0110', '0111', '1010', '1011', '1110', '1111', '0001', '0101')),
      data.frame(i = 3, j = 4, x = 2.5, y = c(rep(-3.66025, 3),rep(13.66, 3)), TF = c('1001', '1010', '1011', '0001', '0010', '0011'))
    )

    # Midpoints for hard-coded crossover regions
    midpoints <- data.frame(
      x = c(-4.37352, -6.04042, -6.04042, -0.43516, -0.32341, -2.19859, -2.19859, 5.65737, 6.16114, 7.03031, 6.04042, 2.54627, 2.54627, 2.53192, 2.53192),
      y = c(-0.65737, 5.43516, 2.45373, 11.04042, -2.03031, 7.19859, 2.46808, 9.36716, -1.16074, 5.31864, 2.45373, 11.04042, -1.04042, 7.19859, 2.46808),
      TF = c('1000', '0100', '1100', '0010', '1010', '0110', '1110', '0001', '1001', '0101', '1101', '0011', '1011', '0111', '1111'))

    # Data for 4 main ellipses
    ellipses <- labelPars <- list()
    offset <- ifelse(rot < 15, 1, 1.5)
    for (i in 1:4) {
      ellipses[[i]] <- do.call(cbind, draw1ellipse(h[i], k[i], a[i], b[i]))
      labelPars[[i]] <- matrix(c(x = ifelse(i < 3, min(ellipses[[i]][,1])-offset, mean(ellipses[[i]][,1]) + offset),
        y = ifelse(i < 3, mean(ellipses[[i]][,2]) + .5, max(ellipses[[i]][,2]) + .5)), ncol = 2)

      if (rot != 0) {
        rotMatrix <- matrix(c(cos(tiltpi <- rot * pi / 180), -sin(tiltpi), sin(tiltpi), cos(tiltpi)), byrow=TRUE, ncol=2)
        ellipses[[i]] <- ellipses[[i]] %*% rotMatrix
        labelPars[[i]] <- labelPars[[i]] %*% rotMatrix
      }
    }

    # Ranges that will help normalize data for viewport
    rng <- apply(do.call(rbind, lapply(ellipses, function(i) apply(i, 2, range))), 2, range)

    # 'Normalize' labels' parameters' values to fit into viewport
    labelPars <- lapply(labelPars, function(x)
      cbind((x[,1] + max(0, -rng[1,1]))/(diff(rng[,1])), (x[,2] + max(0, -rng[1,2]))/(diff(rng[,2]))))

    # Align labels if rot is 45 or -45
    if (rot %in% c(45, 135))
    {
      labelPars[[2]][,2] <- labelPars[[3]][,2] <- max(c(labelPars[[2]][,2], labelPars[[3]][,2]))
      labelPars[[1]][,2] <- labelPars[[4]][,2] <- max(c(labelPars[[1]][,2], labelPars[[4]][,2]))
    }

  ### Plot

    # Plot title
    if (!is.null(Title))
      grid.text(Title, gp = gpar(fontsize=25*shrink, fontface="bold"),x = .5, y = 0.97)

    # Compute width & length that would most likely fit all the labels, without plotting outside boundaries.
    width <- ifelse(rot < 15, .96,.92)-.01 * sum(nchar(labels))
    height <- ifelse(is.null(Title),0.93,.9) + ifelse(rot > 180, 0.03, -0.03)
    pushViewport(viewport(x = 0.03 + .01 * sum(nchar(labels[1:2])), y = height,
      width = width, height = height - ifelse(rot < 180, .05, 0.1), just = c("left", "top")))

    # Plot 4 main ellipses
    for (i in seq_along(ellipses))
    { # 'Normalize' eclipse values to fit into viewport
      e1 <- cbind((ellipses[[i]][,1] + max(0, -rng[1,1]))/(diff(rng[,1])), (ellipses[[i]][,2] + max(0, -rng[1,2]))/(diff(rng[,2])))
      grid.polygon(x = e1[,1], y = e1[,2], gp = gpar(fill=(Colors[c(1,2,4,8)])[i]))
    }

    # Plot crossovers and fill-in values for each part
    for (ii in seq(nrow(midpoints)))
    { pars <- region.polygon(ii, regions[which(regions$TF %in% midpoints$TF[ii]),1:4],
        midpoints[which(midpoints$TF %in% midpoints$TF[ii]),], col=Colors[ii], rot, rng, doColor = !(ii %in% c(1,2,4,8)))
      grid.text(x[ii], x = mean(pars[,1]), y = mean(pars[,2]))
    }

    # Plot labels
    for (i in seq_along(ellipses))
      grid.text(labels[i], x = labelPars[[i]][,1], y = labelPars[[i]][,2], gp = gpar(fontsize = 18 * shrink, fontface = "bold"))

    upViewport()

}
