### $Id: plot.kohonen.R 23 2013-10-11 13:31:05Z ron.wehrens@gmail.com $
### Version 2.0.5: added parameter heatkeywidth (suggestion by Henning
### Rust). Especially useful for multiple plots in one figure.

### Parameter 'main' is mentioned explicitly since in most cases the
### default leaves open ugly space above the plot. We explicitly put
### it 1.2 units above the top row using 'text', if this is within the
### par("usr") range. Else, we use the standard 'title' command.

### Addition 04/11/07: keepMargins, codeRendering and whatmap arguments.
### Adapted for version 2.0: April 11, 2007
### Added quality plot, August 30 2007.
### Added default titles, August 31 2007.

"plot.kohonen" <- function (x,
                            type = c("codes", "changes", "counts",
                              "dist.neighbours", "mapping", "property",
                              "quality"),
                            classif = NULL,
                            labels = NULL, pchs = NULL, main = NULL,
                            palette.name = NULL, ncolors,
                            bgcol=NULL, zlim = NULL, heatkey = TRUE,
                            property, contin, whatmap = NULL,
                            codeRendering = NULL, keepMargins = FALSE,
                            heatkeywidth = .2, ...)
{
  type <- match.arg(type)

  switch(type,
         mapping = plot.kohmapping(x = x, classif = classif,
           main = main, labels = labels, pchs = pchs,
           bgcol = bgcol, keepMargins = keepMargins, ...),
         property = plot.kohprop(x = x, property, main = main,
           palette.name = palette.name, ncolors = ncolors,
           zlim = zlim, heatkey = heatkey,
           contin = contin, keepMargins = keepMargins, 
           heatkeywidth = heatkeywidth, ...),
         codes = plot.kohcodes(x = x, main = main,
           palette.name = palette.name, bgcol = bgcol,
           whatmap = whatmap, codeRendering = codeRendering,
           keepMargins = keepMargins, ...),
         quality = plot.kohquality(x = x, classif = classif, main = main,
           palette.name = palette.name, ncolors = ncolors,
           zlim = zlim, heatkey = heatkey, keepMargins = keepMargins, 
           heatkeywidth = heatkeywidth, ...),
         counts = plot.kohcounts(x = x, classif = classif, main = main,
           palette.name = palette.name, ncolors = ncolors,
           zlim = zlim, heatkey = heatkey, keepMargins = keepMargins,
           heatkeywidth = heatkeywidth, ...),
         changes = plot.kohchanges(x = x, main = main,
           keepMargins = keepMargins, ...),
         dist.neighbours = plot.kohUmatrix(x = x, main = main,
           palette.name = palette.name, ncolors = ncolors,
           zlim = zlim, heatkey = heatkey, keepMargins = keepMargins,
           heatkeywidth = heatkeywidth, ...))
}


### Overwrite the original plot.somgrid in the class library since
### that leaves open an ugly space at the top of the plot in case of
### hexagonal grids

### Unchanged in version 2.0

plot.somgrid <- function(x, xlim, ylim, ...)
{
  ## Following two lines leave equal amounts of space on both
  ## sides of the plot if no xlim or ylim are given
  if (missing(xlim)) xlim <- c(0, max(x$pts[,1]) + min(x$pts[,1]))
  if (missing(ylim)) ylim <-  c(max(x$pts[,2]) + min(x$pts[,2]), 0)
  eqscplot(xlim, ylim, axes = FALSE,
           type = "n", xlab = "", ylab = "", ...)
}


### Adapted for version 2.0: April 11.

plot.kohmapping <- function(x, classif, main, labels, pchs, bgcol,
                            keepMargins, ...)
{
  if (is.null(main)) main <- "Mapping plot"
  
  margins <- rep(0.6, 4)
  if (main != "") margins[3] <- margins[3] + 2
  if (!keepMargins) {
    opar <- par("mar")
    on.exit(par(mar = opar))
  }
  par(mar=margins)
    
  if (is.null(classif) & !is.null(x$unit.classif)) {
    classif <- x$unit.classif
  } else {
    if (is.list(classif) && !is.null(classif$unit.classif))
      classif <- classif$unit.classif
  }
  if (is.null(classif))
    stop("No mapping available")
  
  plot(x$grid, ...)
  title.y <- max(x$grid$pts[,2]) + 1.2
  if (title.y > par("usr")[4] - .2){
    title(main)
  } else {
    text(mean(range(x$grid$pts[,1])),
         title.y,
         main, adj = .5, cex = par("cex.main"),
         font = par("font.main"))
  }

  if (is.null(bgcol)) bgcol <- "transparent"
  symbols(x$grid$pts[, 1], x$grid$pts[, 2],
          circles = rep(0.5, nrow(x$grid$pts)),
          inches = FALSE, add = TRUE, bg = bgcol)
  if (is.null(labels) & !is.null(pchs))
    points(x$grid$pts[classif, 1] + rnorm(length(classif), 0, 0.12),
           x$grid$pts[classif, 2] + rnorm(length(classif), 0, 0.12),
           pch = pchs, ...)
  if (!is.null(labels))
    text(x$grid$pts[classif, 1] + rnorm(length(classif), 0, 0.12),
         x$grid$pts[classif, 2] + rnorm(length(classif), 0, 0.12),
         labels, ...)


  invisible()
}


### Adapted for version 2.0: April 11.

plot.kohprop <- function(x, property, main, palette.name, ncolors,
                         zlim, heatkey, contin, keepMargins,
                         heatkeywidth, ...)
{
  if (is.null(main)) main <- "Property plot"
  if (is.null(palette.name)) palette.name <- heat.colors
  
  margins <- rep(0.6, 4)
  if (heatkey) margins[2] <- margins[2] + 4
  if (main != "") margins[3] <- margins[3] + 2
  if (!keepMargins) {
    opar <- par("mar")
    on.exit(par(mar = opar))
  }
  par(mar = margins)
  
  plot(x$grid, ...)
  title.y <- max(x$grid$pts[,2]) + 1.2
  if (title.y > par("usr")[4] - .2){
    title(main)
  } else {
    text(mean(range(x$grid$pts[,1])),
         title.y,
         main, adj = .5, cex = par("cex.main"),
         font = par("font.main"))
  }
  
  if (is.null(zlim))
    zlim <- range(property, finite = TRUE)

  if (missing(ncolors)) 
    ncolors <- min(length(unique(property[!is.na(property)])), 20)
  bgcol <- palette.name(ncolors)

  bgcolors <- rep("gray", nrow(x$grid$pts))
  showcolors <- as.integer(cut(property,
                               seq(zlim[1], zlim[2],
                                   length = ncolors + 1),
                               include.lowest = TRUE))
  bgcolors[!is.na(showcolors)] <- bgcol[showcolors[!is.na(showcolors)]]

  symbols(x$grid$pts[, 1], x$grid$pts[, 2],
          circles = rep(0.5, nrow(x$grid$pts)), inches = FALSE,
          add = TRUE, fg = "black", bg = bgcolors)

  ## if contin, a pretty labelling of z colors will be used; if not,
  ## all colours will have their own label. The latter only if 
  ## property is a factor, unless explicitly given.
  if (missing(contin))
    contin <- !is.factor(property)
  
  if (heatkey) {
    if (length(unique(property)) < 10 & !contin) {
      plot.heatkey(x, zlim, bgcol, labels = levels(as.factor(property)),
                   contin = contin, heatkeywidth = heatkeywidth, ...)
    } else {
      plot.heatkey(x, zlim, bgcol, labels = NULL, contin = contin,
                   heatkeywidth = heatkeywidth, ...)
    }
  }

  invisible()
}


### Adapted for version 2.0, April 11.
### Whatmap argument left out because of trouble with indexing:
### changes is a matrix with a number of columns that is equal to the
### layers in whatmap; so there is no safe way to distinguish between,
### say, a whatmap of c(1,2,4) and 2:4 in a six-layer map where only
### the first four layers are used. Anyway.
### Checked: April 13.

plot.kohchanges <- function(x, main, keepMargins, ...)
{
  if (is.null(main)) main <- "Training progress"
  
  nmaps <- ncol(x$changes)
  
  ## check whether a legend is necessary and what names should be used
  if (nmaps > 1) {
    if (!is.null(colnames(x$changes))) {
      varnames <- colnames(x$changes)
    } else {
      varnames <- paste("Matrix", 1:ncol(x$changes))
    }
  }

  ## prepare a second y-axis in case of two maps
  if (nmaps == 2) {
    if (!keepMargins) {
      opar <- par("mar")
      on.exit(par(mar = opar))
    }
    par(mar=c(5.1, 4.1, 4.1, 4.1)) # axis scale to the right as well

    ## scale so that both have the same max value; assume only
    ## positive values.
    huhn <- x$changes
    huhn[,2] <- max(x$changes[,1]) * huhn[,2] / max(x$changes[,2])
    ticks <- pretty(x$changes[,2], length(axTicks(2)))
  } else {
    huhn <- x$changes
  }

  ## plot the plot!
  matplot(huhn, type = "l", lty = 1, main = main, 
          ylab = "Mean distance to closest unit", xlab = "Iteration", ...)
  abline(h=0, col="gray")

  ## plot the second axis
  if (nmaps == 2)
    axis(4, col.axis=2, at=ticks * max(x$changes[,1]) / max(x$changes[,2]),
         labels=ticks)

  ## plot the legend
  if (nmaps > 1)
    legend("topright", legend = varnames, lty=1, col = 1:nmaps, bty="n") 

  invisible()
}


### Adapted for version 2.0: April 11.
### Checked: April 13.

plot.kohcounts <- function(x, classif, main, palette.name, ncolors,
                           zlim, heatkey, keepMargins, heatkeywidth, ...)
{
  if (is.null(main)) main <- "Counts plot"
  if (is.null(palette.name)) palette.name <- heat.colors
  
  if (is.null(classif) & !is.null(x$unit.classif)) {
    classif <- x$unit.classif
  } else {
    if (is.list(classif) && !is.null(classif$unit.classif))
      classif <- classif$unit.classif
  }
  if (is.null(classif))
    stop("No mapping available")

  counts <- rep(NA, nrow(x$grid$pts))
  huhn <- table(classif)
  counts[as.integer(names(huhn))] <- huhn

  contin <- FALSE
  if (max(counts, na.rm = TRUE) > 10) contin <- TRUE
  
  plot.kohprop(x, property = counts, main = main,
               palette.name = palette.name, ncolors = ncolors,
               zlim = zlim, heatkey = heatkey, contin = contin,
               keepMargins = keepMargins, heatkeywidth = heatkeywidth, ...)

  invisible(counts)
}

### Introduced for version 2.0.5: Jan 16, 2009

plot.kohUmatrix <- function(x, classif, main, palette.name, ncolors,
                            zlim, heatkey, keepMargins, heatkeywidth, ...)
{
  if (x$method != "som" & x$method != "supersom")
    stop("Neighbour distance plot only implemented for (super)som")
  
  if (is.null(main)) main <- "Neighbour distance plot"
  if (is.null(palette.name)) palette.name <- heat.colors

  nhbrdist <- unit.distances(x$grid, x$toroidal)
  nhbrdist[nhbrdist > 1.05] <- NA
  if (x$method == "som") {
    for (i in 2:nrow(nhbrdist)) {
      for (j in 1:(i - 1)) {
        if (!is.na(nhbrdist[i,j]))
          nhbrdist[i,j] <- nhbrdist[j,i] <- dist(x$codes[c(i,j),])
      }
    }
  } else {
    if (x$method == "supersom") { # superfluous check, really
      nhbrdist[!is.na(nhbrdist)] <- 0
      for (k in 1:length(x$data)) {
        for (i in 2:nrow(nhbrdist)) {
          for (j in 1:(i - 1)) {
            if (!is.na(nhbrdist[i,j]))
              nhbrdist[i,j] <- nhbrdist[i,j] +
                x$weights[k] * dist(x$codes[[k]][c(i,j),])
          }
        }
        
        nhbrdist[j,i] <- nhbrdist[i,j]
      }
    }
  }

  neigh.dists <- colSums(nhbrdist, na.rm = TRUE)
  plot.kohprop(x, property = neigh.dists, main = main,
               palette.name = palette.name, ncolors = ncolors,
               zlim = zlim, heatkey = heatkey, contin = TRUE,
               keepMargins = keepMargins, heatkeywidth = heatkeywidth, ...)

  invisible(neigh.dists)
}

### Newly written for version 2.0, August 30 2007.
### Revised as a property plot: August 31 2007.

plot.kohquality <- function(x, classif, main, palette.name, ncolors,
                            zlim, heatkey, keepMargins, ...)
{
  if (is.null(main)) main <- "Distance plot"
  if (is.null(palette.name)) palette.name <- heat.colors

  distances <- NULL
  if (is.null(classif) & !is.null(x$unit.classif)) {
    classif <- x$unit.classif
    distances <- x$distances
  } else {
    if (is.list(classif) &&
        !is.null(classif$unit.classif) &&
        !is.null(classif$distances)) {
      classif <- classif$unit.classif
      distances <- classif$distances
    }
  }
  if (is.null(distances))
    stop("No mapping or mapping distances available")

  similarities <- rep(NA, nrow(x$grid$pts))
  hits <- as.integer(names(table(classif)))
  similarities[hits] <- sapply(split(distances, classif), mean)
                      
  plot.kohprop(x, property = similarities, main = main,
               palette.name = palette.name, ncolors = ncolors,
               zlim = zlim, heatkey = heatkey, contin = TRUE,
               keepMargins = keepMargins, ...)

  invisible(similarities)
}


### Adapted for version 2.0: April 11.
### Checked: April 13.
### New elements: whatmap, codeRendering, keepMargins, legend
### Added palette.name for version 2.0.6. Aug 3, 2010.

plot.kohcodes <- function(x, main, palette.name, bgcol, whatmap,
                          codeRendering, keepMargins, ...)
{
  if (!keepMargins) {
    opar <- par(c("mar", "ask"))
    on.exit(par(opar))
  }

  if (is.null(palette.name)) palette.name <- terrain.colors
  
  whatmap <- check.whatmap(x, whatmap)
  nmaps <- length(whatmap)
  
  ## check if x$codes is a list; if so, call this function for every
  ## list element separately.
  if (is.list(x$codes)) {
    if (prod(par("mfrow")) < nmaps) par(ask = TRUE)

    for (i in 1:nmaps) {
      ## make a new object that only has one set of codebook vectors
      huhn <- x
      huhn$codes <- huhn$codes[[whatmap[i]]]

      ## allow a different title for every plot
      if (length(main) == length(x$codes)) {
        main.title <- main[whatmap[i]]
      } else {
        if (length(main) == nmaps) {
          main.title <- main[i]
        } else {
          if (length(main) == 1) {
            main.title <- main
          } else {
            if (is.null(main)) {
              if (!is.null(names(x$codes))) {
                main.title <- names(x$codes)[whatmap[i]]
              } else {
                main.title <- "Codes plot"
              }
            }
          }
        }
      }

      ## allow a different codeRendering for every plot
      if (length(codeRendering) == length(x$codes)) {
        cR <- codeRendering[whatmap[i]]
      } else {
        if (length(codeRendering) == nmaps) {
          cR <- codeRendering[i]
        } else {
          cR <- codeRendering
        }
      }

      plot.kohcodes(huhn, main = main.title, palette.name = palette.name,
                    bgcol=bgcol, whatmap = NULL,
                    codeRendering = cR, keepMargins = TRUE, ...)
    }
  } else {
    codes <- x$codes
    nvars <- ncol(codes)
    
    maxlegendcols <- 3  ## nr of columns for the legend
    if (is.null(codeRendering))  ## use default
      codeRendering <- ifelse(nvars < 15, "segments", "lines")
    
    margins <- rep(0.6, 4)  # no text annotation anywhere
    if (!is.null(main))
      margins[3] <- margins[3] + 2
    par(mar = margins)
    
    if (codeRendering == "segments" & # we need space for the legend here...
        ##        nvars < 15 &
        !is.null(colnames(codes))) {
      plot(x$grid, 
           ylim = c(max(x$grid$pts[,2]) + min(x$grid$pts[,2]), -2))
      current.plot <- par("mfg")
      plot.width <- diff(par("usr")[1:2])

      cex <- 1 # First see if the legend fits
      leg.result <- legend(x = mean(x$grid$pts[,1]), xjust = 0.5,
                           y = 0, yjust = 1,
                           legend = colnames(codes),
                           cex=cex, plot=FALSE,
                           ncol = min(maxlegendcols, nvars),
                           fill = palette.name(nvars))
      while (leg.result$rect$w > plot.width) {
        cex <- cex*0.9 # if too large, decrease text size
        leg.result <- legend(x = mean(x$grid$pts[,1]), xjust = 0.5,
                             y = 0, yjust = 1,
                             legend = colnames(codes),
                             cex=cex, plot=FALSE,
                             ncol = min(maxlegendcols, nvars),
                             fill = palette.name(nvars))
      } # until it fits!

      leg.result <- legend(x = mean(x$grid$pts[,1]), xjust = 0.5,
                           y = 0, yjust = 1, cex=cex,
                           legend = colnames(codes), plot=FALSE,
                           ncol = min(maxlegendcols, nvars),
                           fill = palette.name(nvars), ...)

      par(mfg = current.plot)
      plot(x$grid, 
           ylim = c(max(x$grid$pts[,2]) + min(x$grid$pts[,2]),
             -leg.result$rect$h))

      legend(x = mean(x$grid$pts[,1]), xjust = 0.5,
             y = 0, yjust = 1, cex=cex, plot = TRUE,
             legend = colnames(codes),
             ncol = min(maxlegendcols, nvars),
             fill = palette.name(nvars), ...)
    } else {
      plot(x$grid, ...)
    }
    
    title.y <- max(x$grid$pts[,2]) + 1.2
    if (title.y > par("usr")[4] - .2){
      title(main)
    } else {
      text(mean(range(x$grid$pts[,1])),
           title.y,
           main, adj = .5, cex = par("cex.main"),
           font = par("font.main"))
    }
    
    if (is.null(bgcol)) bgcol <- "transparent"
    symbols(x$grid$pts[, 1], x$grid$pts[, 2],
            circles = rep(0.5, nrow(x$grid$pts)), inches = FALSE,
            add = TRUE, bg = bgcol)
    
    if (codeRendering == "lines") {
      yrange <- range(codes)
      codes <- codes - mean(yrange)
    } else {
      codemins <- apply(codes, 2, min)
      codes <- sweep(codes, 2, codemins)
    }

    switch(codeRendering,
           segments = {
             stars(codes, locations = x$grid$pts,
                   labels = NULL, len = 0.4,
                   add=TRUE, col.segments=palette.name(nvars),
                   draw.segments=TRUE)
           },             
           lines = {
             for (i in 1:nrow(x$grid$pts)) { # draw baseline
               if (yrange[1]<0 & yrange[2] > 0) {
                 lines(seq(x$grid$pts[i, 1] - 0.4,
                           x$grid$pts[i, 1] + 0.4,
                           length = 2),
                       rep(x$grid$pts[i, 2], 2),
                       col = "gray")
               }
               lines(seq(x$grid$pts[i, 1] - 0.4,
                         x$grid$pts[i, 1] + 0.4,
                         length = ncol(codes)),
                     x$grid$pts[i, 2] + codes[i, ] * 0.8/diff(yrange),
                     col = "red")
             }
           },
           stars = stars(codes, locations = x$grid$pts,
             labels = NULL, len = 0.4, add=TRUE)
           )
    
  }

  invisible()
}


### Added heatkeywidth parameter in version 2.0.5 (contribution by
### Henning Rust)

plot.heatkey <- function (x, zlim, bgcol, labels, contin, heatkeywidth, ...)  
{
  ncolors <- length(bgcol)
  
  yrange <- range(x$grid$pts[, 2])
  smallestx <- min(x$grid$pts[,1])
  ## A width of .2 looks OK on my screen
  xleft <- c(smallestx - heatkeywidth, smallestx) - 1
  yleft <- seq(yrange[1] - 0.5,
               yrange[2] + 0.5,
               length = ncolors + 1)
  rect(xleft[1], yleft[1:ncolors],
       xleft[2], yleft[2:(ncolors + 1)],
       border = "black", col = bgcol,
       xpd = TRUE)

  cex <- list(...)$cex

  if (contin) {
    zvals <- pretty(zlim)
    zvals <- zvals[zvals <= max(zlim) & zvals >= min(zlim)]
    yvals <- yrange[1] - .5 + (diff(yrange) + 1)*(zvals - zlim[1])/diff(zlim)

    text(xleft[2] - 1.3*diff(xleft),
         yvals,
         formatC(zvals),
         xpd=TRUE, adj=1, cex=cex)
  } else {
    if (is.null(labels))
      labels <- 1:ncolors
    
    text(xleft[2] - 1.3 * diff(xleft),
         yleft[-1] - 0.5*diff(yleft[1:2]),
         sort(labels),
         xpd = TRUE, adj=1, cex=cex)
  }
}

### Show cluster boundaries additional to one of the map plots
### Additional arguments may be col. Based on code from Leo Lopes.
### Oct 8 - rewritten the function:
### * neighbours of different classes are much more easily found, and
###   now are correct, too, also for toroidal maps (bug fix)
### * additional lines on top of the map are drawn if necessary for
###   toroidal maps (bug fix)
### * in rectangular maps some superfluous lines are no longer drawn
###   (purely esthetic improvement)
### Dec 10, 2013: additional bug fix (bug noted by Thomas Campagne),
### and an additional change for toroidal maps: boundaries at the
### outside of the map are now drawn even when they are already
### present at the other end of the map (easier interpretation)

add.cluster.boundaries <- function(x, clustering, lwd = 5, ...)
{
  grd <- x$grid
  if (x$toroidal) {
    ydiff <- diff(grd$pts[1 + c(0, grd$xdim),2])

    botrow <- 1:grd$xdim
    toprow <- grd$xdim*grd$ydim + 1 - (grd$xdim:1)
    rightcol <- (1:grd$ydim)*grd$xdim
    leftcol <- (1:grd$ydim)*grd$xdim + 1 - grd$xdim
    
    newpts <- rbind(cbind(grd$pts[botrow, 1], max(grd$pts[,2]) + ydiff),
                    cbind(grd$pts[toprow, 1], min(grd$pts[,2]) - ydiff),
                    cbind(grd$pts[leftcol, 1] - 1, grd$pts[leftcol, 2]),
                    cbind(grd$pts[rightcol, 1] + 1, grd$pts[rightcol,2]))
    cluster <- c(clustering, clustering[c(botrow, toprow, rightcol, leftcol)])
    if (x$grid$topo == "hexagonal") {
      ## we need to add two extra points, one in the bottom right, the
      ## other top left - explicitly add the clustering of these two
      ## points, too. These are the cluster ids of the top left and
      ## bottom right corners of the _original_ map
      newpts <- rbind(newpts,
                      c(grd$pts[toprow[grd$xdim],1]+1,
                        min(grd$pts[,2]) - ydiff),
                      c(grd$pts[botrow[1], 1]-1,
                        max(grd$pts[,2]) + ydiff))
      cluster <- c(cluster, clustering[toprow[1]], clustering[grd$xdim])
    }
    
    grd$pts <- rbind(grd$pts, newpts)
  } else {
    cluster <- clustering
  }
  
  nhbrdist <- unit.distances(grd, FALSE) ## new grd is treated as non-toroid
  nhbrdist[col(nhbrdist) >= row(nhbrdist)] <- 2
  neighbours <- which(nhbrdist > .95 & nhbrdist < 1.05, arr.ind = TRUE)

  diffclass.idx <-
      sapply(1:nrow(neighbours),
             function(ii)
             cluster[neighbours[ii, 1]] != cluster[neighbours[ii, 2]])
  neighbours <- neighbours[diffclass.idx,]
  ## final step: remove rows in neighbours that are completely outside the
  ## original grid (only relevant for the toroidal case)
  if (x$toroidal) {
    idx <- apply(neighbours, 1, function(x) all(x > grd$xdim*grd$ydim))
    neighbours <- neighbours[!idx,]
  }
  
  ## Function to actually plot the boundaries. For clarity, we 
  ## draw boundaries at the edges on both sides of the map, which is
  ## achieved simply by ignoring double lines - just plot'em all.
  plot.hex.boundary <- function(nb, grd, lwd, ...) {
    radius <- .5/cos(pi/6)             ## horizontal unit distance always 1

    ## for debugging...
    ## text(grd$pts, labels = 1:nrow(grd$pts))
    ## browser()
    for (i in 1:nrow(nb)) {
      u1 <- nb[i,1]
      u2 <- nb[i,2]
      
      dloc <- grd$pts[u1,] - grd$pts[u2,]

      if (abs(dloc[2]) < .1) {         # vertical line segments
        angle <- pi                    # left
        if (dloc[1] > .9) angle <- 0   # right
      } else {
        if (dloc[2] > .1) {            
          angle <- pi/3                # NE
          if (dloc[1] < -.1)            
              angle <- 2*pi/3          # NW
        } else {                       # dloc[2] < -.1
          if (dloc[1] > .1) {
            angle <- -pi/3             # SE
          } else {
            angle <- -2*pi/3           # SW
          }
        }
      }
              
      segments(grd$pts[u2,1]+radius*cos(angle-pi/6),
               grd$pts[u2,2]+radius*sin(angle-pi/6), 
               grd$pts[u2,1]+radius*cos(angle+pi/6),
               grd$pts[u2,2]+radius*sin(angle+pi/6),
               lwd = lwd, xpd = NA, ...)  
    }
  }

  plot.rect.boundary <- function(nb, grd, ...) {
    verticals <- which(apply(nb[,2:1],
                             1,
                             function(idx)
                             diff(grd$pts[idx,1]) == 1 &
                             diff(grd$pts[idx,2]) == 0))
                       
    for (i in verticals) {
      segments(x0 = mean(grd$pts[nb[i,],1]),
               y0 = grd$pts[nb[i,1],2] - .5,
               x1 = mean(grd$pts[nb[i,],1]),
               y1 = grd$pts[nb[i,1],2] + .5,
               ...)
    }

    horizontals <- which(apply(nb[,2:1],
                               1,
                               function(idx)
                               diff(grd$pts[idx,2]) == 1 &
                               diff(grd$pts[idx,1]) == 0))
    for (i in horizontals) {
      segments(x0 = grd$pts[nb[i,1],1] - .5,
               y0 = mean(grd$pts[nb[i,],2]),
               x1 = grd$pts[nb[i,1],1] + .5,
               y1 = mean(grd$pts[nb[i,],2]),
               ...)
    }

  }

  switch(grd$topo,
         rectangular = 
         plot.rect.boundary(neighbours, grd, lwd = lwd, ...),
         plot.hex.boundary(neighbours, grd, lwd = lwd, ...))

  invisible()
}


identify.kohonen <- function(x, ...) {
  ## map units have a radius of 1, so this is the tolerance we would
  ## like to have when pointing at map units
  tol <- par("pin")[1] / diff(par("usr")[1:2])
  identify(x$grid$pts[,1], x$grid$pts[,2], 1:nrow(x$grid$pts),
           tolerance = tol, ...)
}
