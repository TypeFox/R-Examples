"plot.wccsom" <- function (x,
                           type = c("changes", "codes", "counts",
                             "mapping", "prediction", "property",
                             "quality"),
                           classif=NULL, labels = NULL,
                           pchs = NULL, main = NULL, 
                           palette.name = heat.colors, ncolors,
                           unit.colors, unit.bgcol = NULL,
                           zlim = NULL, property = NULL,
                           heatkey = TRUE, contin, ...)
{
  type <- match.arg(type)
  
  switch(type,
         prediction = plot.wccpred(x, property, main, palette.name,
           ncolors, unit.colors, zlim, heatkey, labels, contin, ...),
         mapping = plot.wccmapping(x, classif, main, labels, pchs,
           unit.bgcol, ...),
         property = plot.wccprop(x, property, main, palette.name,
           ncolors, unit.colors, unit.bgcol, zlim, heatkey, contin, ...),
         codes = plot.wcccodes(x, main, unit.bgcol, ...),
         counts = plot.wcccounts(x, classif, main, palette.name,
           ncolors, unit.colors, unit.bgcol, zlim, heatkey, ...),
         changes = plot.wccchanges(x, main, ...),
         quality = plot.wccquality(x, classif, main, palette.name,
           ncolors, unit.colors, unit.bgcol, zlim, heatkey, ...))

  invisible()
}


### Overwrite the original plot.somgrid in the class library since
### that leaves open an ugly space at the top of the plot in case of
### hexagonal grids

"plot.somgrid" <- function(x, type = "p", xlim, ylim, ...)
{
  ## Following two lines leave equal amounts of space on both
  ## sides of the plot if no xlim or ylim are given
  if (missing(xlim)) xlim <- c(0, max(x$pts[,1]) + min(x$pts[,1]))
  if (missing(ylim)) ylim <-  c(max(x$pts[,2]) + min(x$pts[,2]), 0)
  eqscplot(xlim, ylim,                
           axes = FALSE, type = "n", xlab = "", ylab = "", ...)
  invisible()
}



### No unit.bgcol argument because there should not be any NAs

"plot.wccpred" <- function(x, Y, main, palette.name, ncolors,
                           unit.colors, zlim, heatkey,
                           labels, contin, ...)
{
  if (is.null(x$predict.type))
    stop("Prediction plot only available for supervised SOMs...")
  
  if (is.null(Y)) {
    if (is.null(x$codeYs))
      stop("No predictions available")
    
    if (x$predict.type == "continuous") {
      Y <- x$codeYs
      
      if (ncol(Y) > 1) {       # make separate plots for each variable, use
                                        # drop=FALSE to keep column names
        for (i in 1:ncol(Y))
          plot.wccpred(x, Y[,i,drop=FALSE], main,
                       palette.name, ncolors, unit.colors, 
                       zlim, heatkey, labels, contin, ...)
        
        return()
      }
    } else { # classification
      Y <- classmat2classvec(x$codeYs)
    }
  }
  
  if (missing(contin))
    contin <- (x$predict.type == "continuous")
  if (!contin) Y <- as.factor(Y)

  if (!missing(unit.colors)) {
    ncolors <- length(unit.colors)
  } else {
    if (missing(ncolors)) {
      ncolors <- ifelse(contin, 20, min(nlevels(factor(Y)), 20))
    }
    unit.colors <- palette.name(ncolors)
  }
  
   margins <- rep(0.6, 4)
  if (heatkey) 
    margins[2] <- margins[2] + 4
  if (is.null(main))
    main <- colnames(Y)
  margins[3] <- margins[3] + 2

  opar <- par("mar")
  on.exit(par(mar = opar))
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
  
  if (is.null(zlim)) {
    if (contin) {
      zlim <- range(Y, finite = TRUE)
    } else {
      zlim <- c(1, nlevels(Y))
    }
  }

  symbols(x$grid$pts[, 1], x$grid$pts[, 2],
          circles = rep(0.5, nrow(x$grid$pts)), inches = FALSE,
          add = TRUE, fg = "black",
          bg = unit.colors[as.integer(cut(as.numeric(Y), seq(zlim[1], zlim[2],
            length = ncolors + 1), include.lowest = TRUE))])
  
  if (is.null(labels))
    labels <- levels(Y)
  
  if (heatkey) {
    plot.heatkey(x, zlim, unit.colors, labels, contin = contin, ...)
  }
}


plot.wccmapping <- function(x, classif, main, labels, pchs,
                            unit.bgcol, ...)
{
  opar <- par("mar")
  on.exit(par(mar = opar))
  ifelse(is.null(main),
         par(mar = c(0.6, 0.6, 0.6, 0.6)),
         par(mar = c(0.6, 0.6, 2.6, 0.6)))
    
  if (is.null(classif) & !is.null(x$unit.classif)) {
    classif <- x$unit.classif
  } else {
    if (is.list(classif) && !is.null(classif$classif))
      classif <- classif$classif
  }
  if (is.null(classif))
    stop("no classif argument")

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
  
  if (is.null(unit.bgcol)) unit.bgcol <- "gray"
  symbols(x$grid$pts[, 1], x$grid$pts[, 2],
          circles = rep(0.5, nrow(x$grid$pts)),
          inches = FALSE, add = TRUE, bg = unit.bgcol)
  if (is.null(labels)) {
    if (is.null(pchs)) pchs <- 1
    points(x$grid$pts[classif, 1] + rnorm(length(classif), 0, 0.12),
           x$grid$pts[classif, 2] + rnorm(length(classif), 0, 0.12),
           pch = pchs, ...)
  }
  if (!is.null(labels))
    text(x$grid$pts[classif, 1] + rnorm(length(classif), 0, 0.12),
         x$grid$pts[classif, 2] + rnorm(length(classif), 0, 0.12),
         labels, ...)
}


"plot.wccprop" <- function(x, property, main, palette.name, ncolors,
                           unit.colors, unit.bgcol, zlim, heatkey,
                           contin, ...)
{
  margins <- rep(0.6, 4)
  if (heatkey) 
    margins[2] <- margins[2] + 4
  if (!is.null(main)) 
    margins[3] <- margins[3] + 2
  opar <- par("mar")
  on.exit(par(mar = opar))
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

  ## if contin, a pretty labelling of z colors will be used; if not,
  ## all colours will have their own label. The latter only if the
  ## number of categories is smaller than 10, unless explicitly
  ## given.
  if (missing(contin))
    contin <- !(length(unique(property)) < min(10, nrow(x$grid$pts)))
  
  if (!contin) property <- factor(property)
  if (is.factor(property)) {
    contin <- FALSE
    labels <- levels(property)
    property <- as.integer(property)
  } else {
    labels <- NULL
  }
  
  if (is.null(zlim))
    zlim <- range(property, finite = TRUE)

  if (!missing(unit.colors)) {
    ncolors <- length(unit.colors)
  } else {
    if (missing(ncolors)) {
      ncolors <- ifelse(contin, 20, min(nlevels(factor(property)), 20))
    }
    unit.colors <- palette.name(ncolors)
  }

  if (is.null(unit.bgcol)) unit.bgcol <- "gray"
  bgcolors <- rep(unit.bgcol, nrow(x$grid$pts))
  showcolors <- as.integer(cut(property,
                               seq(zlim[1], zlim[2],
                                   length = ncolors + 1),
                               include.lowest = TRUE))
  bgcolors[!is.na(showcolors)] <- unit.colors[showcolors[!is.na(showcolors)]]

  symbols(x$grid$pts[, 1], x$grid$pts[, 2],
          circles = rep(0.5, nrow(x$grid$pts)), inches = FALSE,
          add = TRUE, fg = "black", bg = bgcolors)

  if (heatkey)
    plot.heatkey(x, zlim, unit.colors, labels = labels, contin = contin, ...)
}


"plot.wccchanges" <- function(x, main, ...)
{
  if (is.matrix(x$changes)) { # for supervised networks
    opar <- par("mar")
    on.exit(par(mar = opar))
    par(mar=c(5.1, 4.1, 4.1, 4.1)) # axis scale to the right as well
    
    ## scale so that both have the same max value; assume only
    ## positive values.
    huhn <- x$changes
    huhn[,2] <- max(x$changes[,1]) * huhn[,2] / max(x$changes[,2])
    ticks <- pretty(x$changes[,2], length(axTicks(2)))
    
    matplot(huhn, type = "l", lty = 1, col=c(1,2), main = main, 
            ylab = "Mean change", xlab = "Iteration", ...)
    axis(4, col.axis=2, at=ticks * max(x$changes[,1]) / max(x$changes[,2]),
         labels=ticks)
    legend("topright", legend = c("X update", "Y update"),
           lty=c(1,1), col=c(1,2), bty="n") 
  } else {
    plot(x$changes, type = "l", ylab = "Mean change", main = main,
         xlab = "Iteration", ...)
  }
}


"plot.wcccounts" <- function(x, classif, main,
                             palette.name, ncolors, unit.colors,
                             unit.bgcol, zlim, heatkey, ...) {
  margins <- rep(0.6, 4)
  if (heatkey) 
    margins[2] <- margins[2] + 4
  if (!is.null(main)) 
    margins[3] <- margins[3] + 2
  opar <- par("mar")
  on.exit(par(mar = opar))
  par(mar = margins)

  if (is.null(classif) & !is.null(x$unit.classif)) {
    classif <- x$unit.classif
  } else {
    if (is.list(classif) && !is.null(classif$classif))
      classif <- classif$classif
  }
  if (is.null(classif))
    stop("no classif argument")
  
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

  if (is.null(unit.bgcol)) unit.bgcol <- "gray"
  bgcolors <- rep(unit.bgcol, nrow(x$grid$pts))
  hits <- as.integer(names(table(classif)))
  counts <- as.integer(table(classif))
  if (is.null(zlim))
    zlim <- c(0, max(counts))

  if (!missing(unit.colors)) {
    ncolors <- length(unit.colors)
  } else {
    if (missing(ncolors)) {
      ncolors <- min(max(counts), 20)
    }
    unit.colors <- palette.name(ncolors)
  }

  showcolors <- as.integer(cut(counts,
                               seq(zlim[1], zlim[2],
                                   length = ncolors + 1),
                               include.lowest = TRUE))
  bgcolors[hits] <- unit.colors[showcolors]
  
  symbols(x$grid$pts[, 1], x$grid$pts[, 2],
          circles = rep(0.5, nrow(x$grid$pts)), inches = FALSE,
          add = TRUE, fg = "black", bg = bgcolors)
  
  if (heatkey) plot.heatkey(x, zlim, unit.colors, contin = TRUE, ...)
}



"plot.wcccodes" <- function(x, main, unit.bgcol, ...)
{
  opar <- par("mar")
  on.exit(par(mar = opar))
  ifelse(is.null(main),
         par(mar = c(0.6, 0.6, 0.6, 0.6)),
         par(mar = c(0.6, 0.6, 2.6, 0.6)))
  
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

  if (is.null(unit.bgcol)) unit.bgcol <- "transparent"
  symbols(x$grid$pts[,1], x$grid$pts[,2],
          circles=rep(.5, nrow(x$grid$pts)),
          inches=FALSE, add=TRUE, bg=unit.bgcol)
  for (i in 1:nrow(x$grid$pts)) {
    lines(seq(x$grid$pts[i,1]-.4,
              x$grid$pts[i,1]+.4,
              length=ncol(x$codes)),
          x$grid$pts[i,2] - .2 + x$codes[i,]*.5/max(x$codes[i,]),
          col="red")
  }
}


"plot.wccquality" <- function(x, classif, main, palette.name, ncolors,
                              unit.colors, unit.bgcol, zlim, heatkey, ...)
{
  if (is.null(classif) & !is.null(x$unit.classif) & !is.null(x$wccs)) {
    classif <- x$unit.classif
    wccs <- x$wccs
  } else {
    if (is.list(classif) &&
        !is.null(classif$classif) &&
        !is.null(classif$wccs)) {
      classif <- classif$classif
      wccs <- classif$wccs
    }
  }
  if (is.null(classif) | is.null(wccs))
    stop("no classif argument")

  nc <- nrow(x$codes)
  
  sdwccs <- wccs2 <- rep(NA, nc)
  for (i in 1:nc) {
    allis <- which(classif == i)
    if (length(allis) > 0) {
      wccs2[i] <- mean(wccs[allis])
      if (length(allis) > 1) 
        sdwccs[i] <- sd(wccs[allis])
    }
  }
  
  if (is.null(zlim))
    zlim=range(wccs2, na.rm=TRUE)
  
  if (!missing(unit.colors)) {
    ncolors <- length(unit.colors)
  } else {
    if (missing(ncolors)) {
      ncolors <- 20
    }
    unit.colors <- palette.name(ncolors)
  }

  bgcolors <- unit.colors[as.integer(cut(wccs2,
                                         seq(zlim[1], zlim[2],
                                             length=ncolors+1),
                                         include.lowest=TRUE))]
  if (is.null(unit.bgcol)) unit.bgcol <- "gray"
  bgcolors[is.na(wccs2)] <- unit.bgcol
  
  margins <- rep(0.6, 4)
  if (heatkey) 
    margins[2] <- margins[2] + 4
  
  margins[3] <- margins[3] + 2

  opar <- par("mar")
  on.exit(par(mar = opar))
  par(mar = margins)
  
  plot(x$grid, ...)
  if (!is.null(main)) {
    title.y <- max(x$grid$pts[,2]) + 1.2
    if (title.y > par("usr")[4] - .2){
      title(main)
    } else {
      text(mean(range(x$grid$pts[,1])),
           title.y,
           main, adj = .5, cex = par("cex.main"),
           font = par("font.main"))
    }
  }
  
  symbols(x$grid$pts[,1], x$grid$pts[,2],
          circles=rep(.5, nrow(x$grid$pts)),
          inches=FALSE, add=TRUE, fg="black",
          bg=bgcolors)
  sdangles <-
    pi*(sdwccs - min(sdwccs, na.rm=TRUE)) / max(sdwccs, na.rm=TRUE)
  sdangles <- sdangles - 0.5*pi
  segments(x$grid$pts[,1], x$grid$pts[,2],
           x$grid$pts[,1] + .5*cos(sdangles),
           x$grid$pts[,2] + .5*sin(sdangles),
           col="blue")
  
  if (heatkey)
    plot.heatkey(x, zlim, unit.colors, labels = NULL, contin = TRUE, ...)
}


"plot.heatkey" <- function (x, zlim, unit.colors, labels, contin, ...)
{
  ncolors <- length(unit.colors)
  
  yrange <- range(x$grid$pts[, 2])
  smallestx <- min(x$grid$pts[,1])
  ## A width of .2 looks OK on my screen for small nets but is too
  ## small for large nets. Take the widest of either .2 units or
  ## 0.125 inches
  xleftin <- max(diff(par("usr")[1:2]) / (par("pin")[1] * 8), .2)
  xleft <- c(smallestx - 1 - xleftin, smallestx - 1)
  
  yleft <- seq(yrange[1] - 0.5,
               yrange[2] + 0.5,
               length = ncolors + 1)
  rect(xleft[1], yleft[1:ncolors],
       xleft[2], yleft[2:(ncolors + 1)],
       border = "black", col = unit.colors,
       xpd = TRUE)

  cex <- list(...)$cex

  if (contin) {
    zvals <- pretty(zlim)
    zvals <- zvals[zvals <= max(zlim) & zvals >= min(zlim)]
    yvals <- yrange[1] - .5 + (diff(yrange) + 1)*(zvals - zlim[1])/diff(zlim)

    text(xleft[2] - 1.4 * diff(xleft),
         yvals,
         formatC(zvals),
         xpd=TRUE, adj=1, cex=cex)
  } else {
    if (is.null(labels))
      labels <- 1:ncolors
    
    text(xleft[2] - 1.4 * diff(xleft),
         yleft[-1] - 0.5 * diff(yleft[1:2]),
         sort(labels),
         xpd = TRUE, adj=1, cex=cex)
  }
}
