repList <- function(x, times) {
  if(times == 1)
    l <- x
  else {
    l <- list()
    l <- lapply(1:times, function(i) x)
    names(l) <- paste("l", sapply(1:times, function(i) i), sep = "")
  }
  return(l)
}


.proportional_map <- function(z, maxz) {
  ## Proportional Symbol Mapping in R
  ## Susumu Tanimura, Chusi Kuroiwa, Tsutomu Mizota
  ## Journal of Sratistical Software January 2006
  sizes <- (abs(z) / maxz) ^ 0.57
  return(sizes) 
}


.symbol2pch <- function(symbol){
    ## give the pch associated to some symbol names (used in *.value)
    res <- 22 ## square by default
    if(symbol == "circle"){
        res <- 21
    }  else if(symbol == "diamond"){
        res <- 23
    } else if(symbol == "uptriangle"){
        res <- 24
    } else if(symbol == "downtriangle"){
        res <- 25
    }
    return(res)
}


.textpos <- function(xx, yy, origin = c(0, 0), n = length(xx)) {
  ## justification for labels and positions used in s.arrow and s.corcircle
   
  if(is.vector(origin) & length(origin) == 2) {
    xx <- xx - origin[1]
    yy <- yy - origin[2]
  } else
    stop("Invalid argument 'origin'")
  
  justif <- matrix(0, nrow = 2, ncol = n)
  for(i in 1:n)
    ## text justification
    ## move labels (w/2 for left/right or h/2 for bottom/top)
    if(yy[i] > 0) {
      if(abs(xx[i]) < yy[i])
        justif[, i] <- c(0, 1)
      else if(xx[i] < 0)
        justif[, i] <- c(-1, 0)
      else
        justif[, i] <- c(1, 0)
    } else { ## y<=0     
      if(abs(xx[i]) < abs(yy[i]))
        justif[, i] <- c(0, -1)
      else if(xx[i] < 0)
        justif[, i] <- c(-1, 0)
      else
        justif[, i] <- c(1, 0)
    }
  return(justif)
}


.setposition <- function(position) {
  ## specify where to draw grid text
  if(is.character(position)) {
    if(position == "bottomleft") {
      posi <- c(unit(0.05, "npc"), unit(0.02, "npc"))
      just <- c("left", "bottom")
    } else if(position == "bottomright") {
      posi <- c(unit(0.95, "npc"), unit(0.02, "npc"))
      just <- c("right", "bottom")
    } else if(position == "topleft") {
      posi <- c(unit(0.05, "npc"), unit(0.98, "npc"))
      just <- c("left", "top")
    } else if(position == "topright") {
      posi <- c(unit(0.95, "npc"), unit(0.98, "npc"))
      just <- c("right", "top")
    } else
      stop("Wrong position")
  } else {
    posi <- position
    just <- c("left", "bottom")
  }
  return(list(posi = posi, just = just))
}


.getgrid <- function(xlim, ylim, nbgrid = 5, origin, asp) {
  ## specify where to draw grid lines
  if(missing(origin)) {
    ## i.e. porigin.include = FALSE
    origin <- c(pretty(xlim, n = nbgrid)[1], pretty(ylim, n = nbgrid)[1])
  } 
  minX <- xlim[1]
  minY <- ylim[1]
  maxX <- xlim[2]
  maxY <- ylim[2]
  
  origin <- rep(origin, le = 2)
  cgrid.y <- diff(pretty(ylim, n = nbgrid))[1]
  cgrid.x <- diff(pretty(xlim, n = nbgrid))[1]
  
  if(asp == "iso") {
    if(diff(xlim) > diff(ylim))
      cgrid.x <- cgrid.y
    else
      cgrid.y <- cgrid.x
  }
  
  if(is.na(cgrid.x) || is.na(cgrid.y))
    stop("error while calculating grid")
  
  v0 <- origin[1]
  if((origin[1] + cgrid.x) <= maxX)
    v0 <- c(v0, seq(origin[1] + cgrid.x, maxX, by = cgrid.x))
  if((origin[1] - cgrid.x >= minX))
    v0 <- c(v0, seq(origin[1] - cgrid.x, minX, by = -cgrid.x))
  v0 <- sort(v0[v0 >= minX & v0 <= maxX])
  
  h0 <- origin[2]
  if((origin[2] + cgrid.y) <= maxY)
    h0 <- c(h0, seq(origin[2] + cgrid.y, maxY, by = cgrid.y))
  if((origin[2] - cgrid.y >= minY))
    h0 <- c(h0, seq(origin[2] - cgrid.y, minY, by = -cgrid.y))
  h0 <- sort(h0[h0 >= minY & h0 <= maxY])

  ## clean near-zero values
  delta <- diff(range(v0))/nbgrid
  if (any(small <- abs(v0) < 1e-14 * delta)) 
      v0[small] <- 0
  delta <- diff(range(h0))/nbgrid
  if (any(small <- abs(h0) < 1e-14 * delta)) 
      h0[small] <- 0
  
  res <- list(x0 = c(v0, rep(NA, length.out = length(h0))), x1 = c(v0, rep(NA, length.out = length(h0)))
    , y0 = c(rep(NA, length.out = length(v0)), h0), y1 = c(rep(NA, length.out = length(v0)), h0), d = signif(cgrid.x, 3))
  return(res)
}


setlimits1D <- function(mini, maxi, origin, includeOr) {
  ## computes limits for 1D plots
  if(includeOr) {
    newvalu <- .includeorigin(origin, mini, maxi)
    mini <- newvalu[1L]
    maxi <- newvalu[2L]
  }
  ## add 10% in both directions
  if(abs(diff(c(mini, maxi))) > .Machine$double.eps ^ 2)
    res <- c(mini, maxi) + c(-1, 1) * diff(c(mini, maxi)) / 10
  else { ## if there is only one value
    if(mini < .Machine$double.eps ^ 2)
  	  res <- mini + 0.02 * c(-1, 1)
    else
      res <- mini + c(-1, 1) * abs(mini) / 10 
  }
  return(res)
}


## if aspect.ratio == "iso", we must have identical limits range in x and y
setlimits2D <- function(minX, maxX, minY, maxY, origin = c(0, 0), aspect.ratio = "iso", includeOr) {
  origin <- rep(origin, length.out = 2)
  if(includeOr) { ## to include origin
    newvalu <- list(.includeorigin(origin[1], minX, maxX), .includeorigin(origin[2], minY, maxY))
    minX <- newvalu[[1L]][1L]
    minY <- newvalu[[2L]][1L]
    maxX <- newvalu[[1L]][2L]
    maxY <- newvalu[[2L]][2L]
  }
  
  ## interval sizes
  interX <- diff(c(minX, maxX))
  interY <- diff(c(minY, maxY))
  if(aspect.ratio == "iso") { ## same limits (to have iso square)
    biggest <- max(c(max(interX, interY)))
    if(which(c(interX, interY) == biggest)[1] == 1) { ## biggest is in X
      minY <- minY - (interX - interY) / 2
      maxY <- maxY + (interX - interY) / 2
    } else { ## biggest is in Y
      minX <- minX - (interY - interX) / 2
      maxX <- maxX + (interY - interX) / 2
    }
  }
  
  if(interX > .Machine$double.eps ^ 2 || interY > .Machine$double.eps ^ 2) {
  	xvalu <- c(minX, maxX) + c(-1, 1) * diff(c(minX, maxX)) / 10
  	yvalu <- c(minY, maxY) + c(-1, 1) * diff(c(minY, maxY)) / 10
  } else {
    xvalu <- c(minX, maxX) + c(-1, 1) * abs(max(minX, minY)) / 10
  	yvalu <- c(minY, maxY) + c(-1, 1) * abs(max(minX, minY)) / 10
  }
  
  return(list(xlim = xvalu, ylim = yvalu))
}


.includeorigin <- function(origin, value1, value2) {
  ## compute limits including origin
  return(range(c(origin, value1, value2)))
}


## separates a list of parameters to obtain 4 lists
## the first corresponds to 'padegraphic'
## the second to lattice parameters
## the third to graphics arguments
## the last to unused parameters
sortparamADEg <- function(...) {

  if(try(is.list(...), silent = TRUE) == TRUE) 
    dots <- as.list(...)
  else
    dots <- list(...)
  
  classtest <- try(list(...), silent = TRUE)
  if(class(classtest) == "try-error")
    stop("wrong parameters list, error in sortparamADEg")
  
  trellis <- list()
  adegpar <- list()
  g.args <- list()
  stats <- list()
  s.misc <- list()
  rest <- list()

  if(length(dots)) {
    
    ## compare to trellis parameters
    select <- separation(... , pattern = 1)
    trellis <- select[[1L]]
    rest <- select[[2L]]
    
    ## removing sp.layout items
    if(length(rest)) {
      indix2 <- pmatch(names(rest), "sp.layout")
      if(any(!is.na(indix2))) {
        whereis <- which(!is.na(indix2))
        g.args <- list("sp.layout" = rest[[whereis]])
        rest <- rest[-whereis]
      }
    }

    ## compare to adegpar parameters (pattern = 0 by default)
    if(length(rest)) {
      select <- separation(rest) 
      adegpar <- select[[1L]]
      rest <- select[[2L]]
    }
  
    ## removing g.args items
    if(length(rest)) {
      pattern.g.args <- c("xlim", "ylim", "main", "sub", "xlab", "ylab", "Sp", "nbobject", "samelimits", "scales", "key", "colorkey")
      pmatch.g.args <- pmatch(names(rest), pattern.g.args)
      indix <- which(!is.na(pmatch.g.args))
      pmatch.g.args <- pmatch.g.args[!is.na(pmatch.g.args)]
      if(length(indix)) {
        g.args <- c(g.args, rest[indix])
        names(g.args)[(1 + length(g.args) - length(pmatch.g.args)):length(g.args)] <- c(pattern.g.args[pmatch.g.args])
        rest <- rest[-indix]
      }
    }
  }
  return(list(adepar = adegpar, trellis = trellis, g.args = g.args, rest = rest))
}


########################################################################
###                         FROM MAPTOOLS                        #######
########################################################################

.pointLabel <- function(x, y = NULL, labels, width, height, limits, xyAspect, allowSmallOverlap = FALSE, trace = FALSE) {
  ## xyAspect: width_in/height_inch of the current panel
  ## limits would have been setted before (in ADEg.S2 prepare)
  ## width and height de rectangle en 'npc' (fig in original code__maptools package)
  ## labels <- graphicsAnnot(labels)
  ## to do before
  ## TODO redo
  boundary <- c(limits$xlim, limits$ylim)
  
  toUnityCoords <- function(xy) {
    return(list(x = (xy$x - boundary[1]) / (boundary[2] - boundary[1]) * xyAspect, y = (xy$y - boundary[3]) / (boundary[4] - boundary[3]) / xyAspect))
  }
  
  toUserCoords <- function(xy) {
    return(list(x = boundary[1] + xy$x / xyAspect * (boundary[2] - boundary[1]), y = boundary[3] + xy$y * xyAspect * (boundary[4] - boundary[3])))
  }
  
  z <- xy.coords(x, y, recycle = TRUE)
  z <- toUnityCoords(z)
  x <- z$x
  y <- z$y
  if(allowSmallOverlap) 
    nudgeFactor <- 0.02
  n_labels <- length(x)
  gen_offset <- function(code) {
    c(-1, -1, -1, 0, 0, 1, 1, 1)[code] * (width / 2) + (0 + 1i) * c(-1, 0, 1, -1, 1, -1, 0, 1)[code] * height / 2
  }
  
  rect_intersect <- function(xy1, offset1, xy2, offset2) {  ##intersections calculations
    w <- pmin(Re(xy1 + offset1 / 2), Re(xy2 + offset2 / 2)) - pmax(Re(xy1 - offset1 / 2), Re(xy2 - offset2 / 2))
    h <- pmin(Im(xy1 + offset1 / 2), Im(xy2 + offset2 / 2)) - pmax(Im(xy1 - offset1 / 2), Im(xy2 - offset2 / 2))
    w[w <= 0] <- 0
    h[h <= 0] <- 0
    w * h
  }
  
  nudge <- function(offset) {
    doesIntersect <- rect_intersect(xy[rectidx1] + offset[rectidx1], rectv[rectidx1], xy[rectidx2] + offset[rectidx2], rectv[rectidx2]) > 0
    pyth <- abs(xy[rectidx1] + offset[rectidx1] - xy[rectidx2] - offset[rectidx2]) / nudgeFactor
    eps <- 1e-10
    for (i in which(doesIntersect & pyth > eps)) {
      idx1 <- rectidx1[i]
      idx2 <- rectidx2[i]
      vect <- (xy[idx1] + offset[idx1] - xy[idx2] - offset[idx2]) / pyth[idx1]
      offset[idx1] <- offset[idx1] + vect
      offset[idx2] <- offset[idx2] - vect
    }
    offset
  }
  
  objective <- function(gene) { ## score calculations
    offset <- gen_offset(gene)
    if(allowSmallOverlap) 
      offset <- nudge(offset)
    if(!is.null(rectidx1)) 
      area <- sum(rect_intersect(xy[rectidx1] + offset[rectidx1], rectv[rectidx1], xy[rectidx2] + offset[rectidx2], rectv[rectidx2]))
    else 
      area <- 0
    n_outside <- sum(Re(xy + offset - rectv / 2) < 0 | Re(xy + offset + rectv / 2) > xyAspect | Im(xy + offset - rectv / 2) < 0 | Im(xy + offset + rectv / 2) > 1 / xyAspect)
    if(is.na(n_outside))
      n_outside <- 0 ## TODO: to correct, n_outside sometimes NA
    res <- 1000 * area + n_outside
    res
  }
  
  xy <- x + (0 + 1i) * y
  rectv <- width + (0 + 1i) * height
  rectidx1 <- rectidx2 <- array(0, (length(x)^2 - length(x)) / 2)
  k <- 0
  for(i in 1:length(x)) {
    for(j in seq(len = (i - 1))) {
      k <- k + 1
      rectidx1[k] <- i
      rectidx2[k] <- j
    }
  }
  
  canIntersect <- rect_intersect(xy[rectidx1], 2 * rectv[rectidx1], xy[rectidx2], 2 * rectv[rectidx2]) > 0
  rectidx1 <- rectidx1[canIntersect]  ## which intersect with those in rectidx2
  rectidx2 <- rectidx2[canIntersect]
  if(trace)
    cat("possible intersects = ", length(rectidx1), "\n")
  if(trace) 
    cat("portion covered = ", sum(rect_intersect(xy, rectv, xy, rectv)), "\n")
  ## simulated annealing
  SANN <- function() {
    ## initialisation
    gene <- rep(8, n_labels) ## 'rep' is best to begin at center
    score <- objective(gene) ## initial score 
    bestgene <- gene
    bestscore <- score
    T <- 2.5 ## pseudo initial temperature
    for (i in 1:50) {
      k <- 1
      for (j in 1:50) {
        newgene <- gene
        newgene[sample(1:n_labels, 1)] <- sample(1:8, 1)
        newscore <- objective(newgene)  ## score
        
        if(newscore <= score || runif(1) < exp((score - newscore) / T)) {
          ## empirical law to accept differences: if 'newscore' is better or with a proba exp(Dscorce/T)
          k <- k + 1
          score <- newscore
          gene <- newgene
        }
        if(score <= bestscore) {
          bestscore <- score
          bestgene <- gene
        }
        if(bestscore == 0 || k == 10) 
          break
      }
      if(bestscore == 0) ## no variation
        break
      if(trace)
        cat("overlap area =", bestscore, "\n")
      T <- 0.9 * T ## the temperature regularly decreases to become stable
    }
    
    if(trace) 
      cat("overlap area =", bestscore, "\n")
    nx <- Re(xy + gen_offset(bestgene))
    ny <- Im(xy + gen_offset(bestgene))
    return(list(x = nx, y = ny))
  }
  xy <- SANN()
  xy <- toUserCoords(xy)
  return(xy)
}



## check if z is included in breaks
## no default value
breakstest <- function(breaki, zi, n) {
  breaki <- sort(breaki, decreasing = TRUE)
  if(max(breaki) < max(zi) | min(breaki) > min(zi)) {
    zbis <- pretty(zi, n)
    if(max(breaki) < max(zi)) {
      warning(paste("breaks given does not include z limits,  break added ", max(zbis), sep = " "), call. = FALSE)
      breaki <- c((max(zbis)), breaki)
    }
    if(min(breaki) > min(zi)) {
      warning(paste("breaks given does not include z limits, break added ", min(zbis), sep = " "), call. = FALSE)
      breaki <- c(breaki, min(zbis))
    }
  }
  return(breaki)
}


################ for axis.....
## extract from
## Lattice Graphs { Control of Panel of Panel & Strip Borders
## J H Maindonald
## http://www.maths.anu.edu.au/~johnm
axis.L <- function(side, ..., line.col) {
  col <- trellis.par.get("axis.text")$col
  axis.default(side, ..., line.col = col)
}


.textsize <- function(labels, plabels) {
  ## can be improved see s1d.barchart for non-trivial rotation
  srt <- 0
  if(is.numeric(plabels$srt)) 
    srt <- plabels$srt[1]
  else {
    if(plabels$srt[1] == "horizontal") 
      srt <- 0
    else if(plabels$srt[1] == "vertical") 
      srt <- 90
  }
  
  if(srt == 90) { 
    h <- (convertHeight(stringWidth(labels), unitTo = "native", valueOnly = TRUE) + convertHeight(stringWidth("h"), unitTo = "native", valueOnly = TRUE) / 2) * rep(plabels$cex, length.out = length(labels))
    w <- (convertWidth(stringHeight(labels), unitTo = "native", valueOnly = TRUE) + convertWidth(stringHeight("m"), unitTo = "native", valueOnly = TRUE) / 2) * rep(plabels$cex, length.out = length(labels))
  } else { ## if 0 or an other angle
    w <- (convertWidth(stringWidth(labels), unitTo = "native", valueOnly = TRUE) + convertWidth(stringWidth("m"), unitTo = "native", valueOnly = TRUE) / 2) * rep(plabels$cex, length.out=length(labels))
    h <- (convertHeight(stringHeight(labels), unitTo = "native", valueOnly = TRUE) + convertHeight(stringHeight("h"), unitTo = "native", valueOnly = TRUE) / 2) * rep(plabels$cex, length.out=length(labels))
  }
  return(list(w = w, h = h, srt = srt))
}


.expand.call <- function(thecall, eval.formals = TRUE) {
  ## takes a call as argument and return a "cleaned" call where argument names are filled, and unset non empty formals are added and eventually evaluated using the call as environment
  ## supplied args:
  ans <- as.list(thecall)
  
  ## possible args:
  frmls <- formals(as.character(ans[[1]]))
  ## remove formal args with no presets:
  frmls <- frmls[!sapply(frmls, is.symbol)]
  add <- which(!(names(frmls) %in% names(ans)))
  frmls <- frmls[add]
  
  if(eval.formals) {
    ## evaluate the call locally and recursively
    frmls.new <- lapply(frmls, function(x) do.call("substitute", list(x, c(ans[-1], frmls))))
    while(!isTRUE(all.equal(frmls, frmls.new))) {
      frmls <- frmls.new
      frmls.new <- lapply(frmls, function(x) do.call("substitute", list(x, c(ans[-1], frmls))))
    }
  }
  return(c(ans, frmls))
}

