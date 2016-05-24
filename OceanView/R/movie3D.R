## =============================================================================
## moving slices in Open GL graphics
## =============================================================================

FindInterval <- function(x, vec, ...) {

  if (all(diff(vec) < 0)) { 
    vec <- rev(vec)
    res <- c(length(vec):1) [findInterval(x, vec, ...)]-1
  } else 
    res <- findInterval(x, vec, ...)
  res [ res == 0] <- 1
  res
}

movieslice3D <- function(x, y, z, colvar = NULL, xs = NULL,
  ys = NULL, zs = NULL, along = NULL, col = jet.col(100), NAcol = "white",
  breaks = NULL, colkey = FALSE, clim = NULL, clab = NULL,
  wait  = NULL, ask = FALSE, add = FALSE,
  basename = NULL, ...) {

  dots <- list(...)
  Lengths <- c(length(xs), length(ys), length(zs))

  if (sum(Lengths) == 0) {
    xs <- x
    Lengths[1] <- length(xs)
  }
  
  if (is.null(along))
    along <- which.max(Lengths)
  else if(Lengths[along] == 0)
    stop("'along' should refer to a dimension over which the slice is to move")

 if (is.null(clim))  {
    clim <- range(colvar, na.rm = TRUE)
  } else {
    colvar[colvar < min(clim)] <- NA
    colvar[colvar > max(clim)] <- NA
  }
  crange <- diff(clim)
  cmin <- min(clim)
  N      <- length(col) -1


#  if (! add)
  xsN <- xs
  ysN <- ys
  zsN <- zs
  if (along == 1)
    xsN <- NULL
  else if (along == 2)
    ysN <- NULL
  else if (along == 3)
    zsN <- NULL

  do.call("slice3Drgl", c(alist(x, y, z, colvar, xs = xsN, ys = ysN, zs = zsN, 
     add = add, clim = clim, col = col, NAcol = NAcol, breaks = breaks,
     colkey = colkey, clab = clab), dots))

  # the colors, 1.000..1 to avoid that trunc(1) = 0
  Col <- colvar
  if (is.null(breaks))  {
    Col[is.na(Col)] <- cmin
    Col[] <- col[1 + trunc((Col - cmin)/crange*1.00000000001*N)]
  } else {
    zi <- .bincode(colvar, breaks, TRUE, TRUE)
    Col[] <- col[zi]
  }
  Col[is.na(colvar)] <- NAcol

  popit <- FALSE
  if (!interactive())
    ask  <- FALSE
    
 # if main is passed...
  plist <- getplist()
  if (!is.null(dots$main))  {
    if (is.null(plist$dot$main))
      plist$dot$main <- dots$main
    ids <- plist$rgl$D
    if (! is.null(ids))
      if (!is.null(ids$main)) 
        rgl.pop(type = "shapes", id = ids["main"])
    M <- mtext3d(dots$main, "x++", line = 2)
    dots$main <- NULL
    plist$rgl$D$main <- M
  }

 # Function to add images on a plane as polygons
  image.plane <- function(xs, ys, zs, paint = FALSE, i = 0) {
    if (i == 1)   # all xs are same
      ix <-  FindInterval(xs[1], x, all.inside = TRUE)
    else      
      ix <- FindInterval(xs, x, all.inside = TRUE)
    if (i == 2) 
      iy <- FindInterval(ys[1], y, all.inside = TRUE)
    else 
      iy <- FindInterval(ys, y, all.inside = TRUE)
    if (i == 3) 
      iz <- FindInterval(zs[1], z, all.inside = TRUE)
    else 
      iz <- FindInterval(zs, z, all.inside = TRUE)
    
   # colorvar 
    cv <- matrix(nrow = nrow(xs), ncol = ncol(xs), data = Col[cbind(ix, iy, iz)])
    if (popit) 
      rgl.pop(id = 0)
    popit <<- TRUE
    persp3d(xs, ys, zs, col = cv, add = TRUE, #aspect = FALSE, 
            front = "filled", back = "filled")#, smooth = smooth, alpha = fb$alpha)
    
  } # end function imageplane

 # Function to first create a plane and then draw an image on it
  add.plane <- function(xs, ys, zs, i = 0) {

   if (ask) 
     readline("Hit enter to change slice")
   else if (! is.null(wait))
     Sys.sleep(wait)
    M <- mesh(xs, ys, zs)
    image.plane (M$x[,,], M$y[,,], M$z[,,], i = i) # [,,] to make sure it is an array
   if (! is.null(basename)) {
      filename <- paste(basename, formatC(as.integer(fi), width = 4, format = "d", flag = "0"),".png", sep="")
      rgl.snapshot(filename) 
      fi <<- fi + 1
   } 

  }  # end function addplane

  fi <- 0
  if (along == 1)
      for (x.s in xs[!is.na(xs)])
        add.plane(x.s, y, z, 1)
     
  else if (along == 2)
      for (y.s in ys[!is.na(ys)]) 
        add.plane(x, y.s, z, 2)
    
  else if (along == 3)
      for (z.s in zs[!is.na(zs)]) 
        add.plane(x, y, z.s, 3)
        
  if (popit) 
    rgl.pop(id = 0)

  setplist(plist)
}


moviepoints3D <- function(x, y, z, colvar, t, by = 1, 
  col = jet.col(100), NAcol = "white", breaks = NULL,
  clim = NULL, wait  = NULL, ask = FALSE,
  add = FALSE, basename = NULL, ...) {

    x <- as.vector(x)
    y <- as.vector(y)
    z <- as.vector(z)
    len <- length(x)
    if (length(y) != len) 
        stop("'y' should be of same length as 'x'")
    if (length(z) != len) 
        stop("'z' should be of same length as 'x'")
#    save <- par3d(skipRedraw = TRUE, ignoreExtent = TRUE)
#    on.exit(par3d(save))
    dots <- list(...)
    if (is.null(col)) 
        col <- "black"
    breaks <- check.breaks(breaks, col)

    if (!is.null(colvar)) {
        if (length(colvar) != len) 
            stop("dimension of 'colvar' should be equal to dimension of 'x', 'y' and 'z'")
        if (is.null(clim)) 
            clim <- range(colvar, na.rm = TRUE)
        Col <- variablecol(colvar, col, NAcol, clim, breaks)
    }
    else {
        Col <- col
    }

  if (!interactive()) 
    ask  <- FALSE
  tt <- unique(t)
  popit <- FALSE

  main <- paste("time ", tt)
  if (!is.null(dots$main))
    main <- rep(dots$main, length.out = length(tt))
  dots$main <- NULL
  
  for (i in seq(1, length(tt), by = by)) {
    T <- tt[i]
    if (ask) 
      readline("Hit enter to add new points")
    else if (! is.null(wait))
      Sys.sleep(wait)
    ii <- which(t == T)
    
    do.call ("addpoints", c(alist(x[ii], y[ii], z[ii], Col = Col[ii], popit, 
      main = main[i]), dots))
    if (! is.null(basename)) {
      filename <- paste(basename, formatC(as.integer(i), width = 4, format = "d", flag = "0"),".png", sep="")      
      rgl.snapshot(filename) 
    } 
    popit <- TRUE  
  }
}
  
addpoints <- function (x, y, z, Col, popit, main = NULL, ...) 
{
    dots <- list(...) 
    cex <- dots$cex
    pch <- dots$pch
    if (is.null(pch))
      pch <- "."
    if (popit) 
      rgl.pop(id = 0)
    if (is.null(cex)) 
        cex <- 1
    alpha <- dots$alpha
    plist <- getplist()
    if (is.null(alpha)) 
        alpha <- material3d()$alpha
    if (!is.null(main)) {
        ids <- plist$rgl$D$main
        if (!is.null(ids)) 
           rgl.pop(type = "shapes", id = ids)
        M <- mtext3d(main, "x++", line = 2)
        plist$rgl$D$main <- M
    }

    if (pch == ".")   
      size <- cex
    else
      size <- 6*cex  
    plot3d(x = x, y = y, z = z, size = size, col = Col, add = TRUE, 
        alpha = alpha)
    plist$pt <- list(x.mid = x, y.mid = y, z.mid = z, col = Col, 
        pch = rep(1, length(x)), bg = rep(1, length(x)), cex = rep(cex, 
            length(x)), alpha = rep(alpha, length(x)), proj = rep(NA, 
            length(x)))
    setplist(plist)
}
