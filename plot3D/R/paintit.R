check.breaks <- function (breaks, col) {  # adapted from image.default (graphics package)
  if (! is.null(breaks)) {
     nbreaks <- length(breaks)
     if (length(col) != nbreaks-1)
       stop("must have one more break than col - suggest to use jet.col(", nbreaks-1, ")")

  if (any(!is.finite(breaks)))
            stop("'breaks' must all be finite")
  if (is.unsorted(breaks)) {
            warning("unsorted 'breaks' will be sorted before use")
            breaks <- sort(breaks)
         }
  }
  return (breaks)
}

## =============================================================================
## Add image polygons, inputs are all matrices
## =============================================================================

paintit  <- function (colvar, x, y, z, plist, col, NAcol, clim,
                     border, facets, lwd, lty, dot, 
                     Extend = FALSE, Polar = FALSE, 
                     breaks = NULL) {

  dotshade <- dot$shade  
  if (! is.null(clim)) {
    if (length(clim) != 2)
        stop("'clim' should be a two-valued vector with the ranges of 'colvar'")
    colvar[colvar < min(clim)] <- NA
    colvar[colvar > max(clim)] <- NA
  }  
            
 # Check the plotting arguments x and y
  if (! is.matrix(x))
    stop("'x' should be a matrix")
  if (! is.matrix(y))
    stop("'y' should be a matrix")

 # adapt color palette and range for values = NA
  if (!ispresent(colvar)) {
    if (ispresent(col))
      Col <- col[1]  # take first color
    else
      Col <- "grey"
    if (Extend)
      Col <- rep(Col, length(x))  
    else  
      Col <- rep(Col, length(x[-1,-1]))  

  } else if (any (is.na(colvar)) & ! is.null(NAcol)& is.null(breaks)) {
    CC <- checkcolors(colvar, col, NAcol, clim)
    col <- CC$col
    colvar <- CC$colvar
    clim <- CC$lim  
  }
                             
  cmin   <- clim[1]
  crange <- diff(clim)
  N      <- length(col) -1

  # the colors, 1.000..1 to avoid that trunc(1) = 0  
  if (ispresent(colvar))
    if (is.null(breaks))
      Col <- col[1 + trunc((colvar - cmin)/crange*1.00000000001*N)]
    else {
      zi <- .bincode(colvar, breaks, TRUE, TRUE)
      Col <- col[zi]
      Col[is.na(Col)] <- NAcol
    }
  if (Extend) {
    x <- extend(x)
    y <- extend(y)
    z <- extend(z)
  }
  
  sl <- Sortlist(x, y, z, plist, Polar)

  if (dotshade$type != "none") 
    Col <- facetcols (x, y, z, Col, dotshade, Extend = FALSE)

  imgcol <- matrix(nrow = nrow(x) - 1, ncol = ncol(x) -1, data = Col)
  
  alpha <- dot$alpha; if (is.null(alpha)) alpha <- NA

  img <- list(list(x = x, y = y, z = z, col = imgcol,
    NAcol = NAcol, breaks = breaks, sl = sl, facets = facets, border = border,
    lwd = lwd, lty = lty, alpha = alpha, mapped = FALSE))  

  poly <- list(img = img)
  class(poly) <- "poly"
  invisible(poly)
}

## =============================================================================

mapimg <- function (plist) {
  img <- plist$img
  if (length(img) > 0) { 
    poly <- plist$poly
    for (i in 1:length(img)) {
 
      if (is.null(img[[i]]$mapped)) 
        img[[i]]$mapped <- TRUE
      if (!img[[i]]$mapped) {
        Poly <- with (img[[i]], polyfill(x, y, z, col[sl$list], NAcol, facets, border, sl,        
            lwd, lty, sl$Proj[sl$list], alpha = alpha))
        poly <- addPoly(poly, Poly)
        img[[i]]$mapped <- TRUE       
      }
    }
    plist$img <- img
    plist$poly <- poly
  }
  return(plist)
}

addPoly <- function (poly, Poly) {
  if (is.null(poly) | is.null(poly$x)) {
    poly <- Poly      
  } else if (! is.null(Poly)){
    if (!is.null(Poly$x)) {
      nR1 <- nrow(Poly$x) 
      nR2 <- nrow(poly$x)
      if (nR1 > nR2) {
        nR <- matrix(nrow = nR1 - nR2, ncol = ncol(poly$x), data = NA)
        poly$x <- rbind(poly$x, nR)
        poly$y <- rbind(poly$y, nR)
        poly$z <- rbind(poly$z, nR)
      } else if (nR2 > nR1) {
        nR <- matrix(nrow = nR2 - nR1, ncol = ncol(Poly$x), data = NA)
        Poly$x <- rbind(Poly$x, nR)
        Poly$y <- rbind(Poly$y, nR)
        Poly$z <- rbind(Poly$z, nR)
      }

      poly$x      <- cbind(poly$x,  Poly$x)
      poly$y      <- cbind(poly$y,  Poly$y)
      poly$z      <- cbind(poly$z,  Poly$z)
      poly$proj   <- c(poly$proj,   Poly$proj)
      poly$lwd    <- c(poly$lwd,    Poly$lwd)
      poly$lty    <- c(poly$lty,    Poly$lty)
      poly$border <- c(poly$border, Poly$border)
      poly$col    <- c(poly$col,    Poly$col)
      poly$alpha  <- c(poly$alpha,  Poly$alpha)
      poly$isimg  <- c(poly$isimg,  Poly$isimg)
    }
  }
  poly$img <- NULL
  return(poly)
}

## =============================================================================
## sort facets to draw according to view
## =============================================================================

sortlistvec <- function (x, y, z, plist, ignorez = TRUE) {                      
 
  Proj    <- project(x, y, z, plist, ignorez)
  sortlist <- sort.int(Proj, index.return = TRUE)$ix
  list(list = sortlist, Proj = Proj)
}

Sortlist <- function (x, y, z, plist, Polar = FALSE) {

  if (Polar)
    zz <- meangrid(z)
  else 
    zz <- 0

  xx <- meangrid(x)
  yy <- meangrid(y)

  sl <- sortlistvec(as.vector(xx), as.vector(yy), as.vector(zz), plist, !Polar)
  ind <- expand.sort(sl$list, dim(x)-1) 

  ix <- ind$x; iy <- ind$y

  if (Polar) {
    NN <- length(ix) * 0.5
    ix <- ix[- (1:NN)]
    iy <- iy[- (1:NN)]
    sl$list <- sl$list[- (1:NN)]
    maxProj <- sl$Proj[sl$list[length(sl$list)]]
    sl$Proj <- sl$Proj[sl$Proj <= maxProj] #sl$Proj[sl$list]
  }
  list(ix = ix, iy = iy, list = sl$list, 
    Proj = sl$Proj)
}


## =============================================================================
## Create polygons
## =============================================================================

createpoly <- function (x, y, z, ix, iy, Extend = TRUE) {
    
  if (Extend) {
    xx <- extend(x)
    yy <- extend(y)
    zz <- extend(z)
  } else {
    xx <- x
    yy <- y
    zz <- z
  }
 # the polygons
  PolyX <- rbind(xx[cbind(ix,     iy    )],
                 xx[cbind(ix + 1, iy    )],
                 xx[cbind(ix + 1, iy + 1)],
                 xx[cbind(ix,     iy + 1)], NA)
  PolyY <- rbind(yy[cbind(ix,     iy    )],
                 yy[cbind(ix + 1, iy    )],
                 yy[cbind(ix + 1, iy + 1)],
                 yy[cbind(ix,     iy + 1)], NA)
  PolyZ <- rbind(zz[cbind(ix,     iy    )],
                 zz[cbind(ix + 1, iy    )],
                 zz[cbind(ix + 1, iy + 1)],
                 zz[cbind(ix,     iy + 1)], NA)

  list(X = PolyX, Y = PolyY, Z = PolyZ, xx = xx, yy = yy, zz = zz)

}

## =============================================================================
## Draw polygons
## =============================================================================

polyfill <- function(x, y, z, Col, NAcol, facets, border, sl,
                     lwd, lty, proj = NULL, alpha = NA) {

  Poly <- createpoly(x, y, z, sl$ix, sl$iy, Extend = FALSE) 
  
  if (any (is.na(x) | is.na(y) | is.na(z))) {
    i1 <- which(is.na(Poly$X[-5, ]))
    i2 <- which(is.na(Poly$Y[-5, ]))
    i3 <- which(is.na(Poly$Z[-5, ]))
    ii <- unique(c(i1, i2, i3))
    Poly$X[-5, ][ii] <- NA
    Poly$Y[-5, ][ii] <- NA
    Poly$Z[-5, ][ii] <- NA

    ina <- apply (Poly$X[-5, ], MARGIN = 2, FUN = function(x) 
      any(is.na(x)) & !all(is.na(x)))
      
    for (i in (1:ncol(Poly$X)) [ina]){
      ii <- which(!is.na(Poly$X[1:4, i]))
      Poly$X[,i] <- c(Poly$X[ii,i], rep(NA, 5-length(ii)))
      Poly$Y[,i] <- c(Poly$Y[ii,i], rep(NA, 5-length(ii)))
      Poly$Z[,i] <- c(Poly$Z[ii,i], rep(NA, 5-length(ii)))
    }  

 # remove columns with only NAs or with all but one NA
    notNA  <- ! (is.na(Poly$X[2,]))
    Poly$X <- Poly$X[, notNA]
    Poly$Y <- Poly$Y[, notNA]
    Poly$Z <- Poly$Z[, notNA]
    
    if (length(Col) == length(notNA))
      Col <- Col[notNA]
    if (length(border) == length(notNA))
      border <- border[notNA]
    if (length(lwd) == length(notNA))
      lwd <- lwd[notNA]
    if (length(lty) == length(notNA))
      lty <- lty[notNA]
    proj <- proj[notNA]  
  } 

    
# The colors
  Col <- createcolors(facets, border, Col[])  

  if (is.null(lwd))
    lwd <- 1
  if (is.null(lty))
    lty <- 1
  if (is.null(alpha)) 
    alpha <- NA
    
 # update and return polygons.
  poly <- list(
       x      = Poly$X,
       y      = Poly$Y,
       z      = Poly$Z,                                  
       col    = Col$facet,
       border = Col$border,
       lwd    = rep(lwd, length.out = ncol(Poly$X)),
       lty    = rep(lty, length.out = ncol(Poly$X)),
       isimg  = rep(1, length.out = ncol(Poly$X)), 
       alpha  = rep(alpha, length.out = ncol(Poly$X)), 
       proj   = proj)
  class(poly) <- "poly"
  return(poly)
    
}

## =============================================================================
## Same for 2D plots
## =============================================================================

polyfill2D <- function(x, y, Col, facets, border, lwd, lty, Extend = TRUE) {

 # polygons are painted
  if (Extend) {
   nr <- nrow(x)
   nc <- ncol(x)
   x <- extend(x)
   y <- extend(y)
  } else {
   nr <- nrow(x)-1
   nc <- ncol(x)-1
  }
   ix <- rep(1:nr, nc)
   iy <- as.vector(matrix(nrow = nr, ncol = nc,
                   data = 1:nc, byrow =TRUE))

 # the polygons
  PolyX <- rbind(x[cbind(ix,     iy    )],
                 x[cbind(ix + 1, iy    )],
                 x[cbind(ix + 1, iy + 1)],
                 x[cbind(ix,     iy + 1)], NA)
  PolyY <- rbind(y[cbind(ix,     iy    )],
                 y[cbind(ix + 1, iy    )],
                 y[cbind(ix + 1, iy + 1)],
                 y[cbind(ix,     iy + 1)], NA)

# The colors
  Col <- createcolors(facets, border, Col)

  if (is.null(lwd))
    lwd <- 1
  if (is.null(lty))
    lty <- 1

 # update and return polygons.
  poly <- list(
       x      = PolyX,
       y      = PolyY,
       col    = Col$facet,
       border = Col$border,
       lwd    = rep(lwd , length.out = ncol(PolyX)),
       lty    = rep(lty , length.out = ncol(PolyX)))
  class(poly) <- "poly"
  return(poly)
}

