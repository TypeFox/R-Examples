
## =============================================================================
## RGL versions of the functions in package plot3D
## =============================================================================
par3dpars <- c("antialias","cex","family","font","useFreeType","fontname",
  "FOV","ignoreExtent","modelMatrix", "mouseMode", "projMatrix","bg",
  "scale","skipRedraw","userMatrix","viewport","zoom","bbox", "windowRect")

plotrgl <- function(lighting = FALSE, new = TRUE, add = FALSE, smooth = FALSE, 
   ...) {

  plist <- getplist() 

  if (sum(plist$setlim) > 0) {
    xlim <- plist$xlim
    ylim <- plist$ylim
    zlim <- plist$zlim
    plist <- clipplist(plist$xlim, plist$ylim, plist$zlim)
    plist$xlim <- xlim
    plist$ylim <- ylim
    plist$zlim <- zlim
  }
  if (plist$type =="2D") {
    plotrgl2D(plist, new = new, add = add, smooth = smooth, plot = FALSE, ...)
    pp <- getplist()
    pp$twoD <- plist$twoD
    pp$converted <- TRUE
    plist <- pp
  } else
    plist <- plotrglplist (plist, lighting, new = new, add = add, smooth = smooth, 
         update = TRUE, scale = TRUE, ...) 

  plist$rgl$userMatrix <- par3d("userMatrix")
  setplist(plist)
  dots <- list(...)
  if (is.null(dots$mouseMode)) 
    pan3d(3)           # from the help of rgl.setMouseCallbacks
  invisible(plist)  
}

clipplist <- function(xlim, ylim, zlim, plist = NULL){  
  SS <- function (x, y, z) {
    val <- rep(TRUE, length(x))
    if (! is.null(xlim)) 
      val[x < min(xlim) | x > max(xlim)] <- FALSE
    if (! is.null(ylim)) 
      val[y < min(ylim) | y > max(ylim)] <- FALSE
    if (! is.null(zlim)) 
      val[z < min(zlim) | z > max(zlim)] <- FALSE
    return(val)  
  }
  if (is.null(plist))
    plist <- getplist()
  plist <- selectplist(plist, SS)
  return(plist)
}

## =============================================================================
## Plot all items of a plist
## =============================================================================

plotrglplist <- function(plist, lighting = FALSE, new = TRUE, smooth = FALSE,
  add = FALSE, update = TRUE, scale = TRUE, ...) {

  if (add)
    new <- FALSE
  dots <- list(...)
  clipit <- FALSE
  if (! is.null(dots$xlim))  {
    xlim <- dots$xlim        
    if (xlim[1] > min(plist$xlim) | xlim[2] < max(plist$xlim))
      clipit <- TRUE
  } else 
    xlim <- plist$xlim
    
  if (! is.null(dots$ylim)) {
    ylim <- dots$ylim
    if (ylim[1] > min(plist$ylim) | ylim[2] < max(plist$ylim))
      clipit <- TRUE
  } else 
    ylim <- plist$ylim

  if (! is.null(dots$zlim)){
    zlim <- dots$zlim
    if (zlim[1] > min(plist$zlim) | zlim[2] < max(plist$zlim))
      clipit <- TRUE
    
  } else 
    zlim <- plist$zlim
  
  if (clipit)
    plist <- clipplist(xlim, ylim, zlim, plist)
     
  dots$xlim <- dots$ylim <- dots$zlim <- NULL
  
  materialnames <- names(formals(rgl.material))
  material <- dots[names(dots) %in% materialnames] 
  material$lit <- lighting
  if (new) 
      do.call("open3d", dots[names(dots) %in% par3dpars]) #[!names(dots) %in% materialnames])
  else {
    if (! add & is.null(plist$clearit))
      rgl.clear()
    else if (! is.null(plist$clearit))
      if (plist$clearit) rgl.clear()
    
    dotmat <- dots[names(dots) %in% par3dpars]  
    if (length(dotmat) > 0)
      do.call("material3d", dotmat) #dots[!names(dots) %in% materialnames])
  }
  if (! add & ! is.null(plist$colkey))
    colkey3D(plist$colkey[[1]]$par, plist$colkey[[1]]$col, plist$colkey[[1]]$clim,
     plist$colkey[[1]]$clab, plist$colkey[[1]]$clog,
     plist$colkey[[1]]$New, dots$alpha)

  if (length(material) > 0)
    do.call("material3d", material)

  if (length(plist) == 0)
    stop("nothing to draw")

  if (scale & !add) 
    aspect3d(plist$scalefac$x, plist$scalefac$y, plist$scalefac$z)

  # images
  Poly <- plist$poly

  if (is.null(plist$imgnr)) 
    plist$imgnr <- 0
  if (plist$imgnr > 0) {
   # make col the correct size for persp3d
    changedimg <- FALSE
    for (i in 1:plist$imgnr) {
      img <- plist$img[[i]]
      if (is.null(plist$img[[i]]$col.full)) {
       
        changedimg <- TRUE
        if (any(dim(img$col) - dim(img$x)) != 0) {
          i1 <- cbind(img$col, img$col[,ncol(img$col)])
          CC <- rbind(i1[1,], i1)
          plist$img[[i]]$col.full <- CC 
        } else
          plist$img[[i]]$col.full <- img$col 
      }  
    }
    if (changedimg) setplist(plist)
    # plot the images
    for (i in 1:plist$imgnr) {
      img <- plist$img[[i]]
      fb <- Createcolors(img$facets, img$border, img$col.full, img$alpha)
      if (!all(is.na(fb$facets))) 
          persp3d(img$x, img$y, img$z, col = fb$facets, add = TRUE, aspect = FALSE,
            front = "filled", back = "filled", smooth = smooth, alpha = fb$alpha)
      if (!all(is.na(fb$border))) 
        persp3d(img$x, img$y, img$z, col = fb$border, add = TRUE, aspect = FALSE,
          front = "lines", back = "lines", smooth = smooth)
    }
    
   # polygons that are not images
    if ( length(Poly) > 0) {
      p <- !Poly$isimg
      poly <- list(x = Poly$x[,p], y = Poly$y[,p], z = Poly$z[,p], 
        col = Poly$col[p], border = Poly$border[p],
        lwd = Poly$lwd[p], lty = Poly$lty[p], alpha = Poly$alpha[p])
    } else poly <- NULL    
        
  } else poly <- Poly


 # two types of polygons
  if (length(poly$x) > 0) {
    if (is.null(poly$alpha))
      poly$alpha <- rep(material3d()$alpha, length(poly$col))
    else  
      poly$alpha[is.na(poly$alpha)] <- material3d()$alpha  
    if (nrow(poly$x) > 5)
      stop ("cannot handle polygons with more than 4 nodes in plotrgl")

    ifill <- which(! is.na(poly$col))
    if (length(ifill) > 0)
      rglpoly(poly, ifill, "filled")
    iline <- which(is.na(poly$col))
    if (length(iline) > 0)
      rglpoly(poly, iline, "line")
    
  }
 # line segments
  segm <- plist$segm
  if (length(segm$x.from) > 0) {
    ilwd <- unique(segm$lwd) 
    for (i in 1:length(ilwd)) {
      ii <- which(segm$lwd == ilwd[i])
      alpha <- segm$alpha[ii]
      alpha[is.na(alpha)] <- material3d()$alpha
        segments3d(x = rbind(as.vector(segm$x.from)[ii], as.vector(segm$x.to)[ii]),
                 y = rbind(as.vector(segm$y.from)[ii], as.vector(segm$y.to)[ii]),
                 z = rbind(as.vector(segm$z.from)[ii], as.vector(segm$z.to)[ii]),
                 color = matrix (nrow = 2, byrow = TRUE, data = rep(segm$col[ii], 2)),
                 lwd = ilwd[i], lty = segm$lty[ii[1]], 
                 alpha =  matrix (nrow = 2, byrow = TRUE, data = rep(alpha, 2))
        )
    }           
  }
 # arrows 
  arr <- plist$arr
  if (length(arr$x.from) > 0) {
    ilwd <- unique(arr$lwd) 
    for (i in 1:length(ilwd)) {
      ii <- which(arr$lwd == ilwd[i])
      alpha <- arr$alpha[ii]
      alpha[is.na(alpha)] <- material3d()$alpha
      
      arrfun(as.vector(arr$x.from)[ii], as.vector(arr$y.from)[ii], as.vector(arr$z.from)[ii], 
              as.vector(arr$x.to)[ii], as.vector(arr$y.to)[ii], as.vector(arr$z.to)[ii],
              code = arr$code[ii], col = arr$col[ii], lty = arr$lty[ii], 
              lwd = ilwd[i], length = arr$length[ii], angle = arr$angle[ii],
              type = arr$type[ii], alpha = alpha, 
              sx = 1/plist$scale$x, sy = 1/plist$scale$y, sz = 1/plist$scale$z)      
    }
  }
  pts <- plist$pt 
  if (length(pts$x.mid) > 0) {
    ii <- which(pts$pch == ".") 
    if (length (ii) > 0) {
      unsize <- unique(pts$cex[ii])
      for (j in unsize) {
        ij <- ii[which(pts$cex[ii] == j)]
      alpha <- pts$alpha[ij]
      alpha[is.na(alpha)] <- material3d()$alpha
        plot3d(x = pts$x.mid[ij], y = pts$y.mid[ij], z = pts$z.mid[ij],
               size = 1 *j, col = pts$col[ij], alpha = alpha, add = TRUE)
      }
    }
    ii <- which(pts$pch != ".") 
    
    if (length (ii) > 0) {
      unsize <- unique(pts$cex[ii])
      for (j in unsize) {
        ij <- ii[which(pts$cex[ii] == j)]
        alpha <- pts$alpha[ij]
        alpha[is.na(alpha)] <- material3d()$alpha
        
        plot3d(x = pts$x.mid[ij], y = pts$y.mid[ij], z = pts$z.mid[ij],
               size = 6 *j, col = pts$col[ij], alpha = alpha, add = TRUE)
      }
    } 
  }
  
  pts <- plist$CIpt
  if (length(pts$x.mid) > 0) {
    ii <- which(pts$pch == ".") 
    if (length (ii) > 0) {
      unsize <- unique(pts$cex[ii])
      for (j in unsize) {
        ij <- ii[which(pts$cex[ii] == j)]
        alpha <- pts$alpha[ij]
        alpha[is.na(alpha)] <- material3d()$alpha
        
        plot3d(x = pts$x.mid[ij], y = pts$y.mid[ij], z = pts$z.mid[ij],
              size = 1 *j, col = pts$col[ij], alpha = alpha, add = TRUE)
      }
    }
    ii <- which(pts$pch != ".") 
    if (length (ii) > 0) {
      unsize <- unique(pts$cex[ii])
      for (j in unsize) {
        ij <- ii[which(pts$cex[ii] == j)]
        alpha <- pts$alpha[ij]
        alpha[is.na(alpha)] <- material3d()$alpha

        plot3d(x = pts$x.mid[ij], y = pts$y.mid[ij], z = pts$z.mid[ij],
             size = 6 *j, col = pts$col[ij], alpha = alpha, add = TRUE)
      }
    }
    nCImax <- max(pts$nCI)
    
    for (i in 1: nCImax) {             
      ii <- which(i <= pts$nCI)
      alpha <- pts$alpha[ii]
      alpha[is.na(alpha)] <- material3d()$alpha
      
       segments3d(x = rbind(pts$x.from[ii,i], pts$x.to[ii,i]),
                 y = rbind(pts$y.from[ii,i], pts$y.to[ii,i]),
                 z = rbind(pts$z.from[ii,i], pts$z.to[ii,i]),
                 color = matrix (nrow = 2, byrow = TRUE, data = rep(pts$CIpar$col[ii], 2)),
                 lwd = pts$CIpar$lwd[1], lty = pts$CIpar$lty[1], alpha = alpha)
    }           
  }
  
  labs <- plist$labels
  if (length(labs) > 0) {
      alpha <- labs$alpha 
      alpha[is.na(alpha)] <- material3d()$alpha
    i1 <- which (labs$adj == 0) 
    if (length(i1))
    text3d(x = labs$x[i1], y = labs$y[i1], z = labs$z[i1],
        color = labs$col[i1], texts = labs$labels[i1], font = labs$font[i1],
        cex = labs$cex[i1], adj = 0, alpha = alpha[i1])
    i1 <- which (labs$adj == 0.5) 
    if (length(i1))
    text3d(x = labs$x[i1], y = labs$y[i1], z = labs$z[i1],
        color = labs$col[i1], texts = labs$labels[i1], font = labs$font[i1],
        cex = labs$cex[i1], adj = 0.5, alpha = alpha[i1])
    i1 <- which (labs$adj == 1) 
    if (length(i1))
    text3d(x = labs$x[i1], y = labs$y[i1], z = labs$z[i1],
        color = labs$col[i1], texts = labs$labels[i1], font = labs$font[i1],
        cex = labs$cex[i1], adj = 1, alpha = alpha[i1])
  }
  D <- NULL
  if (plist$persp$drawbox & !add) {
    axes <- FALSE
    box <- TRUE

    if (! is.null(plist$dot$ticktype))
      if (plist$dot$ticktype == "detailed")
        axes <- TRUE
    if (! is.null(plist$dot$bty)){
     if (axes & plist$dot$bty %in% c("b", "b2"))
       box <- FALSE
     if (plist$dot$bty == "n")
       box <- FALSE
    }
    D <- decorate3d(xlim = xlim, ylim = ylim, zlim = zlim, 
	     xlab = plist$dot$xlab, ylab = plist$dot$ylab, zlab =  plist$dot$zlab, 
	     main = plist$dot$main, sub = plist$dot$sub, axes = axes, box = box)
	  if (!axes)
      box3d()   
  } else if (!add)
    D <- title3d(xlab = "", ylab = "", zlab = "", 
	     main = plist$dot$main, sub = plist$dot$sub)	

  pp <- getplist()
  pp$rgl$lighting <- lighting 
  pp$rgl$smooth <- smooth
  pp$rgl$alpha <- material3d()$alpha
	if (! is.null(D)) 
    pp$rgl$D <- as.list(D)
  setplist(pp)
  invisible(pp)
}


## =============================================================================
## function as from plot3D to see whether and which colors to use
## =============================================================================
Createcolors <- function(facets, border, Cols, Alpha) {

  isfacets <- facets
  if (is.null(facets)) isfacets <- NA
  if (is.null(border)) border <- NA
  if (is.null(Cols)) Cols <- NA

  isnaborder <- FALSE
  if (length(border) > 0)
    isnaborder <- is.na(border)

  if (is.na(isfacets)) {
    if (isnaborder)        
      border <- Cols
    Cols[] <- NA

  } else if (isfacets) {    # facets added
    if (isnaborder) {
      border <- Cols
      border[] <- NA 
    }
  } else {
    if (is.na(border))      # no facets
      border <- Cols
    Cols[] <- "white"
  }
  if (length(border ) == 1) {
    bb <- border
    border <- Cols
    border[] <- bb
  }
  alpha <- material3d()$alpha
  if (all(!is.na(Alpha)))
    alpha <- Alpha
  if (! all(is.na(Cols)))
    if (any (Cols == "transparent")) {
      alpha <- matrix(nrow = nrow(Cols), 
                      ncol = ncol(Cols), 
                      data = alpha)
      alpha[Cols == "transparent"] <- 0
    }
  list(border = border, facets = Cols, alpha = alpha)
}

## =============================================================================
## Plot polygons of a plist
## =============================================================================
  
rglpoly <- function(poly, il, front) {
   
  if (is.null(il)) {    # choose all
    i.Tri  <- which(is.na (poly$x[4, ]))
    i.Quad <- which(!is.na(poly$x[4, ]))
  
  } else {
    i.Tri  <- il[which(is.na (poly$x[4, il]))]
    i.Quad <- il[which(!is.na(poly$x[4, il]))]
  }

    
  plotpoly <- function(ipol, func, ir) {
    if (length(poly$alpha) !=  length(poly$lwd))
      poly$alpha <- rep(poly$alpha, length.out = length(poly$lwd)) # to overcome an error in older versions plot3D
    if (front == "filled")  {
      it <- ipol[poly$col[ipol] != "transparent" ]
      F <- func(x = poly$x[1:ir, it], y = poly$y[1:ir, it], z = poly$z[1:ir, it],
        col = matrix (nrow = ir, byrow = TRUE, data = rep(poly$col[it], ir)), 
        front = front, lwd = poly$lwd[it[1]], alpha = poly$alpha[it]) 
      ii <- which(poly$border[ipol] != poly$col[ipol])
      if (length(ii) > 0) {
        is <- ipol[ii]
        irr <- c(1:ir, 1, ir+1)
          lines3d(x = poly$x[irr, is], y = poly$y[irr, is], z = poly$z[irr, is],
          col = matrix (nrow = ir+2, byrow = TRUE, data = rep(poly$border[is], ir+2)), 
          lty = poly$lty[is[1]], lwd = poly$lwd[is[1]], alpha = poly$alpha[it])
      }
    } else  
       F <- func(x = poly$x[1:ir, ipol], y = poly$y[1:ir, ipol], z = poly$z[1:ir, ipol],
          col = matrix (nrow = ir, byrow = TRUE, data = rep(poly$border[ipol], ir)), 
          front = front, back = front, lwd = poly$lwd[ipol[1]], alpha = poly$alpha[ipol]) 
  }
  if (length(i.Tri) > 0)
    plotpoly(i.Tri, triangles3d, 3)
  if (length(i.Quad) > 0) 
    plotpoly(i.Quad, quads3d, 4)
}

## =============================================================================
## An arrow function, not perfect at all
## =============================================================================
  
hypoth <- function(a, b) sqrt(a^2 + b^2)

arrfun <- function (x.from, y.from, z.from, x.to, y.to, z.to, code, 
    col, lty, lwd, ..., length, angle = 20, type = "simple",
    sx, sy, sz, alpha = NA)  # scales       
{
#  if (is.null(angle)) angle <- 30
  length <- length/2
  if (any(length < 0 | ! is.finite(length   )))
    stop ("invalid arrow head length")

  if (any(!is.finite(angle)))
    stop ("invalid arrow head angle")

  if (any(is.na(code) | code < 0  | code > 3))
    stop ("invalid arrow head specification")

  N <- length(x.from)
  Code   <- rep(code, length.out = N)
  Length <- rep(length, length.out = N)
  Col    <- rep(col, length.out = N)
  
  Angle  <- rep(angle*pi/180, length.out = N)
  Type   <- rep(type, length.out = N)
  Alpha <- rep(alpha, length.out = N)
  Alpha[is.na(Alpha)] <- material3d()$alpha
  
  alpha <- matrix(nrow = 2, byrow = TRUE, data = rep(Alpha, 2))
  color <- matrix (nrow = 2, byrow = TRUE, data = rep(Col, 2))
  segments3d(x = rbind(x.from, x.to), y = rbind(y.from, y.to),
             z = rbind(z.from, z.to), color = color,
             lwd = lwd[1], lty = lty[1], alpha = alpha)   

## ? something with scale???
  
  eps <- 1e-8
    
  arrhead <- function (code) {
  
    i <- which (Code == code)
    if (length(i) == 0) return()
    length <- Length[i]
    col    <- rbind(Col[i], Col[i], Col[i], Col[i])
    alpha   <- rbind(Alpha[i], Alpha[i], Alpha[i], Alpha[i])
    angle  <- Angle[i]
    type   <- Type[i]

    if (code %in% c(1, 3)) {
    	x1 <- x.from[i]
    	y1 <- y.from[i]
  	  z1 <- z.from[i]
  	  x0 <- x.to[i]
  	  y0 <- y.to[i]
  	  z0 <- z.to[i]
    } else {
    	x1 <- x.to[i]
    	y1 <- y.to[i]
    	z1 <- z.to[i]
  	  x0 <- x.from[i]
  	  y0 <- y.from[i]
  	  z0 <- z.from[i]
    }
  	xc <- x0 - x1
   	yc <- y0 - y1
    zc <- z0 - z1

    r <- sqrt((yc/sy)^2 + (xc/sx)^2 + (zc/sz)^2)
    phi  <- atan2(yc/sy, xc/sx) 
    thet <- acos((zc/sz)/r)  

    is <- which (type == "cone")
    if (length(is) > 0) {
      N <- 30
      pseq <- seq(0, 2*pi, length.out = N+1)
# the long version, with a loop      
#      for (i in 1: N) {
#        a.xy.1 <- angle[is]*cos(pseq[i])
#        a.xy.2 <- angle[is]*cos(pseq[i+1])
#        a.z.1 <- angle[is]*sin(pseq[i])
#        a.z.2 <- angle[is]*sin(pseq[i+1])
       
#       	x <- rbind (x1[is] + length[is] * sx * cos(phi[is] + a.xy.1)*sin(thet[is] + a.z.1), x1[is], 
#                    x1[is] + length[is] * sx * cos(phi[is] + a.xy.2)*sin(thet[is] + a.z.2))
#       	y <- rbind (y1[is] + length[is] * sy * sin(phi[is] + a.xy.1)*sin(thet[is] + a.z.1), y1[is],  
#                    y1[is] + length[is] * sy * sin(phi[is] + a.xy.2)*sin(thet[is] + a.z.2))
#     	  z <- rbind (z1[is] + length[is] * sz * cos(thet[is] + a.z.1), z1[is], 
#     	              z1[is] + length[is] * sz * cos(thet[is] + a.z.2))
#       }

      Nx <- length(is)
      x <- as.vector(outer (x1[is], 1:N, FUN = function (x,i)
             x + length[is] * sx * cos(phi[is] + angle[is]*cos(pseq[i]))*
                                   sin(thet[is] + angle[is]*sin(pseq[i]))))
      x <- rbind(x, rep(x1[is], times = N), c(x[-(1:Nx)], x[1:Nx]))

      y <- as.vector(outer (y1[is], 1:N, FUN = function (y,i)
             y + length[is] * sy * sin(phi[is] + angle[is]*cos(pseq[i]))*
                                   sin(thet[is] + angle[is]*sin(pseq[i]))))
      y <- rbind(y, rep(y1[is], times = N), c(y[-(1:Nx)], y[1:Nx]))

      z <- as.vector(outer (z1[is], 1:N, FUN = function (z,i)
             z + length[is] * sz * cos(thet[is] + angle[is]*sin(pseq[i]))))
      z <- rbind(z, rep(z1[is], times = N), c(z[-(1:Nx)], z[1:Nx]))

      triangles3d(x = x, y = y, z = z, col = col[1:3, is], lty = lty[is[1]], 
           front = "filled", lwd = lwd[1], alpha = alpha[1:3, is])                    
    }
    ii <- which (type != "cone") 
    if (length(ii) == 0) return()

   	x <- rbind (x1 + length * sx * cos(phi+angle)*sin(thet), x1, 
                x1 + length * sx * cos(phi-angle)*sin(thet), NA)
   	y <- rbind (y1 + length * sy * sin(phi+angle)*sin(thet), y1,  
                y1 + length * sy * sin(phi-angle)*sin(thet), NA)
   	z <- rbind (z1 + length * sz * cos(thet), z1, 
     	          z1 + length * sz * cos(thet), NA)
    is <- which (type == "simple")
    if (length(is) > 0)  	        
     lines3d(x[,is], y[,is], z[,is], col = col[,is], lwd = lwd[1], lty = lty[is[1]], alpha = alpha[,is])                    

    is <- which (type %in% c("triangle", "curved"))
    if (length(is) > 0)  	        
      triangles3d(x = x[1:3,is], y = y[1:3,is], z = z[1:3,is], col = col[1:3,is], lty = lty[is[1]], 
        front = "filled", lwd = lwd[1], alpha = alpha[1:3, is])                    
  }
  arrhead(1)
  arrhead(3)
  Code[Code == 3] <- 2
  arrhead(2)
}

## =============================================================================


plotrgl2D <- function(plist, new, add, smooth, plot = FALSE, ...) {
  
  checkdots <-  function(pdots, dots, add) {
    if (! add) {
      if (! is.null(dots$xlim)) 
        pdots$xlim <- dots$xlim
      if (! is.null(dots$ylim)) 
        pdots$ylim <- dots$ylim
    }
    pdots$colkey <- pdots$clab <- NULL
    pdots$resfac <- pdots$theta <- pdots$add <- NULL
    pdots$rasterImage <- NULL 
    pdots$new <- new
    pdots
  }
  
  img2Dnr <- cont2Dnr <- scat2Dnr <- arr2Dnr <- segm2Dnr <- rect2Dnr <- poly2Dnr <- text2Dnr <- 0
  dots <- list(...)

  p <- plist$twoD
  for (i in 1:length(p$order)) {
    plt <- p$order[i]
    if (plt  == "image") {
      img2Dnr <- img2Dnr + 1
      Dots <- checkdots(p$img2D[[img2Dnr]], dots, add)
      do.call ("image2Drgl", c(alist(add = add, smooth = smooth), Dots))
    } else if (plt  == "contour") {
      cont2Dnr <- cont2Dnr + 1
      Dots <- checkdots(p$cont2D[[cont2Dnr]], dots, add)
      do.call ("contour2Drgl", c(alist(add = add), Dots))
    } else if (plt == "scatter") {
      scat2Dnr <- scat2Dnr + 1
      Dots <- checkdots(p$scat2D[[scat2Dnr]], dots, add)
      do.call ("scatter2Drgl", c(alist(add = add), Dots))
    } else if (plt == "text") {
      text2Dnr <- text2Dnr + 1
      Dots <- checkdots(p$text2D[[text2Dnr]], dots, add)
      do.call ("text2Drgl", c(alist(add = add), Dots))
    } else if (plt %in% c("arrows", "ArrType")) {
      arr2Dnr <- arr2Dnr + 1
      Dots <- checkdots(p$arr2D[[arr2Dnr]], dots, add)
      do.call ("arrows2Drgl", c(alist(add = add), Dots))
    } else if (plt == "segments") {
      segm2Dnr <- segm2Dnr + 1
      Dots <- checkdots(p$segm2D[[segm2Dnr]], dots, add)
      do.call ("segments2Drgl", c(alist(add = add), Dots))
    } else if (plt == "rect") {
      stop("rect not supported in rgl")
    } else if (plt == "polygon") {
      stop("polygons not supported in rgl")
    }
    add <- TRUE
    new <- FALSE
  }
  
  invisible(plist)
 
}

