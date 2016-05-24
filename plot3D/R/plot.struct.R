## =============================================================================
## Plotting the information in the plotting list
## =============================================================================

plotdev <- function(...)  {
  x <- plot.plist(getplist(), ...)
  invisible(x)    
}

## =============================================================================

plot.plist <- function(x, ...)  {

  if (length(x) == 0)
    stop("nothing to draw")

  dot <- list(...)

  if (x$type == "2D") {
    x <- plot2Dplist(x, ...)
    setplist(x)
    return(invisible(x))
  }
  
  if (length(x$img) > 0) { 
    x <- mapimg(x)
    setplist(x)
  }
  
  dot <- list(...)

  if (length(dot) != 0) {
    projectnames <- c("theta", "phi", "xlim", "ylim", "zlim", 
                      "scale", "expand", "r", "d")
    if (sum( projectnames %in% names(dot)) > 0)
      x <- proj3D (x, dot[projectnames])
  
    shadenames <- c("ltheta", "lphi", "shade", "lighting")
    if (sum( shadenames %in% names(dot)) > 0)
      if (! is.null(x$poly))
        x$poly <- color3D(x$poly, x$scalefac, dot[shadenames], dot$lighting, dot[["alpha"]])
  
    if ("alpha" %in% names(dot))
      x <- alpha3D(x, dot[["alpha"]])
  }     
  
  x$persp$box <- FALSE
  x <- plot.struct.3D (x, plot = TRUE)

  if (x$type == "23D") {
    x <- plot2Dplist(x, ...)
    setplist(x)
  }

  invisible(x)  
}

## =============================================================================

plot.struct.3D <- function(plist, pt = NULL, CIpt = NULL, poly = NULL, 
  segm = NULL, labels = NULL, arr = NULL, other = NULL, plot = TRUE) {       

  if (plot) {
    if (is.null(plist$plt)) 
      plist$plt$ori   <- par("plt")

    if (!is.null(plist)) 
      par(plt = plist$plt$main)

    if (!plist$persp$box) 
      plist <- plotbox(plist)
  }
  
  plist <- update.3D(plist, pt, CIpt, poly, segm, labels, arr, other, 
    expand = plot)

  if (plot) {
#    if (length(plist$img) > 0) { 
#      plist <- mapimg(plist)
#      setplist(plist)
#    }
  
    plotlist3D(plist)
    if (plist$persp$bty == "f")
      drawfullbox(plist)
    if (! is.null(plist$colkey))
      drawallcols(plist)

    par(plt = plist$plt$ori)
    par(mar = par("mar"))
  }
  invisible(plist)  
}

## =============================================================================

alpha3D <- function(plist, alpha) { # makes colors transparant

  if (alpha == 1 | alpha == 0) 
    return (plist)
  
  plist$poly$col       <- setalpha(plist$poly$col, alpha)
  plist$poly$border    <- setalpha(plist$poly$border, alpha)
  plist$pt$col         <- setalpha(plist$pt$col, alpha)
  plist$pt$bg          <- setalpha(plist$pt$bg, alpha)
  plist$CIpt$col       <- setalpha(plist$CIpt$col, alpha)
  plist$CIpt$bg        <- setalpha(plist$CIpt$bg, alpha)
  plist$CIpt$CIpar$col <- setalpha(plist$CIpt$CIpar$col, alpha)
  plist$labels$col     <- setalpha(plist$labels$col, alpha)
  plist$segm$col       <- setalpha(plist$segm$col, alpha)
  plist$arr$col        <- setalpha(plist$arr$col, alpha)
  
  if (! is.null(plist$numkeys))
    for (i in 1:plist$numkeys)
      plist$colkey[[i]]$col  <- setalpha(plist$colkey[[i]]$col, alpha)
      
  return(plist)
}

## =============================================================================

color3D <- function(poly, scalefac, shade, lighting, hasalpha = NULL) {  # lighting and shading

  shade$lighting <- NULL
  
  shade <- check.shade(shade, lighting)
  
  if (is.null(hasalpha)) 
   if (any(!is.na(poly$alpha))) {
    shade$alpha <- poly$alpha
    shade$alpha[is.na(shade$alpha)] <- 1
  }
   
  light   <- setuplight(shade$lphi, shade$ltheta) [1:3]              
  px <- poly$x * scalefac$x
  py <- poly$y * scalefac$y
  pz <- poly$z * scalefac$z
  isthree <- which(is.na(px[4,]))
    
  if (length(isthree) == ncol(px))
    Normals <- normal.points.tri(rbind(px[1,], py[1,], pz[1,]) , 
                             rbind(px[2,], py[2,], pz[2,]) , 
                             rbind(px[3,], py[3,], pz[3,]) 
                           )
  
  else if (length(isthree) == 0)
    Normals <- normal.points(rbind(px[1,], py[1,], pz[1,]) , 
                             rbind(px[3,], py[3,], pz[3,]) , 
                             rbind(px[4,], py[4,], pz[4,]) ,
                             rbind(px[2,], py[2,], pz[2,]) 
                             )
  else {
    Normals <- list(u = rep(0, ncol(px)), 
                    v = rep(0, ncol(px)), w = rep(0, ncol(px)))
 
    NN <- normal.points.tri(rbind(px[1,isthree], py[1,isthree], pz[1,isthree]) , 
                             rbind(px[2,isthree], py[2,isthree], pz[2,isthree]) , 
                             rbind(px[3,isthree], py[3,isthree], pz[3,isthree]) 
                       )
    Normals$u[isthree] <- NN$u
    Normals$v[isthree] <- NN$v
    Normals$w[isthree] <- NN$w
    isfour <- which(!is.na(px[4,]))
        
    NN <- normal.points(rbind(px[1,isfour], py[1,isfour], pz[1,isfour]) , 
                             rbind(px[3,isfour], py[3,isfour], pz[3,isfour]) , 
                             rbind(px[4,isfour], py[4,isfour], pz[4,isfour]) ,
                             rbind(px[2,isfour], py[2,isfour], pz[2,isfour]) 
                             )
    Normals$u[isfour] <- NN$u
    Normals$v[isfour] <- NN$v
    Normals$w[isfour] <- NN$w
  
  }                           
    
  ii <- which (! is.na(poly$col) & poly$col != "transparent" )
  if (length(ii) > 0) {   
    pcol <- facetcols.shadelight(light, Normals, poly$col, shade)
    poly$col[ii] <- pcol[ii] 
  }    

  ii <- which (! is.na(poly$border))
  if (length(ii) > 0) {   
    pcol <-  facetcols.shadelight(light, Normals, poly$border, shade)
    poly$border[ii] <- pcol[ii] 
  }    
  
  poly
}

## =============================================================================

proj3D <- function(plist, dot) {

  if (!is.null(dot$theta))
    plist$persp$theta <- dot$theta

  if (!is.null(dot$phi))
    plist$persp$phi <- dot$phi

  lims <- c(!is.null(dot$xlim), !is.null(dot$ylim), !is.null(dot$zlim))
  if (sum(lims) > 0) {
    SS <- function(x, y, z) {
     sel <- rep(TRUE, length.out = length(x))
     if (lims[1])
       sel[x < min(dot$xlim) | x > max(dot$xlim)] <- FALSE
     if (lims[2])
       sel[y < min(dot$ylim) | y > max(dot$ylim)] <- FALSE
     if (lims[3])
       sel[z < min(dot$zlim) | z > max(dot$zlim)] <- FALSE
     return(sel)
    } 
    plist <- selectplist(plist, SS)
  } 
      
  if (!is.null(dot$xlim))
    plist$xlim <- dot$xlim

  if (!is.null(dot$ylim))
    plist$ylim <- dot$ylim

  if (!is.null(dot$zlim))
    plist$zlim <- dot$zlim

  if (!is.null(dot$scale))
    plist$dot$scale <- dot$scale

  if (!is.null(dot$expand))
    plist$dot$expand <- dot$expand

  lim <- setlim (plist$xlim, plist$ylim, plist$zlim, 
    plist$dot$scale, plist$dot$expand) 

  plist$scalefac <- lim

  if (!is.null(dot$r))
    plist$dot$r <- dot$r

  if (!is.null(dot$d))
    plist$dot$d <- dot$d

  plist$mat <- transmat (plist$persp$phi, plist$persp$theta, plist$scalefac, 
    plist$dot$r, plist$dot$d)

 # update projections
  if (!is.null(plist$pt)) 
    plist$pt$proj <- project(plist$pt$x.mid,  plist$pt$y.mid,  plist$pt$z.mid,  plist)

  if (!is.null(plist$CIpt))
    plist$CIpt$proj <- project(plist$CIpt$x.mid,  plist$CIpt$y.mid,  plist$CIpt$z.mid,  plist)

  if (!is.null(plist$labels))
    plist$labels$proj <- project(plist$labels$x, plist$labels$y, plist$labels$z, plist)

  if (!is.null(plist$poly))
    plist$poly$proj <- project(colMeans(plist$poly$x, na.rm = TRUE),
                               colMeans(plist$poly$y, na.rm = TRUE), 
                               colMeans(plist$poly$z, na.rm = TRUE), plist)

  if (!is.null(plist$segm)) 
    plist$segm$proj <- project(0.5*(plist$segm$x.from + plist$segm$x.to), 
                               0.5*(plist$segm$y.from + plist$segm$y.to), 
                               0.5*(plist$segm$z.from + plist$segm$z.to), plist)

  if (!is.null(plist$arr))
    plist$arr$proj <- project(0.5*(plist$arr$x.from + plist$arr$x.to), 
                              0.5*(plist$arr$y.from + plist$arr$y.to), 
                              0.5*(plist$arr$z.from + plist$arr$z.to), plist)
  return(plist)
} 

## =============================================================================
## updates plist
## =============================================================================

update.3D <- function(plist, pt = NULL, CIpt = NULL, poly = NULL, 
       segm = NULL, labels = NULL, arr = NULL, other = NULL, expand = TRUE) {       
# ------------------------------------------------------------------------------
# update plist with new elements     
# ------------------------------------------------------------------------------
  if (any(plist$setlim)) {
    PL <- list(pt = pt, CIpt = CIpt, poly = poly, segm = segm, labels = labels,
    arr = arr, imgnr = plist$imgnr)
  }
  
  if (is.null(plist$pt))
    plist$pt <- pt

  else if (! is.null(pt)) {
    plist$pt$x.mid  <- c(plist$pt$x.mid,  pt$x.mid)
    plist$pt$y.mid  <- c(plist$pt$y.mid,  pt$y.mid)
    plist$pt$z.mid  <- c(plist$pt$z.mid,  pt$z.mid)
    plist$pt$col    <- c(plist$pt$col,    pt$col)
    plist$pt$pch    <- c(plist$pt$pch,    pt$pch)
    plist$pt$bg     <- c(plist$pt$bg,     pt$bg)
    plist$pt$cex    <- c(plist$pt$cex,    pt$cex)   
    plist$pt$lwd    <- c(plist$pt$lwd,    pt$lwd)
    plist$pt$alpha  <- c(plist$pt$alpha,  pt$alpha)
    plist$pt$proj   <- c(plist$pt$proj,   pt$proj)
  }

  if (is.null(plist$CIpt))
    plist$CIpt <- CIpt
  else if (! is.null(CIpt)) {
    plist$CIpt$x.to   <- rbind(plist$CIpt$x.to,  CIpt$x.to)
    plist$CIpt$y.to   <- rbind(plist$CIpt$y.to,  CIpt$y.to)
    plist$CIpt$z.to   <- rbind(plist$CIpt$z.to,  CIpt$z.to)
    plist$CIpt$x.from <- rbind(plist$CIpt$x.from,CIpt$x.from)
    plist$CIpt$y.from <- rbind(plist$CIpt$y.from,CIpt$y.from)
    plist$CIpt$z.from <- rbind(plist$CIpt$z.from,CIpt$z.from)
    plist$CIpt$nCI    <- c(plist$CIpt$nCI,   CIpt$nCI)
    plist$CIpt$x.mid  <- c(plist$CIpt$x.mid, CIpt$x.mid)
    plist$CIpt$y.mid  <- c(plist$CIpt$y.mid, CIpt$y.mid)
    plist$CIpt$z.mid  <- c(plist$CIpt$z.mid, CIpt$z.mid)
    plist$CIpt$len    <- c(plist$CIpt$len,   CIpt$len)
    plist$CIpt$lty    <- c(plist$CIpt$lty,   CIpt$lty)
    plist$CIpt$lwd    <- c(plist$CIpt$lwd,   CIpt$lwd)
    plist$CIpt$col    <- c(plist$CIpt$col,   CIpt$col)
    plist$CIpt$pch    <- c(plist$CIpt$pch,   CIpt$pch)   
    plist$CIpt$bg     <- c(plist$CIpt$bg,    CIpt$bg)
    plist$CIpt$cex    <- c(plist$CIpt$cex,   CIpt$cex)   
    plist$CIpt$proj   <- c(plist$CIpt$proj,  CIpt$proj)
    plist$CIpt$alpha  <- c(plist$CIpt$alpha, CIpt$alpha)

    plist$CIpt$CIpar$col  <- c(plist$CIpt$CIpar$col, CIpt$CIpar$col)
    plist$CIpt$CIpar$lwd  <- c(plist$CIpt$CIpar$lwd, CIpt$CIpar$lwd)
    plist$CIpt$CIpar$lty  <- c(plist$CIpt$CIpar$lty, CIpt$CIpar$lty)
    plist$CIpt$CIpar$alen <- c(plist$CIpt$CIpar$alen, CIpt$CIpar$alen)
  }

  if (is.null(plist$labels))
    plist$labels <- labels
  else if (! is.null(labels)) {
    plist$labels$x      <- c(plist$labels$x,  labels$x)
    plist$labels$y      <- c(plist$labels$y,  labels$y)
    plist$labels$z      <- c(plist$labels$z,  labels$z)
    plist$labels$labels <- c(plist$labels$labels, labels$labels)
    plist$labels$adj    <- c(plist$labels$adj,    labels$adj)
    plist$labels$cex    <- c(plist$labels$cex,    labels$cex)
    plist$labels$col    <- c(plist$labels$col,    labels$col)
    plist$labels$font   <- c(plist$labels$font,   labels$font)
    plist$labels$srt    <- c(plist$labels$srt,    labels$srt)
    plist$labels$proj   <- c(plist$labels$proj,   labels$proj)
    plist$labels$alpha  <- c(plist$labels$alpha,  labels$alpha)
  }
    
  if (! is.null(poly)) {

      if (is.null(plist$imgnr)) {
        plist$imgnr <- 0
        plist$img <- list()
      }
      if (length(poly$img) > 0) {
        Img <- poly$img  
        for (ii in 1: length(Img)) {
          img <- Img[[ii]]
          plist$imgnr <- plist$imgnr + 1
          if (expand) {
            if (is.null(img$mapped)) 
              img$mapped <- TRUE
            if (!img$mapped) {
              Poly <- with (img, polyfill(x, y, z, col[sl$list], NAcol,        
                facets, border, sl, lwd, lty, sl$Proj[sl$list], alpha = alpha))
              poly <- addPoly(poly, Poly)
            }
          }
          plist$img[[plist$imgnr]] <- img
        }
      }
#    if (expand) 
      plist$poly <- addPoly(plist$poly, poly)
  }
  if (expand & length(plist$img) >0 )  
    plist <- mapimg(plist)

  if (is.null(plist$segm))
    plist$segm <- segm
  else if (! is.null(segm)) {
    plist$segm$x.from <- c(plist$segm$x.from, segm$x.from)
    plist$segm$y.from <- c(plist$segm$y.from, segm$y.from)
    plist$segm$z.from <- c(plist$segm$z.from, segm$z.from)
    plist$segm$x.to   <- c(plist$segm$x.to,   segm$x.to)
    plist$segm$y.to   <- c(plist$segm$y.to,   segm$y.to)
    plist$segm$z.to   <- c(plist$segm$z.to,   segm$z.to)
    plist$segm$proj   <- c(plist$segm$proj,   segm$proj)
    plist$segm$lwd    <- c(plist$segm$lwd,    segm$lwd)
    plist$segm$lty    <- c(plist$segm$lty,    segm$lty)
    plist$segm$col    <- c(plist$segm$col,    segm$col)
    plist$segm$alpha  <- c(plist$segm$alpha,  segm$alpha)
  }

  if (is.null(plist$arr))
    plist$arr <- arr
  else if (! is.null(arr)) {
    plist$arr$x.from <- c(plist$arr$x.from, arr$x.from)
    plist$arr$y.from <- c(plist$arr$y.from, arr$y.from)
    plist$arr$z.from <- c(plist$arr$z.from, arr$z.from)
    plist$arr$x.to   <- c(plist$arr$x.to,   arr$x.to)
    plist$arr$y.to   <- c(plist$arr$y.to,   arr$y.to)
    plist$arr$z.to   <- c(plist$arr$z.to,   arr$z.to)
    plist$arr$proj   <- c(plist$arr$proj,   arr$proj)
    plist$arr$code   <- c(plist$arr$code,   arr$code)
    plist$arr$angle  <- c(plist$arr$angle,  arr$angle)
    plist$arr$length <- c(plist$arr$length, arr$length)
    plist$arr$lwd    <- c(plist$arr$lwd,    arr$lwd)
    plist$arr$lty    <- c(plist$arr$lty,    arr$lty)
    plist$arr$col    <- c(plist$arr$col,    arr$col)
    plist$arr$type   <- c(plist$arr$type,   arr$type)
    plist$arr$alpha  <- c(plist$arr$alpha,  arr$alpha)
  }                                                                     
  return (plist)
}

## =============================================================================
## plots points, polygons, segments, labels, arrows
## =============================================================================

plotlist3D <- function(plist) {       

    pt     <- plist$pt
    CIpt   <- plist$CIpt    
    poly   <- plist$poly
    segm   <- plist$segm
    labels <- plist$labels
    arr    <- plist$arr
    other  <- plist$other

# ------------------------------------------------------------------------------
# project (x, y, z) to plane (x, y) and count number of different structures
# ------------------------------------------------------------------------------
    numStruct <- 0  
  
    if (! is.null(poly)) {     
      pol <- trans3D(x = poly$x, y = poly$y, z = poly$z, pmat = plist$mat)
      numStruct <- numStruct + 1          
    }
    if (! is.null(other))  {    
      stop ("structure 'other' undefined")
      numStruct <- numStruct + 1
    }
    if (! is.null(CIpt$x.from)) {   
      CI.to   <- trans3D(x = CIpt$x.to,   y = CIpt$y.to,   z = CIpt$z.to,   pmat = plist$mat)
      CI.from <- trans3D(x = CIpt$x.from, y = CIpt$y.from, z = CIpt$z.from, pmat = plist$mat)
      CI.mid  <- trans3D(x = CIpt$x.mid , y = CIpt$y.mid , z = CIpt$z.mid , pmat = plist$mat)
      numStruct <- numStruct + 1
    }
    if (! is.null(pt)) {
      pt.mid  <- trans3D(x = pt$x.mid , y = pt$y.mid , z = pt$z.mid , pmat = plist$mat)
      numStruct <- numStruct + 1
    }
    if (! is.null(segm)) {   
      segm.to   <- trans3D(x = segm$x.to,   y = segm$y.to,   z = segm$z.to,   pmat = plist$mat)
      segm.from <- trans3D(x = segm$x.from, y = segm$y.from, z = segm$z.from, pmat = plist$mat)
      numStruct <- numStruct + 1
    }
    if (! is.null(arr)) {   
      arr.to   <- trans3D(x = arr$x.to,   y = arr$y.to,   z = arr$z.to,   pmat = plist$mat)
      arr.from <- trans3D(x = arr$x.from, y = arr$y.from, z = arr$z.from, pmat = plist$mat)
      numStruct <- numStruct + 1
    }
    if (! is.null(labels)) {   
      lab <- trans3D(x = labels$x, y = labels$y, z = labels$z, pmat = plist$mat)
      numStruct <- numStruct + 1
    }    

    List <- c(pt$proj, CIpt$proj, poly$proj, segm$proj, 
                           labels$proj, arr$proj, other$proj)
    if (length(List) == 0) return()
    sortlist <- sort.int(List, index.return = TRUE)$ix

# ------------------------------------------------------------------------------
# check if there is only one types and plot if true
# ------------------------------------------------------------------------------
    if (numStruct == 1) {

      if (! is.null(pt)) {  
        points(pt.mid$x[sortlist], pt.mid$y[sortlist], 
               col = pt$col[sortlist], 
               pch = pt$pch[sortlist],
               cex = pt$cex[sortlist], 
               lwd = pt$lwd[sortlist],
               bg  = pt$bg[sortlist])

      } else if (!is.null(poly)) {  # only polygons
        polygon(pol$x[ ,sortlist], pol$y[ ,sortlist], 
                lwd = poly$lwd[sortlist],
                lty = poly$lty[sortlist], 
                border = poly$border[sortlist], 
                col = poly$col [sortlist])

      } else if (!is.null(segm)) {  # only segments
        segments(segm.from$x[sortlist], segm.from$y[sortlist],  
                 segm.to$x  [sortlist], segm.to$y[sortlist], 
                 col = segm$col[sortlist],
                 lwd = segm$lwd[sortlist],
                 lty = segm$lty[sortlist])

      } else if (!is.null(labels)) {    # only labels
        text(x = lab$x[sortlist], y = lab$y[sortlist], 
             labels = labels$labels[sortlist], 
             col = labels$col[sortlist], 
             adj = labels$adj[sortlist], 
             cex = labels$cex[sortlist],
             font = labels$font[sortlist],
             srt = labels$srt[sortlist[1]])

      } else if (!is.null(arr)) {  # only simple arrows
         ArrType (arr.from$x[sortlist], arr.from$y[sortlist], 
            arr.to$x[sortlist], arr.to$y[sortlist], arr$length[sortlist], 
            arr$angle[sortlist], arr$code[sortlist], 
            col = arr$col[sortlist], type = arr$type[sortlist], 
            lwd = arr$lwd[sortlist], lty = arr$lty[sortlist])
                 
      }
      if (is.null(CIpt))
        return()    
    }
 
   # it is a mix of types  
    Lpt     <- length(pt$proj)
    LCIpt   <- length(CIpt$proj)
    Lpoly   <- length(poly$proj)
    Lsegm   <- length(segm$proj)
    Llabels <- length(labels$proj)
    Larr    <- length(arr$proj)
    Lother  <- length(other$proj)
    LCI     <- Lpt + LCIpt
    LCP     <- LCI + Lpoly
    LCPS    <- LCP + Lsegm
    LCPSL   <- LCPS + Llabels
    LCPSLA  <- LCPSL + Larr
    Ltot    <- LCPSLA + Lother

    type <- rep(0, Ltot)                                   # 0 = points
    type[sortlist > Lpt   & sortlist <= LCI]    <- 1       # points + CI
    type[sortlist > LCI   & sortlist <= LCP]    <- 2       # polygons
    type[sortlist > LCP   & sortlist <= LCPS]   <- 3       # segments
    type[sortlist > LCPS  & sortlist <= LCPSL]  <- 4       # labels
    type[sortlist > LCPSL & sortlist <= LCPSLA] <- 5       # arrows
    type[sortlist > LCPSLA                    ] <- 6       # other - undefined

    plotit <- function(ii) {
      i <- sortlist[ii]  
    
      if (type[ii] == 0) { # points   
        points(pt.mid$x[i], pt.mid$y[i], 
               col = pt$col[i], pch = pt$pch[i],
               cex = pt$cex[i], bg = pt$bg[i])
   
     } else if (type[ii] == 1) { # points + CI   
        io  <- i - Lpt
        nCI <- CIpt$nCI[io]
        for(j in 1:nCI)  
          arrows(CI.from$x[io, j], CI.from$y[io, j], 
                 CI.to$x[io, j], CI.to$y[io, j], 
                 angle = 90, length = CIpt$CIpar$alen[io], code = 3, 
                 col = CIpt$CIpar$col[io], 
                 lty = CIpt$CIpar$lty[io], 
                 lwd = CIpt$CIpar$lwd[io])

        if (CIpt$dopoints) 
          points(CI.mid$x[io], CI.mid$y[io], 
                 col = CIpt$col[io], pch = CIpt$pch[io],
                 cex = CIpt$cex[io], bg = CIpt$bg[io])
                                                                
    } else if (type[ii] == 2) {
       io <- i - LCI 
       polygon(pol$x[, io], pol$y[, io], 
               lwd = poly$lwd[io],
               lty = poly$lty[io], 
               border = poly$border[io], 
               col = poly$col[io])

    } else if (type[ii] == 3) {
       io <- i - LCP
       segments(segm.from$x[io], segm.from$y[io],  
                segm.to$x[io], segm.to$y[io], 
                col = segm$col[io], 
                lwd = segm$lwd[io], 
                lty = segm$lty[io])

    } else if (type[ii] == 4) {
       io <- i - LCPS
       text(x = lab$x[io], y = lab$y[io], 
            labels = labels$labels[io], 
            col = labels$col[io],
            adj = labels$adj[io], 
            cex = labels$cex[io],
            font = labels$font[io],
            srt  = labels$srt[io[1]])

    } else if (type[ii] == 5) {
       io <- i - LCPSL
         ArrType (arr.from$x[io], arr.from$y[io], 
            arr.to$x[io], arr.to$y[io], arr$length[io], 
            arr$angle[io], arr$code[io], 
            col = arr$col[io], type = arr$type[io], 
            lwd = arr$lwd[io], lty = arr$lty[io])

   } else if (type[ii] == 6) {
      io <- i - LCPSLA
      dp <- extractdots(other$dot, i)
      stop("structure 'other' not defined")  
      }
    
    }
  
  sapply(FUN = plotit, 1:length(sortlist))
#    for (i in 1:length(sortlist)) plotit(i)

}

selectplist <- function (plist, SS) {
 #
 
  if (! is.null(plist$pt)) {
    pt  <- plist$pt                  
    ipt <- with (pt, SS(x.mid, y.mid, z.mid))
    if (sum(ipt) > 0)
      plist$pt <- with(pt, list(x.mid = x.mid[ipt], y.mid = y.mid[ipt],
        z.mid = z.mid[ipt], col = col[ipt], pch = pch[ipt],
        cex = cex[ipt], bg = bg[ipt], alpha = alpha[ipt], proj = proj[ipt]))
    else
      plist$pt <- NULL
  }
 #
  if (!is.null(plist$CIpt)) {
    CIpt <- plist$CIpt
    ipt <- with (CIpt, SS(x.mid, y.mid, z.mid))
    if (sum(ipt) > 0)
      plist$CIpt <- with (CIpt, list(x.to = x.to[ipt,], y.to = y.to[ipt,],
       z.to = z.to[ipt,], x.from = x.from [ipt,], y.from = y.from[ipt,],
       z.from = z.from[ipt,], nCI = nCI[ipt], x.mid = x.mid[ipt],
       y.mid = y.mid[ipt], z.mid = z.mid[ipt], length = CIpt$length[ipt],
       col= col[ipt], pch = pch[ipt],
       bg = bg[ipt], cex = cex[ipt], alpha = alpha[ipt], proj = proj[ipt],
       CIpar = list (col = CIpar$col[ipt], lwd = CIpar$lwd[ipt],
        lty = CIpar$lty[ipt], alen = CIpar$alen[ipt])))
    else
      plist$CIpt <- NULL
  }

 # polygons
  if (length(plist$poly$x) > 0) {
    xm <- colMeans(plist$poly$x, na.rm = TRUE)
    ym <- colMeans(plist$poly$y, na.rm = TRUE)
    zm <- colMeans(plist$poly$z, na.rm = TRUE)
    remove <- NULL
    if (any(is.nan(xm)))
      remove <- c(remove, which(is.nan(xm)))
    if (any(is.nan(ym)))
      remove <- c(remove, which(is.nan(ym)))
    if (any(is.nan(zm)))
      remove <- c(remove, which(is.nan(zm)))
    if (! is.null(remove)) {
      remove <- unique(remove)
      xm[is.nan(xm)] <- 0
      ym[is.nan(ym)] <- 0
      zm[is.nan(zm)] <- 0
    }      
    ip <- SS(xm, ym, zm)
    if (! is.null(remove)) 
      ip[remove] <- FALSE

    if (sum(ip) > 0)
      plist$poly <- with (plist$poly, list(x = as.matrix(x[,ip]), 
         y = as.matrix(y[,ip]), z = as.matrix(z[,ip]),
        col = col[ip], border = border[ip], lwd = lwd[ip], lty = lty[ip], 
        alpha = alpha[ip], isimg = isimg[ip], proj = proj[ip]))
    else
      plist$poly <- NULL
  } else
      plist$poly <- NULL

 # images
  imgxrange <- NULL
  imgyrange <- NULL
  imgzrange <- NULL
  
  if (length(plist$img)> 0 ) {
    for (i in plist$imgnr:1) {
      img <- plist$img[[i]]
      
     # because col has one row and column less than x, y, z 
      Col <- img$col
      if (nrow(Col) != nrow(img$z))
        Col <- rbind(Col, Col[nrow(Col),])
      if (ncol(Col) != ncol(img$z))
        Col <- cbind(Col, Col[,ncol(Col)])
      
      # expand all values of x and y
      if (is.vector(img$x)) {
        XY <- mesh(img$x, img$y)
        img$x <- XY$x; img$y <- XY$y
      }
      Nx <- nrow(img$z)
      Ny <- ncol(img$z)

      remove <- NULL
      if (any(is.na(img$x)))
        remove <- c(remove, which(is.na(img$x)))
      if (any(is.na(img$y)))
        remove <- c(remove, which(is.na(img$y)))
      if (any(is.na(img$z)))
        remove <- c(remove, which(is.na(img$z)))
      if (! is.null(remove)) {
        remove <- unique(remove)
        img$x[is.na(img$x)] <- 0
        img$y[is.na(img$y)] <- 0
        img$z[is.na(img$z)] <- 0
      }      

      ipt <- SS(img$x, img$y, img$z)
      if (! is.null(remove))  {
        ipt[remove] <- FALSE
        img$x[remove] <- NA
        img$y[remove] <- NA
        img$z[remove] <- NA
      }

      if (sum(ipt) > 0) {
        zset <- min(img$z[ipt], na.rm = TRUE)
        Select <- img$z; Select[] <- ipt
        Noselect <- which(!ipt)
        Col[Noselect] <- "transparent"
        img$z[Noselect] <- NA #zset
        arrsel <- which (Select == 1, arr.ind = TRUE)
        xr <- range(arrsel[,1])
        yr <- range(arrsel[,2])
  
        xsel <- xr[1] : xr[2]
        ysel <- yr[1] : yr[2]
        
        if (is.vector(plist$img[[i]]$x)) 
          plist$img[[i]]$x <- plist$img[[i]]$x[xsel]
        else
          plist$img[[i]]$x <- img$x[xsel, ysel]

        if (is.vector(plist$img[[i]]$x)) 
          plist$img[[i]]$y <- plist$img[[i]]$y[ysel]
        else
          plist$img[[i]]$y <- img$y[xsel, ysel]
        
        plist$img[[i]]$z <- img$z[xsel, ysel]
        plist$img[[i]]$col <- Col[xsel, ysel]
        imgxrange <- range( c(imgxrange, plist$img[[i]]$x), na.rm = TRUE) 
        imgyrange <- range( c(imgyrange, plist$img[[i]]$y), na.rm = TRUE) 
        imgzrange <- range( c(imgzrange, plist$img[[i]]$z), na.rm = TRUE) 


      } else {
        plist$img[[i]] <- NULL    
        plist$imgnr <- plist$imgnr - 1
      }  
    }
  }
  
  if (!is.null(plist$labels)) {
    labels <- plist$labels
    il <- with (labels, SS(x, y, z))
    if (sum(il) > 0)
      plist$labels <- with (labels, list(x = x[il], y = y[il], z = z[il],
      labels = labels[il], adj = adj[il], cex = cex[il],
      col = col[il], font = font[il], alpha = alpha[il], proj = proj[il]))
    else
      plist$labels <- NULL
  }

  if (!is.null(plist$segm)) {
    segm <- plist$segm
    is <- with (segm, SS(colMeans(rbind(as.vector(x.from), as.vector(x.to))),
      colMeans(rbind(as.vector(y.from), as.vector(y.to))), colMeans(rbind(as.vector(z.from), as.vector(z.to)))))
    if (sum(is) > 0)
      plist$segm <- with (segm, list(x.from = as.vector(x.from)[is], 
        y.from = as.vector(y.from)[is], z.from = as.vector(z.from)[is], 
        x.to = as.vector(x.to)[is], y.to = as.vector(y.to)[is], 
        z.to = as.vector(z.to)[is],
        proj = proj[is], lwd  = lwd[is], lty  = lty[is], alpha = alpha[is], col  = col[is]))
    else
      plist$segm <- NULL
  }

  if (!is.null(plist$arr)) {
    arr <- plist$arr
    is <- with (arr, SS(colMeans(rbind(as.vector(x.from), as.vector(x.to))),
      colMeans(rbind(as.vector(y.from), as.vector(y.to))), colMeans(rbind(as.vector(z.from), as.vector(z.to)))))
    if (sum(is) > 0)
      plist$arr <- with (arr, list(x.from = as.vector(x.from)[is], 
        y.from = as.vector(y.from)[is], z.from = as.vector(z.from)[is], 
        x.to = as.vector(x.to)[is], y.to = as.vector(y.to)[is], 
        z.to = as.vector(z.to)[is],
        proj = proj[is], lwd  = lwd[is], lty  = lty[is], col  = col[is],
        code = code[is], angle = angle[is], length = length[is], type = type[is], alpha = alpha[is]))
    else
      plist$arr <- NULL
  }

  # new ranges 
   xs <- c(plist$pt$x.mid, plist$CIpt$x.to, plist$CIpt$x.from,
      plist$poly$x, plist$labels$x, plist$segm$x.from, plist$segm$x.to,
      plist$arr$x.from, plist$arr$x.to, imgxrange)
   if (length(xs) > 0)
     plist$xlim <- newlim(xs)

   ys <- c(plist$pt$y.mid, plist$CIpt$y.to, plist$CIpt$y.from,
      plist$poly$y, plist$labels$y, plist$segm$y.from, plist$segm$y.to,
      plist$arr$y.from, plist$arr$y.to, imgyrange)   
   if (length(ys) > 0)
     plist$ylim <- newlim(ys)

   zs <- c(plist$pt$z.mid, plist$CIpt$z.to, plist$CIpt$z.from,
      plist$poly$z, plist$labels$z, plist$segm$z.from, plist$segm$z.to,
      plist$arr$z.from, plist$arr$z.to, imgzrange)
   if (length(zs) > 0)
     plist$zlim <- newlim(zs)

  return (plist)
}

newlim <- function(xx) {

  lim <- range(xx, na.rm = TRUE)

  if (any(is.infinite(lim))) return(c(-0.1,0.1))
  if (diff(lim) == 0)
    lim <- lim * c(0.8, 1.2)
  if (diff(lim) == 0)
    lim <- lim + c(-0.1, 0.1)
  return(lim)
}
