lines3D <- function(x, y, z, ...) {
  dot <- list(...)
  if (is.null(dot$type)) 
    dot$type <- "l"
  plist <- do.call("scatter3D", c(alist(x, y, z), dot))
  invisible(plist) 
}

points3D <- function(x, y, z,  ...) {
  dot <- list(...)
  if (is.null(dot$type)) 
    dot$type <- "p"
  plist <- do.call("scatter3D", c(alist(x, y, z), dot))
  invisible(plist) 
}

## =============================================================================
## scatterplot in 3D
## =============================================================================
    
scatter3D <- function(x, y, z, ..., colvar = z, 
                      phi = 40, theta = 40,
                      col = NULL, NAcol = "white", breaks = NULL,
                      colkey = NULL, 
                      panel.first = NULL, clim = NULL, clab = NULL, 
                      bty = "b", CI = NULL, surf = NULL, 
                      add = FALSE, plot = TRUE) {

  plist <- initplist(add)

  dot <- splitdotpersp(list(...), bty, NULL, x, y, z, plist = plist, breaks = breaks)

  x <- as.vector(x)
  y <- as.vector(y)
  z <- as.vector(z)
  len <- length(x)

  if (length(y) != len)
    stop("'y' should have same length as 'x'")
  if (length(z) != len)
    stop("'z' should have same length as 'x'")

  if (len > 1 & ispresent(colvar)) {
  
    if (length(colvar) != len)
      stop("'colvar' should have same length as 'x', 'y' and 'z'")

    colvar <- as.vector(colvar)
    
    if (is.null(clim)) 
      clim <- range(colvar, na.rm = TRUE)
    
    if (dot$clog) {
      colvar <- log(colvar)
      clim <- log(clim) 
    }

    if (is.null(col))
      if (is.null(breaks))
        col <- jet.col(100)
      else
        col <- jet.col(length(breaks)-1)

    iscolkey <- is.colkey(colkey, col)
    if (iscolkey) 
      colkey <- check.colkey(colkey)

    if (! is.null(dot$alpha)) 
      col <- setalpha(col, dot$alpha)
    
    Col <- variablecol(colvar, col, NAcol, clim, breaks)
    if (length(Col) == 1)
      Col <- rep(Col, length.out = len)

  } else {
    if (is.null(col))
      col <- "black"
    if (! is.null(dot$alpha)) 
      col <- setalpha(col, dot$alpha)
    Col <- rep(col, length.out = len)  
    iscolkey <- FALSE
  }
  
  if (is.null(plist)) {
    do.call("perspbox", c(alist(x = range(x), y = range(y), 
             z =  range(z, na.rm = TRUE),
             phi = phi, theta = theta, plot = plot, colkey = colkey, col = col), 
             dot$persp))
    plist <- getplist()
  } 
  breaks <- check.breaks(breaks, col)

 # droplines with a fitted surface  
  fit <- NULL
  if (! is.null(surf)) {
    if (! is.list(surf))
      stop("'surf' should be a 'list' or 'NULL'")
    fit <- surf$fit
    surf$fit <- NULL    
  }  
  
  if (! is.null(fit)){
    if (! is.null(CI)) {
      if (! is.null(CI$z))
        stop("cannot combine a confidence interval (CI) on 'z' with 'fit' in 'surf'")
    } else 
      CI <- list()
    
    if (length(fit) != length(z))
      stop("'fit', argument of 'surf'  should be of equal size of 'z'")  
    disttoz <- fit - z
    CIz <- matrix(ncol = 2, data = c(-disttoz, disttoz))
    CIz[CIz > 0] <- 0 
    CI$z = CIz
    CI$alen = 0
  } 

 # confidence intervals
  isCI <- is.list(CI)
  if (isCI) 
    CI <- check.CI(CI, len, 3)
   
  if (is.null(CI) & any (Col == "transparent")){
    ii <- Col != "transparent"
    x <- x[ii]
    y <- y[ii]
    z <- z[ii]
    Col <- Col[ii]
    len <- length(x)
  }
  
  if (is.function(panel.first)) 
    panel.first(plist$mat)  

  if (! is.null(surf)) {
    
    if (is.null(surf$colvar))
      surf$colvar <- surf$z

    if (is.null(surf$breaks))
      surf$breaks <- breaks

    if (is.null(surf[["col"]])) {
      surf$col <- col
      if (is.null(surf$clim))  
        surf$clim <- clim
    }  
    if (is.null(surf$clim))  
      surf$clim <- range(surf$colvar)

    surf$colvar[surf$colvar < min(surf$clim)]  <- NA
    surf$colvar[surf$colvar > max(surf$clim)]  <- NA
 
    surf$z[surf$z < dot$persp$zlim[1]]  <- NA
    surf$z[surf$z > dot$persp$zlim[2]]  <- NA

    spoly <- do.call("addimg", c(alist(poly = NULL, plist = plist), surf))
  
  } else
    spoly <- NULL

  sseg <- NULL  # segments
          
  dtype <- dot$points$type
  dot$points$type <- NULL

  lwd <- dot$points$lwd ; if (is.null(lwd)) lwd <- 1
  lty <- dot$points$lty ; if (is.null(lty)) lty <- 1
  alpha <- dot$alpha; if (is.null(alpha)) alpha <- NA
  alpha <- rep(alpha, length.out = len)

  if (is.null(dtype))
    dtype <- "p"

 # droplines
  if (dtype == "h") {   
    zmin <- dot$persp$zlim[1]
    Proj   <- project(x, y, 0.5 *(z + zmin), plist, FALSE)

    sseg <- list(x.from = x, 
                 x.to   = x,
                 y.from = y, 
                 y.to   = y,      
                 z.from = rep(zmin, len), 
                 z.to   = z,                            
                 col    = Col,
                 lwd    = rep(lwd , length.out = len),
                 lty    = rep(lty , length.out = len),
                 alpha  = alpha,
                 proj   = Proj)
  
    class(sseg) <- "segments"

 # segments between points 
  } else if (dtype %in% c("b", "l", "o")) {  
   # segment color = mean of point colors
    LCol <- Col
    if (length(LCol) > 1) {
      LCol <- cbind(Col[-1], Col[-len])
      LCol <- apply(LCol, MARGIN = 1, FUN = MeanColors)
      if (! is.null(dot$alpha)) 
        LCol <- setalpha(LCol, dot$alpha)
      
    }

    Proj   <- project(0.5*(x[-len]+x[-1]), 0.5*(y[-len]+y[-1]), 
                     0.5*(z[-len]+z[-1]), plist)
    sseg <- list(x.from = c(sseg$x.from, x[-len]), 
                 x.to   = c(sseg$x.to,   x[-1]),
                 y.from = c(sseg$y.from, y[-len]), 
                 y.to   = c(sseg$y.to,   y[-1]),                                  
                 z.from = c(sseg$z.from, z[-len]), 
                 z.to   = c(sseg$z.to,   z[-1]),                                  
                 col    = c(sseg$col, LCol),
                 lwd    = c(sseg$lwd, rep(lwd , length.out = len-1)),
                 lty    = c(sseg$lty, rep(lty , length.out = len-1)),
                 alpha  = c(sseg$alpha, alpha),
                 proj   = c(sseg$proj, Proj))
    class(sseg) <- "segments"
  }
  
  pch <- dot$points$pch 
  if (is.null(pch)) 
    pch <- 1
  bg  <- dot$points$bg  
  if (is.null(bg)) 
    bg <- 1
  cex <- dot$points$cex 
  if (is.null(cex)) 
    cex <- 1

  if (dtype == "l") 
    dopoints <- FALSE
  else 
    dopoints <- TRUE  

  CIpt <- NULL
  pt <- NULL

  if (! is.null(CI)) {
 # points and confidence intervals
    CIpt <- list(x.from = NULL, 
                 y.from = NULL, 
                 z.from = NULL,
                 x.to = NULL, 
                 y.to = NULL, 
                 z.to = NULL,
                 x.mid = x, 
                 y.mid = y, 
                 z.mid = z,
                 col = Col,
                 pch = rep(pch, length.out = len),
                 cex = rep(cex, length.out = len),
                 bg  = rep(bg, length.out = len),
                 alpha  = alpha
                 )
    class(CIpt) <- "CIpt"

    CIpt$CIpar <- CI   #[c("lty", "lwd", "col")]

   # length of arrow head
    CIpt$CIpar$alen <- rep(mean(par("fin")) * CI$alen, len)
    
    if (is.null(CIpt$CIpar$col))
      CIpt$CIpar$col <- Col

    if (is.null(CIpt$CIpar$lwd))
      CIpt$CIpar$lwd <- 1

    if (is.null(CIpt$CIpar$lty))
      CIpt$CIpar$lty <- 1

   # define par settings for all points 
    CIpt$CIpar <- setdots(CIpt$CIpar, len)

    addCI <- function(x.from, y.from, z.from, 
                      x.to, y.to, z.to, tr)  {
      tr$x.from <- cbind(tr$x.from, x.from) 
      tr$y.from <- cbind(tr$y.from, y.from)    
      tr$z.from <- cbind(tr$z.from, z.from)    
      tr$x.to   <- cbind(tr$x.to, x.to)
      tr$y.to   <- cbind(tr$y.to, y.to)
      tr$z.to   <- cbind(tr$z.to, z.to)
      tr
    }

    if (! is.null(CI$x))  # CI in x direction
      CIpt <- addCI(x - CI$x[ ,1], y, z, 
                    x + CI$x[ ,2], y, z, CIpt)
  
    if (! is.null(CI$y)) 
      CIpt <- addCI(x, y - CI$y[ ,1], z, 
                    x, y + CI$y[ ,2], z, CIpt)

    if (! is.null(CI$z))     
      CIpt <- addCI(x, y, z - CI$z[ ,1], 
                    x, y, z + CI$z[ ,2], CIpt)

   # upscale to 3 columns each
    nCI <- ncol(CIpt$x.from)
    CIpt$nCI <- rep(nCI, len)
    if (nCI < 3) {
      n0 <- matrix(nrow = len, ncol = 3 - nCI, data = NA)
      CIpt$x.from <- cbind(CIpt$x.from, n0) 
      CIpt$y.from <- cbind(CIpt$y.from, n0)    
      CIpt$z.from <- cbind(CIpt$z.from, n0)    
      CIpt$x.to   <- cbind(CIpt$x.to, n0)
      CIpt$y.to   <- cbind(CIpt$y.to, n0)
      CIpt$z.to   <- cbind(CIpt$z.to, n0)
    }

    CIpt$dopoints <- TRUE
 
  } else if (dopoints) {
   # points - remove NAs
    ii <- which(!is.na(x) & !is.na(y) & !is.na(z))
    pt <- list(x.mid = x[ii], 
               y.mid = y[ii], 
               z.mid = z[ii],
               col = Col[ii],
               pch = rep(pch, length.out = length(ii)),
               lwd = rep(lwd, length.out = length(ii)),
               cex = rep(cex, length.out = length(ii)),
               bg = rep(bg, length.out = length(ii)),
               alpha = alpha
               )
    
    class(pt) <- "pt"
  }

  Proj   <- project(x, y, z, plist)

  if (! is.null(CI)) 
    CIpt$proj <- Proj
  else if (! is.null(pt))
    pt$proj <- Proj
    
  if (iscolkey) 
    plist <- plistcolkey(plist, colkey, col, clim, clab, dot$clog,
      type = "scatter3D", breaks = breaks)

 # plot it
  plist <- plot.struct.3D(plist, pt = pt, CIpt = CIpt, 
               poly = spoly, segm = sseg, plot = plot)  
  
  setplist(plist)
  invisible(plist$mat)
}

# Note: check.CI is in file scatter.R

