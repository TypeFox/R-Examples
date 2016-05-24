## =============================================================================
## 3-D arrows function
## =============================================================================

arrows3D  <- function(x0, y0, z0, x1 = x0, y1 = y0, z1 = z0,
                    ..., colvar = NULL, phi = 40, theta = 40,
                    col = NULL, NAcol = "white", breaks = NULL,
                    colkey = NULL, panel.first = NULL,
                    clim = NULL, clab = NULL, bty = "b", type = "triangle",
                    add = FALSE, plot = TRUE)  {
  plist <- initplist(add)

  dot  <- splitdotpersp(list(...), bty, NULL, 
    c(x0, x1), c(y0, y1), c(z0, z1), plist = plist, breaks = breaks)

  len <- length(x0)
  if (length(y0) != len)
    stop("'y0' should have same length as 'x0'")
  if (length(z0) != len)
    stop("'z0' should have same length as 'x0'")
  if (length(x1) != len)
    stop("'x1' should have same length as 'x0'")
  if (length(y1) != len)
    stop("'y1' should have same length as 'x0'")
  if (length(z1) != len)
    stop("'z1' should have same length as 'x0'")

  if (is.null(col) & is.null(breaks))
    col <- jet.col(100)
  else if (is.null(col))
    col <- jet.col(length(breaks)-1)

  breaks <- check.breaks(breaks, col)
  if (ispresent(colvar)) {
    if (length(colvar) != len)
      stop("'colvar' should have same length as 'x0', 'y0' and 'z0'")

    if (length(col) == 1)
      col <- c(col, col)

    if (is.null(clim))
      clim <- range(colvar, na.rm = TRUE)

    if (dot$clog) {                    
      colvar <- log(colvar)
      clim <- log(clim)
    }

    iscolkey <- is.colkey(colkey, col) 
    if (iscolkey) 
      colkey <- check.colkey(colkey)
    if (! is.null(dot$alpha)) 
      col <- setalpha(col, dot$alpha)
    Col <- variablecol(colvar, col, NAcol, clim, breaks)

  } else {
    if (is.null(col))
      col <- "black"
    if (! is.null(dot$alpha)) 
      col <- setalpha(col, dot$alpha)
      
    Col <- rep(col, length.out = len)  
    iscolkey <- FALSE
  }   

  if (is.null(plist)) {
    do.call("perspbox", 
       c(alist(x = range(c(x0, x1)), y = range(c(y0, y1)), 
               z = range(c(z0, z1)), 
               phi = phi, theta = theta, plot = plot, col = col), dot$persp))
    plist <-  getplist()
  }  
  
  if (is.function(panel.first)) 
    panel.first(plist$mat)
  
  length <- dot$points$length
  if (is.null(length)) 
    length <- 0.2
  
  angle <- dot$points$angle
  if (is.null(angle)) 
    angle <- 30
  
  code <- dot$points$code
  if (is.null(code)) 
    code <- 2
  
  lwd <- dot$points$lwd
  if (is.null(lwd)) 
    lwd <- 1
  
  lty <- dot$points$lty
  if (is.null(lty)) 
    lty <- 1

  Proj <- project (0.5*(x0 + x1), 0.5*(y0 + y1), 0.5*(z0 + z1), plist)
  alpha <- dot$alpha; if (is.null(alpha)) alpha <- NA
  alpha <- rep(alpha, length.out = len)

  arr  <- list(x.from = x0, 
               x.to   = x1,
               y.from = y0, 
               y.to   = y1,                                  
               z.from = z0, 
               z.to   = z1,                                  
               col    = Col,
               length = rep(length, length.out = len),
               code   = rep(code  , length.out = len),
               angle  = rep(angle , length.out = len),
               lwd    = rep(lwd   , length.out = len),
               lty    = rep(lty   , length.out = len),
               type   = rep(type  , length.out = len),
               alpha  = alpha,
               proj   = Proj)
  class(arr) <- "arr"

  if (iscolkey) 
    plist <- plistcolkey(plist, colkey, col, clim, clab, 
      dot$clog, type = "arrows3D", breaks = breaks)

  plist <- plot.struct.3D(plist, arr = arr, plot = plot)  

  setplist(plist)   
  invisible(plist$mat)
}


## This part is adapted from (my) package "shape"...

ArrType <- function (x0, y0, x1, y1, length = 0.4, angle = 30, 
    code = 2, adj = 1, type = "curved", col = "black", 
    lty = 1, lwd = 1, ...) {
  if (length(unique(type)) > 1) {
    nr <- length(x0)
    length <- rep(length, length.out = nr)
    angle <- rep(angle, length.out = nr)
    code <- rep(code, length.out = nr)
    adj <- rep(adj, length.out = nr)
    col <- rep(col, length.out = nr)
    lty <- rep(lty, length.out = nr)
    lwd <- rep(lwd, length.out = nr)
    for (t in unique(type)) {
      ii <- which(type == t)
       Arrow (x0[ii], y0[ii], x1[ii], y1[ii], length[ii], angle[ii], 
         code[ii], adj[ii], type = t, col[ii], lty[ii], lwd[ii],  ...)     
    }
  } else
     Arrow(x0, y0, x1, y1, length, angle, 
       code, adj, type = type[1], col, lty, lwd, ...)   
}
    
Arrow <- function (x0, y0, x1, y1, length = 0.4, angle = 30, 
    code = 2, adj = 1, type = "curved", col = "black", 
    lty = 1, lwd = 1, ...)  {
  if (all (is.na(x0)) | all(length == 0))
    return()
  
  sel <- function (val, ii) {
    if (length(val) == 1) 
      return(val)
    else
      return(val[ii])
  }    
  ii <- which(length <= 0)
  if (length(ii) > 0) {
    col <- rep(col, length.out = length(x0))
    col[ii] <- "transparent"  # to avoid from seeing a "." (segments, arrows)
  }  
  if (type == "simple") {
    arrows(x0, y0, x1, y1, code = code, length = length/2.54, 
          angle = angle, lty = lty, col = col, lwd = lwd, ...)
    return()
  }
  width <- 2* tan(angle/180) * length
  segments(x0, y0, x1, y1, col = col, lty = lty, lwd = lwd, ...)
  user <- par("usr")
  pin <- par("pin")
  pin <- pin/max(pin)
  sy <- (user[4] - user[3])/pin[2]
  sx <- (user[2] - user[1])/pin[1]
  angle <- atan((y1 - y0)/(x1 - x0) * sx/sy)/pi * 180
  angle[is.nan(angle)] <- 0
  angle[x1 < x0] <- 180 + angle[x1 < x0]
  xx <- x1
  yy <- y1
  if (length(code) > 1) {
    code <- rep(code, length.out = length(xx)) 
    ii <- which (code %in% c(2, 3)) 
    
    if (length(ii) > 0) 
        Arrow.head(x0 = xx[ii], y0 = yy[ii], angle = angle[ii], 
        col = sel(col, ii), adj = sel(adj, ii), lty = sel(lty, ii), 
        length = sel(length, ii), width = sel(width, ii), 
        type = sel(type, ii), lwd = sel(lwd, ii))
    ii <- which (code %in% c(1, 3)) 
    if (length(ii) > 0) {
        angle[ii] <- 180 + angle[ii]
        xx[ii] <- x0[ii]
        yy[ii] <- y0[ii]
      Arrow.head(x0 = xx, y0 = yy, angle = angle, col = col,  
        adj = adj, lty = lty, length = length, 
        width = width, type = type, lwd = lwd)
    }
  } else {
    if (code %in% c(2, 3)) 
        Arrow.head(x0 = xx, y0 = yy, angle = angle, col = col, 
          adj = adj, lty = lty, length = length, 
          width = width, type = type, lwd = lwd)
    if (code %in% c(1, 3)) {
        angle <- 180 + angle
        xx <- x0
        yy <- y0
        Arrow.head(x0 = xx, y0 = yy, angle = angle, col = col,  
          adj = adj, lty = lty, length = length, 
          width = width, type = type, lwd = lwd)
    }
  }
}

Arrow.head <- function (x0, y0, angle = 0, 
    length = 0.4, width = length/2, 
    adj = 0.5, type = "triangle", col = "black", lty = 1, 
    lwd = 2, npoint = 5) {
  
  adj <- adj[1]
    
  if (type == "curved") {
    rad <- 0.7
    len <- 0.25 * pi
    mid <- c(0, rad)
    x <- seq(1.5 * pi + len, 1.5 * pi, length.out = npoint)
    rr <- cbind(mid[1] - rad * cos(x), mid[2] + rad * sin(x))
    mid <- c(0, -rad)
    x <- rev(x)
    rr <- rbind(rr, cbind(mid[1] - rad * cos(x), mid[2] - 
        rad * sin(x)))
    mid <- c(rr[nrow(rr), 1], 0)
    rd <- rr[1, 2]
    x <- seq(pi/2, 3 * pi/2, length.out = 3 * npoint)
    rr <- rbind(rr, cbind(mid[1] - rd * 0.5 * cos(x), mid[2] - 
            0.5 * rd * sin(x)))
    rr[, 1] <- rr[, 1] * 2.6
    rr[, 2] <- rr[, 2] * 3.45
  }
  else if (type %in% c("cone", "triangle")) {
    x <- c(-0.2, 0, -0.2)
    y <- c(-0.2, 0, 0.2)
    rr <- 6.22 * cbind(x, y)
  }
  else if (type %in% c("circle", "ellipse")) {
    if (type == "circle") 
      width <- length
    rad <- 0.1
    mid <- c(-rad, 0)
    x <- seq(0, 2 * pi, length.out = 15 * npoint)
    rr <- 6.22 * cbind(mid[1] + rad * sin(x), mid[2] + rad * 
            cos(x))
  }
  if (adj == 0.5) 
    rr[, 1] <- rr[, 1] - min(rr[, 1])/2
  if (adj == 0) 
    rr[, 1] <- rr[, 1] - min(rr[, 1])
  
  user <- par("usr")
  pcm <- par("pin") * 2.54
  sy <- (user[4] - user[3])/pcm[2]
  sx <- (user[2] - user[1])/pcm[1]
  nr <- max(length(x0), length(y0), length(angle), length(length), 
        length(width), length(lty), length(col))
  if (nr > 1) {
    x0 <- rep(x0, length.out = nr)
    y0 <- rep(y0, length.out = nr)
    angle <- rep(angle, length.out = nr)
    length <- rep(length, length.out = nr)
    width <- rep(width, length.out = nr)
    col <- rep(col, length.out = nr)
    lty <- rep(lty, length.out = nr)
  }
  RR <- rr
  for (i in 1:nr) {
    dx <- rr[, 1] * length[i]
    dy <- rr[, 2] * width[i]
    angpi <- angle[i]/180 * pi
    cosa <- cos(angpi)
    sina <- sin(angpi)
    RR[, 1] <- cosa * dx - sina * dy
    RR[, 2] <- sina * dx + cosa * dy
    RR[, 1] <- x0[i] + RR[, 1] * sx
    RR[, 2] <- y0[i] + RR[, 2] * sy
    polygon(RR, col = col[i], border = col[i], lty = lty[i], 
            lwd = lwd)
  }
}

