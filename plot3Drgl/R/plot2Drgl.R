namespersp <- c("xlim", "ylim", "zlim", "xlab", "ylab", "zlab",
        "main", "sub", "r", "d", "scale", "expand", "box", "axes", "type",
        "nticks", "ticktype", "col.ticks", "lwd.ticks", "log", "shade",
        "cex.axis", "col.axis", "font.axis", "col.panel",
        "bty", "lwd.panel", "col.grid", "lwd.grid")

## =============================================================================

plot2Drgl <- function(func.name, x, y, colvar, col, NAcol = "white", breaks, clim, add,
  clab = NULL, colkey = FALSE, namesextra = NULL, z = NULL, ...) {
  
# ------------------------------------------------------------------------------

  dots <- list(...)
# image has as argument 'zlim' -> is clim here
  if (func.name == "image3D" & ! is.null(dots$zlim)) {
    if (!is.null(clim))
      stop ("'zlim' and 'clim' cannot both be specified in image2Drgl")
    clim <- dots$zlim
    dots$zlim <- NULL
  }
  dots$expand <- dots$ticktype <- dots$zlab <- dots$bty <- NULL
  dots$box <- TRUE
  namesscat <- c(namespersp, namesextra)
  dotpersp <- dots[names(dots) %in% namesscat]
  dots$add <- add
  dots <- testdots(dots, func.name)
  if (is.null(z)) 
    z <- rep(1, length.out = length(x))
  plot3D:::refresh(FALSE)
  do.call(func.name, c(alist(x = x, y = y, z = z, colkey = colkey,
    colvar = colvar, col = col, NAcol = NAcol, breaks = breaks,
    clim = clim, clab = clab, bty = "b",
    plot = FALSE, add = add, zlab = "", ticktype = "simple"), dotpersp))
##  plot3D:::refresh(TRUE)

  dots
}

testdots <- function(dots, func.name) {
  if (! is.null(dots$rgltype)) {
    if (dots$rgltype == "new") {
      dots$new <- TRUE
      dots$add <- FALSE
    } else if (dots$rgltype == "add") {
      dots$add <- TRUE
      dots$new <- FALSE
    } else if (dots$rgltype == "rep") {
      dots$new <- FALSE
      dots$add <- TRUE
      ids <- rgl.ids()
      toreplace <- NULL
      if (func.name%in% c("persp3D", "slice3D","slicecont3D","isosurf3D",
        "surf3D","spheresurf3D","image3D")) toreplace <- "surface"
      else if (func.name%in% c("ribbon3D", "hist3D","box3D","rect3D"))
        toreplace <- "quads"
      else if (func.name%in% c("scatter3D", "points3D", "voxel3D"))
        toreplace <- "points"
      else if (func.name%in% c("lines3D", "segments3D", "border3D", "contour3D"))
        toreplace <- "lines"
      else if (func.name == "text3D")
        toreplace <- "text"
      rgl.pop(type = "shapes", id = ids[ids$type == toreplace, 1])
      plist <- getplist()
      plist$imgnr<-plist$img<-plist$segm<-plist$pt<-plist$CIpt<-plist$labels <- NULL
      setplist(plist)
    }
    dots$rgltype <- NULL
  }
  if (is.null(dots$add))
    dots$add <- FALSE
  if(dots$add)
    dots$new <- FALSE
  dots
}
# ------------------------------------------------------------------------------

plot2Drglbis <- function(func.name, x0, y0, x1, y1, dz, colvar, col,
  NAcol = "white", breaks = NULL, clim, add, clab = NULL,
  colkey = FALSE, namesextra = NULL, z0 = NULL, z1 = NULL, ...) {

# ------------------------------------------------------------------------------
# same but with arguments x0,... z1

  dots <- list(...)
  z0 <- z1 <- rep(1 + dz, length.out = length(x0))

  dots$expand <- dots$ticktype <- dots$zlab <- dots$bty <- NULL
  dots$box <- TRUE
  namesscat <- c(namespersp, namesextra)
  dotpersp <- dots[names(dots) %in% namesscat]
  dots$add <- add
  dots <- testdots(dots, func.name)

  plot3D:::refresh(FALSE)
  do.call(func.name, c(alist(x0 = x0, y0 = y0, z0 = z0,
    x1 = x1, y1 = y1, z1 = z1, colkey = colkey,
    colvar = colvar, col = col, NAcol = NAcol, breaks = breaks,
    clim = clim, bty = "b", clab = clab,
    plot = FALSE, add = add, zlab = "", ticktype = "simple"), dotpersp))
##  plot3D:::refresh(TRUE)

  dots
}

## =============================================================================

finishplotrgl <- function(dots, namesextra = NULL) {
#  new <- dots $new
  plist <- getplist()
  plist$type <- "3D" 
  setplist(plist)
  
  do.call("plotrgl", c(dots[!names(dots) %in% c(namespersp, namesextra)]))

  mouseNULL <- is.null(dots$mouseMode)
  if (mouseNULL) 
    dots$mouseMode <- c("zoom", "zoom", "zoom")

  par3d(mouseMode = dots$mouseMode)
  if (mouseNULL) 
    pan3d(3)           # from the help of rgl.setMouseCallbacks

  view3d(phi = 0, fov = 0)
#  decorate3d(zlab = "")
  
  axis3d("x")
  axis3d("y") 
  pp <- getplist()
#  pp$type <- "2D"
  pp$dot$ticktype <- "simple"
  pp$rgl$userMatrix <- par3d("userMatrix")
  setplist(pp) 
}
                                       
