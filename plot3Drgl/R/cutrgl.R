cutrgl <- function(...) {
  SS <- select3d()
  plist <- selectplist(getplist(), SS)
  rgl.clear()
  dots <- list(...)
  light <- dots$lighting
  if (is.null(light))
    light <- getplist()$rgl$lighting
  smooth <- dots$smooth
  if (is.null(smooth))
    smooth <- getplist()$rgl$smooth
  dots$lighting <- dots$new <- dots$smooth <- NULL
  do.call ("plotrglplist", c(alist(plist, lighting = light, new = FALSE, 
    update = FALSE, scale = FALSE, smooth = smooth), dots))

  decorate3d(axes = FALSE, box = TRUE, zlab = "")

 # zoom and centering restored
  par3d(zoom = 1)
  uM <- par3d("userMatrix")
  uM[ ,4] <- c(0,0,0,1)
  par3d(userMatrix = uM)

  if (plist$type == "2D") {
    axis3d("x")
    axis3d("y") 
  }
  invisible(plist)
}

croprgl <- function(xlim = NULL, ylim = NULL, zlim = NULL, ...) {
  if (is.null(xlim) & is.null(ylim) & is.null(zlim))
    plist <- getplist()
  
  else 
    plist <- clipplist(xlim, ylim, zlim)

  rgl.clear()
  dots <- list(...)
  light <- dots$lighting
  if (is.null(light))
    light <- getplist()$rgl$lighting
  smooth <- dots$smooth
  if (is.null(smooth))
    smooth <- getplist()$rgl$smooth
  dots$lighting <- dots$new <- dots$smooth <- NULL
  do.call ("plotrglplist", c(alist(plist, lighting = light, new = FALSE, 
    update = FALSE, scale = FALSE, smooth = smooth), dots))

  decorate3d(axes = FALSE, box = TRUE, zlab = "")

 # zoom and centering restored
  par3d(zoom = 1)
  uM <- par3d("userMatrix")
  uM[ ,4] <- c(0,0,0,1)
  par3d(userMatrix = uM)

  if (plist$type == "2D") {
    axis3d("x")
    axis3d("y") 
  }
  invisible(plist)
}

uncutrgl <- function(...) {
  rgl.clear()
  dots <- list(...)
  light <- dots$lighting
  if (is.null(light))
    light <- getplist()$rgl$lighting
  if (is.null(light)) 
    light <- FALSE
  smooth <- dots$smooth
  if (is.null(smooth))
    smooth <- getplist()$rgl$smooth
  dots$lighting <- dots$new <- dots$smooth <- NULL

  do.call ("plotrgl", c(alist(lighting = light, new = FALSE, smooth = smooth), dots))
  par3d(zoom = 1)

  plist <- getplist()
#  par3d(userMatrix = plist$rgl$userMatrix)
  if (! is.null(plist$type))
    if (getplist()$type == "2D") {
      axis3d("x")
      axis3d("y") 
    }
  invisible(plist)  
}

uncroprgl <- function(...) 
  invisible(uncutrgl(...))


