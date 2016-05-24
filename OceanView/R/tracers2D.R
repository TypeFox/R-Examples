## CHECK WITH QUIVER
## =============================================================================
## =============================================================================
## Tracers in 2D
## =============================================================================
## =============================================================================

tracers2D <- function(x, y, colvar = NULL, ..., 
                    col = NULL, NAcol = "white", colkey = NULL, 
                    mask = NULL, image = FALSE, contour = FALSE, 
                    clim = NULL, clab = NULL)  {
# ------------------------------------------------------------------------------
# check input
# ------------------------------------------------------------------------------

  dots  <- splitpardots( list(...) )
  dp    <- dots$points
  dm    <- dots$main 

  image   <- check.args(image, NULL)
  contour <- check.args(contour, NULL)
  add <- FALSE 

 # images, contours or masked cells
  if (! is.null(mask)) {
    if (image$add)
      stop ("cannot have both 'image' and 'mask' specified")

    maskNAcol <- mask$NAcol
    if (is.null(maskNAcol)) 
      maskNAcol <- "black" 
    X <- mask$x 
    Y <- mask$y 
    Z <- mask$z 
    Z[!is.na(Z)] <- 0
    Z[is.na(Z)]  <- 1
    do.call ("image2D", c(alist(z = Z, x = X, y = Y, colkey = FALSE,
               col = c("white", maskNAcol)), dm))
    add <- TRUE
  }

  if (image$add) {
    if (is.null(image$args$z))
      stop("image$'z' should be specified if image is a list")
    if (is.null(image$args$col))
      image$args$col <- jet.col(100)
    if (!add) 
      image$args <- c(image$args, dm)
    
    do.call("image2D", c(alist(colkey = image$colkey), image$args))   
    add <- TRUE
  }

  if (contour$add) {
    if (image$add | add) 
      contour$args$add <- TRUE
    else
      contour$args$add <- FALSE
      
    if (contour$args$add) 
      contour$args$colkey <- FALSE

    if (is.null(contour$args$z))
      stop("contour$'z' should be specified if 'contour' is a list")

    if (!contour$args$add) 
      contour$args <- c(contour$args, dm)

    do.call("contour", contour$args)
    add <- TRUE
  }

  if (! add) {
    if (is.null(dm$xlab)) 
      dm$xlab <- "x"
    if (is.null(dm$ylab)) 
      dm$ylab <- "y"

    do.call("plot", c(alist(x, y, type = "n"), dm))
  }
  
  do.call("scatter2D", c(alist(x, y, colvar = colvar, 
     col = col, NAcol = NAcol, clim = clim,
     add = add, colkey = colkey),  dp))
  
}
