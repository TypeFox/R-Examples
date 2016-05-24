plot2d.default <- function(x, residuals, range, col.residuals = "black",
  fill.select = NULL, col.polygons = NULL, col.rug = NULL, pb = FALSE, 
  x.co = NULL, rug = FALSE, jitter = FALSE, specs)
{
  if(residuals && !is.null(pres <- attr(x, "partial.resids")))
    residuals <- TRUE
  else
    residuals <- FALSE
  if(nrow(x) > 1)
    x <- na.omit(x)
  if(!is.matrix(x))
    x <- matrix(x, nrow = 1L)
  if(residuals)
    e <- attr(x, "partial.resids")
  x <- unique(x)
  if(pb) {
    nc <- ncol(x)
    if(length(ux <- unique(x[,2L:nc])) < 3L) {
      fill.select <- NULL
      if(!is.matrix(ux))
        ux <- matrix(ux, nrow = 1L)
    } else ux <- matrix(unique(x[,2L:nc]), ncol = (nc - 1L))
    nux <- nrow(ux)
    if(nux < 2L) {
      nux <- 2L
      ux <- rbind(ux, ux)
    }
    x.co <- seq(x.co + range[1L], x.co - range[2L], length = nux)
    x <- cbind(x.co, ux)
    x <- rbind(x, x, x)
  }
  x <- x[order(x[,1L]),]
  if(!is.null(fill.select)) {      
    ufs <- unique(fill.select)
    ufs <- ufs[ufs != 0]
    nu <- length(ufs)
    if(!is.null(specs$poly.lty))
      specs$poly.lty <- rep(specs$poly.lty, length.out = nu)
    else
      specs$poly.lty <- rep(0, nu)
    if(is.null(specs$angle))
      specs$angle <- rep(45, nu)
    else
      specs$angle <- rep(specs$angle, length.out = nu)
    if(!is.null(specs$density))
      specs$density <- rep(specs$density, length.out = nu)
    else
      specs$density <- NULL
    if(!is.null(specs$border))
      specs$border <- rep(specs$border, length.out = nu)
    if(!is.null(specs$poly.lwd))
      specs$poly.lwd <- rep(specs$poly.lwd, length.out = nu)
    else
      specs$poly.lwd <- rep(1, nu)
    for(k in 1L:nu) {
      check <- fill.select == ufs[k]
      if(length(check) == ncol(x)) {
        poly <- x[,check]
        p1 <- poly[,1L]
        p2 <- poly[,2L]
        y.co <- c(p1, p2[length(p2):1L])
        x.co <- x[,1L]
        x.co <- c(x.co, x.co[length(x.co):1L])
        graphics::polygon(x = x.co, y = y.co, col = col.polygons[k], 
          lty = specs$poly.lty[k], border = specs$border[k], 
          density = specs$density[k], angle = specs$angle[k], 
          lwd = specs$poly.lwd[k])
      }
    }
  }    
  if(residuals) {
    pargs <- list()
    pargs$x <- pres[,1L]
    pargs$y <- pres[,2L]
    pargs$cex <- specs$cex
    pargs$type <- specs$type
    pargs$pch <- specs$pch
    pargs$col <- col.residuals
    do.call(graphics::points, pargs)
  }
  for(k in 2L:ncol(x)) {
    lines(x[,k] ~ x[,1L], lty = specs$lty[k - 1L], lwd = specs$lwd[k - 1L], 
      col = specs$col.lines[k - 1L])
  }
  if(rug) {
    specs$col <- col.rug
    rugp <- if(!is.null(specs$rugp)) specs$rugp else x[,1L]
    if(jitter)      
      specs$x <- jitter(rugp)
    else
      specs$x <- rugp
    do.call(graphics::rug, delete.args(graphics::rug, specs))
  }

  return(invisible(NULL))
}

