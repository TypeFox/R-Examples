panel.hexloess <-
function(bin, w=NULL, span = 2/3, degree = 1, family = c("symmetric",
         "gaussian"), evaluation = 50, lwd = add.line$lwd, lty = add.line$lty,
         col, col.line = add.line$col, ...)
{
    stop("panel.hexloess is no longer available")
    add.line <- trellis.par.get("add.line")

##     x <- bin@xcm
##     y <- bin@ycm
##     if(is.null(w))w <- bin@count
##     control <- loess.control(...)
##     notna <- !(is.na(x) | is.na(y))
##     new.x <- seq(min(x[notna]), max(x[notna]), length = evaluation)
##     family <- match.arg(family)
##     iterations <- if (family == "gaussian") 1 else control$iterations
##     fit <- stats:::simpleLoess(y, x, w, span, degree, FALSE, FALSE,
##                                normalize = FALSE, "none", "interpolate",
##                                control$cell, iterations, control$trace.hat)
##     kd <- fit$kd
##     z <- .C("loess_ifit", as.integer(kd$parameter), as.integer(kd$a),
##         as.double(kd$xi), as.double(kd$vert), as.double(kd$vval),
##         as.integer(evaluation), as.double(x), fit = double(evaluation),
##         PACKAGE = "stats")$fit
##     if (length(x) > 0) {
##          if (!missing(col) && missing(col.line)) {
## 			 col.line <- col
##          }
##       add.line <- trellis.par.get("add.line")
##       panel.lines(new.x, z, col = col.line, lty = lty, lwd = lwd)
##     }
}

panel.hexgrid <- function(h, border=grey(.85))
{
  hexGraphPaper(h,border=border)
}
