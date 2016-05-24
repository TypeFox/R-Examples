"gplot.matrix" <-
  function(x, y, se.y = NULL, xlab, ylab, residuals = NULL, rugplot = FALSE, scale = 
           0, se = FALSE, fit, ...)
{
  if(ncol(x) != 2) {
    warning(paste("A perspective plot was requested for \"", ylab,
                  "\" but the \"x\" variable has dimension other than 2",
                  sep = ""))
    invisible(return(0))
  }
  bivar.dup <- function(x)
    {
      if(is.null(dx <- dim(x)) || dx[2] > 2)
        stop("x must be bivariate")
      duplicated(x[, 1] + (1i) * x[, 2])
    }


#  interp.loaded<-require("akima")
  interp.loaded<-TRUE
#  if(!interp.loaded)
  xname <- dimnames(x)[[2]]
  dups <- bivar.dup(x)
  if (requireNamespace("akima", quietly = TRUE)) {
    xyz <- akima::interp(x[!dups, 1], x[!dups, 2], y[!dups])
   } else {
    stop("You need to install the package 'akima' from the R contributed libraries to use this plotting method for bivariate functions")
   }

  zmin <- min(xyz$z[!is.na(xyz$z)])
  z <- ifelse(is.na(xyz$z), zmin, xyz$z)
  scale2 <- diff(range(z))
                                        # Adjust scale
  scale <- max(scale, scale2)
                                        #	persp(xyz$x, xyz$y, (z - zmin)/scale, xlab = xname[1], ylab = xname[
                                        #		2], zlab = ylab, ...)
  persp(xyz$x, xyz$y, z, xlab = xname[1], ylab = xname[2], zlab = ylab,
        ...)
  invisible(scale)
}
