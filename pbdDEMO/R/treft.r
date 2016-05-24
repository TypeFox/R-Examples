#' Surface Air Temperature at Reference Height (TREFHT)
#' 
#' This is a practical example in NetCDF4 format and for data reading, writing,
#' and transforming. This dataset is a partial output of the Surface Air
#' Temperature at Reference Height (TREFHT) which is monthly averaged of Jan.
#' 2004 from a CAM5 simulation.  This dataset only contains a tiny part of
#' ultra-large simulations conducted by Mr Prabhat and Michael Wehner of
#' Lawrence Berkeley National Laboratory.
#' 
#' Version 5.0 of the Community Atmosphere Model (CAM) is the latest in a
#' series of global atmosphere models developed primarily at the National
#' Center for Atmospheric Research (NCAR).
#' 
#' \code{TREFHT} contains two lists: \code{def} and \code{data}.
#' 
#' \code{def} is a list contains usual definitions of NetCDF4. In this case,
#' they define the variable \dQuote{TREFHT} including 2D dimensions 1152
#' longitudes and 768 latitudes, 1 time step, the unit in Kelvin, \dots{} etc.
#' 
#' \code{data} contains values in matrix with dimension \eqn{1152\times
#' 768}{1152*768}. Note that this matrix stores data in C format (column
#' major), so it needs a transpose to obtains the R/Fortran format (row major).
#' Also, the longitude order is not the same as the \pkg{maps} package. Please
#' see the example below for the adjustment or by calling
#' \code{demo('trefht','pbdDEMO')} inside an R session.
#' 
#' @name Temperature at Reference Height
#' @aliases Practical Example TREFHT
#' @docType data
#' @format An R data file contains two lists: \code{def} for structure
#' definition of \dQuote{TREFHT} in \code{ncvar4} class (see \pkg{pbdNCDF4}
#' package for details), and \code{data} for output values of simulation in a
#' matrix where rows are for 1152 longitudes and columns are for 768 latitudes.
#' @author Mr Prabhat and Michael Wehner.
#' @seealso \code{ncvar_put_2D} and \code{ncvar_get_2D}.
#' @references More datasets are available on ESGF
#' (\url{http://www.earthsystemgrid.org/}) through the C20C project (on the
#' NERSC portal).
#' 
#' CAM5: \url{http://www.cesm.ucar.edu/models/cesm1.0/cam/}
#' 
#' Programming with Big Data in R Website: \url{http://r-pbd.org/}
#' @keywords datasets
#' @examples
#'  \dontrun{ library(maps)
#' library(RColorBrewer)
#' library(pbdDEMO, quiet = TRUE)
#' 
#' lon <- TREFHT$def$dim[[1]]$vals               # longitude
#' lat <- TREFHT$def$dim[[2]]$vals               # latitude
#' da <- TREFHT$data                             # surface temperature
#' 
#' # Define Axes.
#' x <- c(lon[lon > 180] -360, lon[lon <= 195])  # adjustment for maps
#' y <- lat
#' z <- rbind(da[lon > 180,], da[lon <= 195,])   # adjustment for maps
#' xlim <- range(x)
#' ylim <- range(y)
#' zlim <- range(z)
#' col.z <- c(colorRampPalette(c("#0000FF", "#2BFCD3"))(100),
#'            colorRampPalette(c("#2BFCD3", "#5300AB"))(100),
#'            colorRampPalette(c("#5300AB", "#7CFA82"))(100),
#'            colorRampPalette(c("#7CFA82", "#A90055"))(100),
#'            colorRampPalette(c("#A90055", "#D6FC28"))(100),
#'            colorRampPalette(c("#D6FC28", "#FE0001"))(100))
#' 
#' # Plot
#' layout(matrix(c(1, 2), ncol = 1), heights = c(2, 1))
#' par(mar = c(4, 4, 4, 0))
#' plot(NULL, NULL, xlim = xlim, ylim = ylim, type = "n", axes = FALSE,
#'      xlab = "Longitude", ylab = "Latitude", main = "TREFHT (Jan. 2004)")
#' image(x, y, z, zlim = zlim, xlim = xlim, ylim = ylim,
#'       col = col.z, add = TRUE)
#' 
#' # Add Map.
#' map(add = TRUE)
#' abline(h = c(-23.5, 0, 23.5), v = 0, lty = 2)
#' xtickets <- seq(-180, 180, by = 30)
#' ytickets <- seq(-90, 90, by = 30)
#' box()
#' axis(1, at = xtickets, labels = xtickets)
#' axis(2, at = ytickets, labels = ytickets)
#' 
#' # Add Legend.
#' z.temp <- matrix(seq(zlim[1], zlim[2], length = 500), ncol = 1)
#' ztickets <- seq(230, 300, by = 10)
#' par(mar = c(4, 4, 0, 1))
#' plot(NULL, NULL, xlim = zlim, ylim = c(0, 1), type = "n", axes = FALSE,
#'      xlab = "TREFHT (Kelvin)", ylab = "")
#' image(z.temp, 0, z.temp, zlim = zlim, xlim = zlim, ylim = c(0, 1),
#'       col = col.z, add = TRUE)
#' axis(1, at = ztickets, labels = ztickets)
#' }
#' 
NULL

