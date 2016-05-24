#' Plots for exploratory spatial data analysis (ESDA)
#' 
#' This function creates four plots for exploratory spatial data analysis 
#' (ESDA): histogram + density plot, bubble plot, variogram plot, and variogram 
#' map.
#' 
#' @param z Vector of numeric values of the variable for with ESDA plots should 
#' be created.
#' @param lat Vector of numeric values containing the y coordinate (latitude) 
#' of the point locations where the \code{z} variable was observed.
#' @param lon Vector of numeric values containing the x coordinate (longitude) 
#' of the point locations where the \code{z} variable was observed.
#' @param lags Numerical vector; upper boundaries of lag-distance classes. See
#' argument \code{boundaries} of \code{\link[gstat]{variogram}} for more info.
#' @param cutoff Integer value defining the spatial separation distance up to 
#' which point pairs are included in semi-variance estimates. Defaults to the 
#' length of the diagonal of the box spanning the data divided by three.
#' @param width Integer value specifying the width of subsequent distance 
#' intervals into which data point pairs are grouped for semi-variance 
#' estimates. Defaults to \code{width = cutoff / 20}.
#' 
#' @details
#' The user should visit the help pages of \code{\link[gstat]{variogram}},
#' \code{\link[pedometrics]{plotHD}}, \code{\link[sp]{bubble}} and
#' \code{\link[sp]{spplot}} to obtain more details about the main functions 
#' used to built \code{plotESDA}.
#' 
#' @return
#' Four plots: histogram and density plot, bubble plot, empirical variogram, 
#' and variogram map.
#' 
#' @references
#' Cressie, N.A.C. (1993) \emph{Statistics for Spatial Data}. New York: John 
#' Wiley \& Sons, p.900, 1993.
#' 
#' Pebesma, E.J. (2004) Multivariable geostatistics in S: the gstat package. 
#' \emph{Computers \& Geosciences}, 30:683-691, 2004.
#' 
#' Webster, R. \& Oliver, M.A. \emph{Geostatistics for environmental 
#' scientists}. Chichester: John Wiley \& Sons, p.315, 2007.
#' 
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' 
#' @seealso \code{\link[gstat]{variogram}}, \code{\link[pedometrics]{plotHD}},
#' \code{\link[sp]{bubble}}, \code{\link[sp]{spplot}}.
#' @importFrom graphics plot
#' @export
#' @examples
#' # require(gstat)
#' # data(meuse)
#' # plotESDA(z = meuse$zinc, lat = meuse$y, lon = meuse$x)
#' 
#' @keywords dplot
#' 
# FUNCTION #####################################################################
plotESDA <- 
  function (z, lat, lon, lags, cutoff, width = c(cutoff / 20)) {
    
    # Check if suggested packages are installed
    pkg <- c("gstat", "sp")
    id <- !sapply(pkg, requireNamespace, quietly = TRUE)
    if (any(id)) {
      pkg <- paste(pkg[which(id)], collapse = " ")
      stop(paste("Package(s) needed for this function to work but not",
                 "installed: ", pkg, sep = ""), call. = FALSE)
    }
    
    if (missing(z)) {
      stop("<z> is a mandatory argument")
    }
    if (missing(lon)) {
      stop("<lon> is a mandatory argument")
    }
    if (missing(lat)) {
      stop("<lat> is a mandatory argument")
    }
    if (!any(class(z) == c("numeric", "integer"))) {
      stop("<z> should be of class numeric or integer")
    }
    if (!any(class(lat) == c("numeric", "integer"))) {
      stop("<lat> should be of class numeric or integer")
    }
    if (!any(class(lon) == c("numeric", "integer"))) {
      stop("<lon> should be of class numeric or integer")
    }
    if (length(unique(c(length(z), length(lat), length(lon)))) > 1) {
      stop("<z>, <lat> and <lon> must have the same length")
    }
    db <- data.frame(lon = lon, lat = lat, z = z)
    sp::coordinates(db) <- ~ lon + lat
    
    # Estimate the cutoff
    cutoff <- max(gstat::variogram(z ~ 1, loc = db)$dist)
    
    # Bubble plot
    v1 <- sp::bubble(db, zcol = "z", fill = FALSE, main = "", maxsize = 1)
    
    # Variogram map
    v2 <- gstat::variogram(z ~ 1, loc = db, map = TRUE, cutoff = cutoff, 
                           width = width)
    v2 <- sp::spplot(v2$map[2], col.regions = sp::bpy.colors(64))
    
    # Sample variogram
    if (missing(lags)) {
      v3 <- gstat::variogram(z ~ 1, loc = db, cutoff = cutoff, width = width)
    } else {
      v3 <- gstat::variogram(z ~ 1, loc = db, boundaries = lags)
    }
    v3 <- plot(v3, cex = 0.5, type = "b", pch = 20, asp = 1)
    
    # Histogram
    v4 <- plotHD(z, HD = "over", stats = FALSE, asp = 1, xlab = "z",
                 col = c("skyblue", "red"))
    print(v4, split = c(1, 1, 2, 2), more = TRUE)
    print(v3, split = c(1, 2, 2, 2), more = TRUE)
    print(v1, split = c(2, 1, 2, 2), more = TRUE)
    print(v2, split = c(2, 2, 2, 2), more = FALSE)
  }
# End!
