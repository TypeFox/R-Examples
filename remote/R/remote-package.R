#' R EMpirical Orthogonal TEleconnections
#' 
#' A collection of functions to facilitate empirical orthogonal teleconnection analysis. 
#' Some handy functions for preprocessing, such as deseasoning, denoising, lagging 
#' are readily available for ease of usage.
#' 
#' @name remote-package
#' @aliases remote
#' @docType package
#' @title R EMpirical Orthogonal TEleconnections
#' @author Tim Appelhans, Florian Detsch\cr
#' \cr
#' \emph{Maintainer:} Tim Appelhans \email{tim.appelhans@@gmail.com}
#' 
#' @keywords package
#' @references 
#' Empirical Orthogonal Teleconnections\cr
#' H. M. van den Dool, S. Saha, A. Johansson (2000)\cr
#' Journal of Climate, Volume 13, Issue 8 (April 2000) pp. 1421 - 1435\cr
#' 
#' Empirical methods in short-term climate prediction\cr
#' H. M. van den Dool (2007)\cr
#' Oxford University Press, Oxford, New York (2007)\cr
#' @seealso \pkg{remote} is built upon Raster* classes from the \code{\link{raster-package}}.
#' Please see their documentation for data preparation etc.
#' @import Rcpp raster foreach gridExtra latticeExtra mapdata scales
#' 
NULL
#' 
#' @docType data 
#' @name vdendool
#' @title Mean seasonal (DJF) 700 mb geopotential heights
#' @description NCEP/NCAR reanalysis data of mean seasonal (DJF) 700 mb geopotential heights from 1948 to 1998
#' @details NCEP/NCAR reanalysis data of mean seasonal (DJF) 700 mb geopotential heights from 1948 to 1998
#' @format a RasterBrick with the following attributes\cr
#' \cr
#' dimensions  : 14, 36, 504, 50  (nrow, ncol, ncell, nlayers)\cr
#' resolution  : 10, 4.931507  (x, y)\cr
#' extent      : -180, 180, 20.9589, 90  (xmin, xmax, ymin, ymax)\cr
#' coord. ref. : +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0\cr 
#' @references
#' The NCEP/NCAR 40-year reanalysis project\cr
#' Kalnay et al. (1996)\cr
#' Bulletin of the American Meteorological Society, Volume 77, Issue 3, pp 437 - 471\cr
#' \url{http://journals.ametsoc.org/doi/abs/10.1175/1520-0477(1996)077%3C0437%3ATNYRP%3E2.0.CO%3B2
#' }
#' @source
#' \url{http://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.derived.pressure.html}\cr
#' \emph{Original Source:} NOAA National Center for Environmental Prediction
NULL
#' 
#' @docType data 
#' @name australiaGPCP
#' @title Monthly GPCP precipitation data for Australia
#' @description Monthly Gridded Precipitation Climatology Project precipitation data 
#' for Australia from 1982/01 to 2010/12
#' @details Monthly Gridded Precipitation Climatology Project precipitation data 
#' for Australia from 1982/01 to 2010/12
#' @format a RasterBrick with the following attributes\cr
#' \cr
#' dimensions  : 12, 20, 240, 348  (nrow, ncol, ncell, nlayers)\cr
#' resolution  : 2.5, 2.5  (x, y)\cr
#' extent      : 110, 160, -40, -10  (xmin, xmax, ymin, ymax)\cr
#' coord. ref. : +proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs\cr 
#' @references
#' The Version-2 Global Precipitation Climatology Project (GPCP) Monthly Precipitation Analysis (1979 - Present)\cr
#' Adler et al. (2003)\cr
#' Journal of Hydrometeorology, Volume 4, Issue 6, pp. 1147 - 1167\cr
#' \url{http://dx.doi.org/10.1175/1525-7541(2003)004<1147:TVGPCP>2.0.CO;2}
NULL
#' 
#' @docType data 
#' @name pacificSST
#' @title Monthly SSTs for the tropical Pacific Ocean
#' @description Monthly NOAA sea surface temperatures for the tropical Pacific Ocean from 1982/01 to 2010/12
#' @details Monthly NOAA sea surface temperatures for the tropical Pacific Ocean from 1982/01 to 2010/12
#' @format a RasterBrick with the following attributes\cr
#' \cr
#' dimensions  : 30, 140, 4200, 348  (nrow, ncol, ncell, nlayers)\cr
#' resolution  : 1, 1  (x, y)\cr
#' extent      : 150, 290, -15, 15  (xmin, xmax, ymin, ymax)\cr
#' coord. ref. : +proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs\cr 
#' @references
#' Daily High-Resolution-Blended Analyses for Sea Surface Temperature\cr
#' Reynolds et al. (2007)\cr
#' Journal of Climate, Volume 20, Issue 22, pp. 5473 - 5496\cr
#' \url{http://dx.doi.org/10.1175/2007JCLI1824.1}
NULL
