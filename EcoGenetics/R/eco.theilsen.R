#' Theil-sen regression for a raster time series
#' 
#' @description This function computes the theil-sen estimator and 
#' the P-value associated for each pixel over time in a stack of images,
#' writing the values in a raster (one for the estimators and one for 
#' the P-values). It is recommended to use a "RasterBrick", that
#' is more efficient in managing memory.
#' 
#' @param stacked Stacked images ("RasterLayer"  or "RasterBrick").
#' @param date data vector with dates for each image.
#' @param adjust P-values correction method for multiple tests 
#' passed to \code{\link[stats]{p.adjust}}. Defalut is "none".
#' 
#' @seealso \code{\link[rkt]{rkt}}.
#' 
#' @examples
#' \dontrun{
#' require("raster")
#' set.seed(6)
#' 
#' temp <- list()
#' for(i in 1:100) {
#' temp[[i]] <- runif(36,-1, 1)
#' temp[[i]] <- matrix(temp[[i]], 6, 6)
#' temp[[i]] <- raster(temp[[i]])
#'}
#'
#'temp <- brick(temp)
#'
#'
#'writeRaster(temp,"temporal.tif", overwrite=T)
#'rm(temp)
#'ndvisim <- brick("temporal.tif")
#'
#'date <- seq(from = 1990.1, length.out = 100, by = 0.2)
#'
#'eco.theilsen(ndvisim, date)
#'
#'pvalue <- raster("pvalue.tif")
#'slope <- raster("slope.tif")
#'par(mfrow = c(1, 2))
#'plot(pvalue, main = "p-value")
#'plot(slope, main = "slope")
#'}
#'
#' @references 
#' Sen, P. 1968. Estimates of the regression coefficient based on Kendall's tau. 
#' Journal of the American Statistical Association, Taylor and Francis Group, 63: 1379-1389.
#' 
#' Theil H. 1950. A rank-invariant method of linear and polynomial regression analysis, 
#' Part 3 Proceedings of Koninalijke Nederlandse Akademie van Weinenschatpen A, 53: 397-1412.
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @export

setGeneric("eco.theilsen", 
           function(stacked, date, 
                    adjust = "none") {
                    
             adjust <- match.arg(adjust)
             
             esperar <- function(i) {
               cat ("\r", ceiling(100 * i / steps), "% ",
                    "completed", sep = "")
             }
              
             cat("starting...", "\n\n")
             
             fun <- function(date, data)  {
               mod <- rkt::rkt(date, data)
               return(c(as.numeric(mod[3]), as.numeric(mod[1])))
             }
             
             cat("pre-processing data...", "\n\n")
             
             pendiente <- rep(NA, raster::ncell(stacked))
             pvalor <- rep(NA, raster::ncell(stacked))
             df <- raster::as.matrix(stacked)
             steps = nrow(df)
             
             for(i in 1:nrow(df)) {
               resultados <- fun(date,as.numeric(df[i, ]))
               pendiente[i] <- resultados[1]
               pvalor[i] <- resultados[2]
               esperar(i)
             }
             cat("\n\n")
             
             
             if(adjust != "none") {
               cat(paste("adjusting p values with", adjust, "method"), "\n")
               pvalor <- p.adjust(pvalor, "adjust")
             }
             
             cat("writing slope image to workspace...", "\n\n")
             
             pendiente <- matrix(pendiente, nrow = stacked@nrows,
                                 ncol = stacked@ncols, byrow = TRUE)
             pendiente <- raster::raster(pendiente, crs = stacked@crs,
                                         xmn = stacked@extent@xmin,
                                         ymn = stacked@extent@ymin, 
                                         xmx = stacked@extent@xmax,
                                         ymx = stacked@extent@ymax)
             raster::writeRaster(pendiente, "slope.tif", overwrite = T)
             
             cat("writing P-value image to workspace...", "\n")
             
             pvalor <- matrix(pvalor, nrow=stacked@nrows, ncol = stacked@ncols,
                              byrow = TRUE)
             pvalor <- raster::raster(pvalor, crs = stacked@crs,
                                      xmn = stacked@extent@xmin,
                                      ymn = stacked@extent@ymin, 
                                      xmx = stacked@extent@xmax,
                                      ymx = stacked@extent@ymax)
             
             cat("\n","done!","\n\n" )
             
             
             raster::writeRaster(pvalor, "pvalue.tif", overwrite = T)
             
             
           })
