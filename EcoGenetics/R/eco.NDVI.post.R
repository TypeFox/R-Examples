#' Postprocessing for NDVI and MSAVI 2 temporal series of Landsat 5 and 7
#' 
#' @description 
#' This program must be used sequentially after \code{\link{eco.NDVI}}.
#' The inputs required (tab, correct, method) are the same described
#' and used in that function. The algorithm stacks the images and save
#' the stack into the working directory with the name "time.tif". 
#' If the user wishes, the program can also compute images of max, min,
#' mean and var by pixel over the temporal sequence. Default is "mean".
#' 
#' @param tab Table used with \code{\link{eco.NDVI}}.
#' @param correct Correction method used in \code{\link{eco.NDVI}}.
#' @param method The vegetation index used in \code{\link{eco.NDVI}}.
#' @param datatype Type of data, see \code{\link[raster]{dataType}}. 
#' Default "FLT4S".
#' @param what Functions to apply over the created stack. 
#' The values permitted are: "none", "max", "min", "mean" and "var".
#' The functions are implemented with \code{\link[raster]{calc}}.
#' If more that one function would be applied, must be used the following
#' syntax: c("fun_1", "fun:2", "fun_i").
#' 
#' @seealso eco.NDVI
#' @seealso extract
#' 
#' @examples
#' \dontrun{
#' require(raster)
#' 
#' data(tab)
#' data(eco3)
#' temp <- list()
#'
#' # we create 4 simulated rasters for the data included in the object tab:
#'
#' for(i in 1:4) {
#' temp[[i]] <- runif(19800, 0, 254)
#' temp[[i]] <- matrix(temp[[i]], 180, 110)
#' temp[[i]] <- raster(temp[[i]], crs="+proj=utm")
#' extent(temp[[i]])<-c(3770000, 3950000, 6810000, 6920000)
#' }
#'
#' writeRaster(temp[[1]], "20040719b4.tif", overwrite = T) 
#' writeRaster(temp[[2]], "20040719b3.tif", overwrite = T)
#' writeRaster(temp[[3]], "20091106b4.tif", overwrite = T)
#' writeRaster(temp[[4]], "20091106b3.tif", overwrite = T)
#'
#' # Computing NDVI images: 
#'
#' eco.NDVI(tab, "COST", "NDVI", "LT5")
#'
#' # Mean NDVI image computed over the NDVI images that we calculated:
#'
#' eco.NDVI.post(tab, "COST", "NDVI", what = c("mean", "var"))
#' mean.ndvi <- raster("NDVI.COST.mean.tif")
#' plot(mean.ndvi)
#'
#' # Extraction of the mean NDVI for each point in the object eco and plot
#' # of the data:
#' 
#' ndvi <- extract(mean.ndvi, eco3[["XY"]])
#' ndvi<- aue.rescale(ndvi)
#' plot(eco3[["XY"]][, 1], eco3[["XY"]][, 2], col=rgb(ndvi, 0, 0),
#' pch=15, main = "Mean NDVI", xlab = "X", ylab  = "Y")
#' }
#' 
#' @references 
#' Chander G., B. Markham, and D. Helder. 2009. Summary of current radiometric calibration 
#' coefficients for Landsat MSS, TM, ETM+, and EO-1 ALI sensors. Remote sensing of 
#' environment, 113: 893-903. 
#' 
#' Chavez P. 1989. Radiometric calibration of Landsat Thematic Mapper multispectral images. 
#' Photogrammetric Engineering and Remote Sensing, 55: 1285-1294.
#' 
#' Chavez P. 1996. Image-based atmospheric corrections-revisited and improved. 
#' Photogrammetric engineering and remote sensing, 62: 1025-1035. 
#' 
#' Goslee S. 2011. Analyzing remote sensing data in R: the landsat package. 
#' Journal of Statistical Software, 43: 1-25.
#' 
#' Song C., C. Woodcock, K. Seto, M. Lenney and S. Macomber. 2001. 
#' Classification and change detection using Landsat TM data: when and how to correct 
#' atmospheric effects?. Remote sensing of Environment, 75: 230-244.
#' 
#' Tucker C. 1979. Red and photographic infrared linear combinations for monitoring 
#' vegetation. Remote sensing of Environment, 8: 127-150. 
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @export


setGeneric("eco.NDVI.post", 
           function(tab,correct = c("COST", "DOS"), 
                    method = c("NDVI", "MSAVI2"), 
                    datatype =  c("FLT4S", "FLT8S", "INT4U", "INT4S", 
                                  "INT2U", "INT2S", "INT1U", "INT1S", 
                                  "LOG1S"), 
                    what = c("mean", "max", "min", "var", "none")) {
             
             correct <- match.arg(correct)														
             method <- match.arg(method)
             datatype <- match.arg(datatype)
             what <- match.arg(what, several.ok =  TRUE)
             
             y <- pmatch(what, c("mean", "max", "min", "var", "none"))
             
             
             esperar <- function(i){
               cat ("\r", 100 * i / steps, "% ", "complete",  sep = "")
             }
             
             
             steps <- nrow(tab)
             
             message("\n", "loading parameters", "\n")
             
             xmin <- rep(0, nrow(tab))
             xmax <- rep(0, nrow(tab))
             ymin <- rep(0, nrow(tab))
             ymax <- rep(0, nrow(tab))
             
             for(i in 1:nrow(tab)) {
               xmin[i] <- raster::extent(raster::raster(as.character(paste(method, correct,
                                                                           tab[i, 7], ".tif",
                                                                           sep = ""))))@xmin
               xmax[i] <- raster::extent(raster::raster(as.character(paste(method, correct,
                                                                           tab[i, 7], ".tif",
                                                                           sep = ""))))@xmax
               ymin[i] <- raster::extent(raster::raster(as.character(paste(method, correct,
                                                                           tab[i, 7], ".tif",
                                                                           sep = ""))))@ymin
               ymax[i] <- raster::extent(raster::raster(as.character(paste(method, correct,
                                                                           tab[i, 7], ".tif",
                                                                           sep = ""))))@ymax
             }
             
             dimension = c(max(xmin), min(xmax), max(ymin), min(ymax))
             
             message("\n", "stacking images, please wait...", "\n")
             
             raster::writeRaster(raster::crop(raster::raster(as.character(paste(method, correct,
                                                                                tab[1, 7], ".tif",
                                                                                sep = ""))), 
                                              raster::extent(dimension)), paste(method, correct, "time.tif",
                                                                                sep = ""), format = "GTiff",
                                 datatype = datatype, overwrite = T)
             
             
             for(i in 2:nrow(tab)) {
               raster::writeRaster(raster::addLayer(raster::brick(paste(method, correct, "time.tif",
                                                                        sep = "")), 
                                                    raster::crop(raster::raster(as.character(paste(method, correct,
                                                                                                   tab[i, 7], 
                                                                                                   ".tif", sep = ""))),
                                                                 raster::extent(dimension))), format = "GTiff",
                                   datatype = datatype,  "temporal.tif", 
                                   overwrite = T)    
               raster::writeRaster(raster::brick("temporal.tif"), paste(method, correct, "time.tif", 
                                                                        sep = ""), format = "GTiff",
                                   datatype = datatype,  overwrite = T)
               esperar(i)
               
             }
             file.remove("temporal.tif")
             
             
             message("computing...")
             
             time <- raster::brick(paste(method, correct, "time.tif",  sep = ""))
             
             if(any(y == 1)) {
               cat("\n", "computing mean", "\n")
               temp <- raster::calc(time, mean)
               raster::writeRaster(temp, paste(method, ".", correct, ".", "mean.tif",
                                               sep = ""), format="GTiff", 
                                   datatype = datatype,  overwrite = T)
             }
             
             if(any(y == 2)) {
               cat("\n", "computing max", "\n")
               temp <- raster::calc(time, max)
               raster::writeRaster(temp, paste(method, ".", correct, ".", "max.tif",
                                               sep = ""), format="GTiff",
                                   datatype = datatype, overwrite = T)
             }
             
             if(any(y == 3)) {
               cat("\n", "computing min", "\n") 
               temp <- raster::calc(time, min)
               raster::writeRaster(temp, paste(method, ".", correct, ".", "min.tif",
                                               sep = ""), format="GTiff", 
                                   datatype = datatype, overwrite = T)
             }
             
             
             if(any(y == 4)) {
               cat("\n", "computing variance", "\n")
               temp <- raster::calc(time, var)
               raster::writeRaster(temp, paste(method, ".", correct, ".", "var.tif",
                                               sep = ""), format = "GTiff", 
                                   datatype = datatype, overwrite = T)
             }
             
             cat("\n", "done", "\n")

           })
