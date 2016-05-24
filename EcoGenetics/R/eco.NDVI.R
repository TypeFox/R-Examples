#' Generating atmospherically corrected NDVI and MSAVI2 images for 
#' temporal series of Landsat 5 and 7
#' 
#' @description This program generates atmospherically corrected images of
#' NDVI and MSAVI2. The images of multiple dates can be processed in a single
#' run. The bands 4 and 3 of each date (previously subsetted to
#' the region of analysis) must be in the working directory. A table with information
#' for each image as described in the parameter tab (see also the example) 
#' is needed for processing the information.
#' 
#' @param tab data.frame with 7 columns: The date of the images
#' (format: YYYY/MM/DD), the sun elevation (both values could be extracted
#' from Landsat headers), the name of the band 4, the name of the band 3,
#' the starting haze value of the band 4, the starting haze value of
#' the band 3, and the name of the output file. Each row corresponds to
#' an image of different date.
#' @param correct Correction method ("COST", "DOS").
#' @param method Vegetation index ("NDVI", "MSAVI2").
#' @param landsat Satellite data source ("LT5" for Landsat 5, 
#' "LT7.L" for Landsat 7 low gain and "LT7.H" for Landsat 7 high gain).
#' @param datatype type of data, see \code{\link[raster]{dataType}}. 
#' Default "FLT4S".
#' 
#' @examples
#' \dontrun{
#' require(raster)
#' 
#' data(tab)
#' 
#' temp <-list()
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
#' writeRaster(temp[[1]], "20040719b4.tif", overwrite=T)
#' writeRaster(temp[[2]], "20040719b3.tif", overwrite=T)
#' writeRaster(temp[[3]], "20091106b4.tif", overwrite=T)
#' writeRaster(temp[[4]], "20091106b3.tif", overwrite=T)
#'
#' # Computing NDVI images: 
#' 
#' eco.NDVI(tab, "COST", "NDVI", "LT5")
#' 
#' example <- raster("NDVICOST20040719.tif")
#' image(example)
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

setGeneric("eco.NDVI", 
           function(tab, 
                    correct = c("COST", "DOS"), 
                    method = c("NDVI", "MSAVI2"), 
                    landsat = c("LT5", "LT7.L", "LT7.H"),
                    datatype =  c("FLT4S", "FLT8S", "INT4U", "INT4S", 
                                  "INT2U", "INT2S", "INT1U", "INT1S", 
                                  "LOG1S")) {
             
             correct <- match.arg(correct)
             method <- match.arg(method)
             landsat <- match.arg(landsat)
             datatype <- match.arg(datatype)
             
             
             steps <- nrow(tab)
             
             if(landsat == "LT5") {
               bresb3 <- -2.213976
               gresb3 <- 1.043976
               bresb4 <- -2.386024
               gresb4 <- 0.8760236
               esunb4 <- 1031
               esunb3 <- 1536
               
             } else if(landsat == "LT7.L") {
               
               bresb3 <- -5.9425197
               gresb3 <- 0.9425197
               bresb4 <- -6.0692913
               gresb4 <- 0.9692913
               esunb4 <- 1039
               esunb3 <- 1533
               
             } else if(landsat == "LT7.H") {
               
               bresb3 <- -5.6216535
               gresb3 <- 0.6216535
               bresb4 <- -5.7397638
               gresb4 <- 0.6397638
               esunb4 <- 1039
               esunb3 <- 1533
               
             }
             
             radiancia <- function(banda, brescale, grescale) {
               banda[] <- grescale * banda[] + brescale
               banda
             }
             
             reflectancia <- function(banda, edist, Esun, coselev, TAUz) {
               banda[] <-  (pi * edist ^ 2 * banda[]) / (Esun * coselev * TAUz)
               banda
             }
             
             Lhaze <- function(Lhazeban, grescale, brescale, Esun, coselev,
                               edist, TAUz) {
               Lhazeban <-  (grescale * Lhazeban + brescale) - 
                 0.01 * (Esun * coselev * TAUz) / (pi * edist ^ 2)
               Lhazeban
             }
             
             
             
             esperar <- function(i) {
               cat ("\r", ceiling(100 * i / steps), "% ", "completed", sep = "")
             }
             
             
             ESdist <- function (adate) {
               edist <- julian(as.Date(adate),
                               origin = as.Date(paste(substring(adate, 1, 4),
                                                      "12", "31", sep = "-")))[[1]]
               edist <- 1 - 0.016729 * cos((2 * pi) * (0.9856 * (edist - 4) / 360))
               edist
             }
             
             cat("0% completed")
             
             j <- 0
             
             for(i in 1:steps) {
               sunelev <- tab[i, 2]
               edist <- ESdist(as.character(tab[i, 1]))
               suntheta <- (90 - sunelev) * pi / 180           
               suntheta <- cos(suntheta)
               if(correct == "DOS") {
                 TAUz <- 1} else {
                   TAUz <- suntheta
                 }
               
               b3 <- raster::raster(as.character(tab[i, 4]))
               
               
               message("band 3 loaded")
               
               Lhazeb3 <- tab[i, 6]
               
               b3 <- radiancia(b3, bresb3, gresb3)
               Lhazeb3 <- Lhaze(Lhazeb3, gresb3, bresb3, esunb3, suntheta, edist, TAUz)
               b3[] <- b3[] - Lhazeb3
               b3 <- reflectancia(b3, edist, esunb3, suntheta, TAUz)
               b3[b3 < 0] <- 0
               
               j <- j + steps / 4
               esperar(j)
               
               
               b4 <- raster::raster(as.character(tab[i, 3]))
               
               
               message("band 4 loaded")
               Lhazeb4 <- tab[i, 5]
               b4 <- radiancia(b4, bresb4, gresb4)
               Lhazeb4 <- Lhaze(Lhazeb4, gresb4, bresb4, esunb4, suntheta, edist, TAUz)
               b4[] <- b4[] - Lhazeb4
               b4 <- reflectancia(b4, edist, esunb4, suntheta, TAUz)
               b4[b4 < 0] <- 0
               
               j <- j + steps / 4
               esperar(j) 
               
               
               if(raster::extent(b4) != raster::extent(b3)) {  
                 dimension <- raster::intersect(b4, b3)
                 b4 <- crop(b4, dimension)
                 b3 <- crop(b3, dimension)
               }
               
               
               if (method  ==  "NDVI") {
                 NDVI = (b4 - b3) / (b4 + b3)
                 NDVI[NDVI <- 1] <- NA
                 #NDVI=(NDVI+1)*255/2    #OPTIONAL FOR RESCALING TO 8 BITS
                 
                 message("writing NDVI image on disk")
                 
                 raster::writeRaster(NDVI, as.character(paste("NDVI", correct, tab[i, 7],
                                                              sep = "")), format = "GTiff",
                                     datatype = datatype, overwrite = T)
                 
               } else if (method  ==  "MSAVI2") {
                 MSAVI2 = ((2 *b4 + 1) - sqrt((2 * b4 + 1) ^ 2 - 8 * (b4 - b3))) / 2
                 MSAVI2[MSAVI2 <- 1] <- NA
                 #MSAVI2 <- (MSAVI2+1)*255/2       #OPTIONAL FOR RESCALING TO 8 BITS
                 message("writing MSAVI2 image on disk")
                 
                 raster::writeRaster(MSAVI2, as.character(paste("MSAVI2", correct,
                                                                tab[i, 7], sep="")),
                                     format="GTiff", datatype = datatype, overwrite =T)
               }
               
             }
             cat("done!")
           })
