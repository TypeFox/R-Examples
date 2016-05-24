
#' Create aoi polygon from topleft and bottomright coordinates
#' @description
#' An AOI (Area of Interest) is created based on two points (topleft and bottomright) using a coordinate reference system.
#' @param topleft     a vector with topleft x,y coordinates 
#' @param bottomright a vector with bottomright x,y coordinates
#' @param EPSG        Coordinate reference system EPSG code
#' @return object of class SpatialPolygons
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David
#' @examples 
#' tl <- c(493300, -3592700)
#' br <- c(557200, -3700000) 
#' aoi <- createAoi(topleft = tl, bottomright=br, EPSG=32619)
#' plot(aoi)
#' @import rgdal raster sp 
#' @export
# Maybe i can provide some points and use CHull like in QGIS-Geostat
createAoi <- function(topleft, bottomright, EPSG){
  aoi <- SpatialPolygons(
    list(Polygons(list(Polygon(coords = matrix(
      c(topleft[1],bottomright[1], bottomright[1],topleft[1],topleft[1],
        topleft[2], topleft[2], bottomright[2], 
        bottomright[2],topleft[2]), ncol=2, nrow= 5))), ID=1)))
  if(!missing(EPSG)){aoi@proj4string <- CRS(paste0("+init=epsg:", EPSG))}
  return(aoi)
}

#' Load Landsat data from folder
#' @description
#' This function loads Landsat bands from a specific folder. 
#' @param path  folder where band files are stored
#' @param sat   "L7" for Landsat 7, "L8" for Landsat 8 or "auto" to guess from filenames
#' @param aoi   area of interest to crop images, if waterOptions("autoAoi") == 
#' TRUE will look for any object called aoi on .GlobalEnv
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
loadImage <-  function(path = getwd(), sat="auto", aoi){
  if(sat=="auto"){sat = getSat(path)} #DRY!
  if(sat=="L8"){bands <- 2:7}
  if(sat=="L7"){bands <- c(1:5,7)}
  
  stack1 <- list()
  for(i in 1:6){
    stack1[i] <- raster(list.files(path=path, 
      pattern = paste0("^L[EC]\\d+\\w+\\d+_(B|band)", bands[i] ,".(TIF|tif)$")))
  }
  raw.image <- do.call(stack, stack1)
  raw.image <- aoiCrop(raw.image, aoi)                               
  raw.image <- saveLoadClean(imagestack = raw.image, 
                       stack.names = c("B", "G", "R", "NIR", "SWIR1", "SWIR2"), 
                       file = "imageDN", 
                       overwrite=TRUE)
  return(raw.image) 
}  

#' Calculates Top of atmosphere reflectance
#' @description
#' This function calculates the TOA (Top Of Atmosphere) reflectance considering only the image metadata.
#' @param image.DN      raw image in digital numbers
#' @param sat           "L7" for Landsat 7, "L8" for Landsat 8 or "auto" to guess from filenames 
#' @param ESPA          Logical. If TRUE will look for espa.usgs.gov related products on working folder
#' @param aoi           area of interest to crop images, if waterOptions("autoAoi") == TRUE will look for any object called aoi on .GlobalEnv
#' @param incidence.rel solar incidence angle, considering the relief
#' @param MTL           Landsat Metadata File
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#'
#' LPSO. (2004). Landsat 7 science data users handbook, Landsat Project Science Office, NASA Goddard Space Flight Center, Greenbelt, Md., (http://landsathandbook.gsfc.nasa.gov/) (Feb. 5, 2007) \cr
#' @export
calcTOAr <- function(image.DN, sat="auto", 
                     ESPA=FALSE, aoi, incidence.rel, MTL){
  path = getwd()
  if(sat=="auto"){sat = getSat(path)}
  if(sat=="L8"){bands <- 2:7}
  if(sat=="L7"){bands <- c(1:5,7)}
  if(ESPA==TRUE & sat=="L8"){
    files <- list.files(path = path, pattern = "_toa_band+[2-7].tif$", full.names = T)
    stack1 <- list()
    for(i in 1:6){
      stack1[i] <- raster(files[i])
    }
    image_TOA <- do.call(stack, stack1)
    image_TOA <- aoiCrop(image_TOA, aoi)
    image_TOA <- image_TOA / 10000
  }
  ### Ro TOA L7
  if(sat=="L7"){
    if(missing(MTL)){MTL <- list.files(path = path, pattern = "MTL.txt", full.names = T)}
    MTL <- readLines(MTL, warn=FALSE)
    ESUN <- c(1997, 1812, 1533, 1039, 230.8, 84.90) # Landsat 7 Handbook
    ## On sect 11.3, L7 handbook recommends using this formula for Ro TOA
    ## O using DN - QCALMIN with METRIC 2010 formula.
    Gain <- c(1.181, 1.210, 0.943, 0.969, 0.191, 0.066)
    Bias <- c(-7.38071, -7.60984, -5.94252, -6.06929, -1.19122, -0.41650)
    if(missing(image.DN)){image.DN <- loadImage(path = path)}
    time.line <- grep("SCENE_CENTER_TIME",MTL,value=TRUE)
    date.line <- grep("DATE_ACQUIRED",MTL,value=TRUE)
    sat.time <-regmatches(time.line,regexec(text=time.line,
                                            pattern="([0-9]{2})(:)([0-9]{2})(:)([0-9]{2})(.)([0-9]{2})"))[[1]][1]
    sat.date <-regmatches(date.line,regexec(text=date.line,
                                            pattern="([0-9]{4})(-)([0-9]{2})(-)([0-9]{2})"))[[1]][1]
    sat.datetime <- strptime(paste(sat.date, sat.time), 
                             format = "%Y-%m-%d %H:%M:%S", tz="GMT")
    DOY <-  sat.datetime$yday +1
    d2 <- 1/(1+0.033*cos(DOY * 2 * pi/365))
    dr <- 1 + 0.033 * cos(DOY * (2 * pi / 365))
    Ro.TOAr <- list()
    for(i in 1:6){
      Ro.TOAr[i] <- (pi * (Gain[i] * image.DN[[i]] + Bias[i])) / (ESUN[i] * cos(incidence.rel) * dr)
    }
    image_TOA <- do.call(stack, Ro.TOAr)
  }
  #### 
  image_TOA <- saveLoadClean(imagestack = image_TOA, 
                             stack.names = c("B", "G", "R", "NIR", "SWIR1", "SWIR2"), 
                             file = "image_TOAr", 
                             overwrite=TRUE)
  return(image_TOA)
}  

#' Calculates surface reflectance for L7
#' @description
#' Calculates surface reflectance from top of atmosphere radiance using the model developed by Tasumi et al. (2008) and Allen et al. (2007), which considers a band-by-band basis.
#' @param image.TOAr      raster stack. top of atmosphere reflectance image
#' @param sat             "L7" for Landsat 7, "L8" for Landsat 8 or "auto" to guess from filenames 
#' @param ESPA            Logical. If TRUE will look for espa.usgs.gov related products on working folder
#' @param aoi             area of interest to crop images, if waterOptions("autoAoi") == TRUE will look for any object called aoi on .GlobalEnv
#' @param incidence.hor   solar incidence angle, considering plain surface
#' @param WeatherStation  Weather Station data
#' @param surface.model   rasterStack with DEM, Slope and Aspect. See surface.model()
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David 
#' @references 
#' Tasumi M.; Allen R.G. and Trezza, R. At-surface albedo from Landsat and MODIS satellites for use in energy balance studies of evapotranspiration Journal of Hydrolog. Eng., 2008, 13, (51-63) \cr
#'
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
# incidence hor from TML?? 
calcSR <- function(image.TOAr, sat="auto", ESPA=FALSE, aoi, incidence.hor, 
                   WeatherStation, surface.model){
  if(class(WeatherStation)== "waterWeatherStation"){
    WeatherStation <- getDataWS(WeatherStation)
  }
  path <- getwd()
  if(sat=="auto"){sat = getSat(path)}
  if(sat=="L8"){bands <- 2:7}
  if(sat=="L7"){bands <- c(1:5,7)}
  if(ESPA==TRUE & sat=="L8"){
    files <- list.files(path = path, pattern = "_sr_band+[2-7].tif$", full.names = T)
    stack1 <- list()
    for(i in 1:6){
      stack1[i] <- raster(files[i])
    }
    image_SR <- do.call(stack, stack1)
    image_SR <- aoiCrop(image_SR, aoi) 
    image_SR <- image_SR / 10000
  }
  if(sat=="L7"){
    if(missing(image.TOAr)){image.TOAr <- calcTOAr()}
    P <- 101.3*((293-0.0065 * surface.model$DEM)/293)^5.26
    ea <- (WeatherStation$RH/100)*0.6108*exp((17.27*WeatherStation$temp)/(WeatherStation$temp+237.3))
    W <- 0.14 * ea * P + 2.1
    Kt <- 1
    Cnb <- matrix(data=c(0.987, 2.319, 0.951, 0.375, 0.234, 0.365,
                         -0.00071, -0.00016, -0.00033, -0.00048, -0.00101, -0.00097,
                         0.000036, 0.000105, 0.00028, 0.005018, 0.004336, 0.004296,
                         0.0880, 0.0437, 0.0875, 0.1355, 0.056, 0.0155,
                         0.0789, -1.2697, 0.1014, 0.6621, 0.7757, 0.639,
                         0.64, 0.31, 0.286, 0.189, 0.274, -0.186), byrow = T, nrow=6, ncol=6)
    tau_in <- list()
    tau_out <- list()
    for(i in 1:6){
      tau_in[i] <- Cnb[1,i] * exp((Cnb[2,i]*P/(Kt*cos(incidence.hor)))-
                                    ((Cnb[3,i]*W+Cnb[4,i])/cos(incidence.hor)))+Cnb[5,i]
    }
    eta = 0 # eta it's the satellite nadir angle
    for(i in 1:6){
      tau_out[i] <- Cnb[1,i] * exp((Cnb[2,i]*P/(Kt*cos(eta)))-
                                     ((Cnb[3,i]*W+Cnb[4,i])/cos(eta)))+Cnb[5,i]
    }
    path_refl <- list()
    for(i in 1:6){
      path_refl[i] <- Cnb[6,i] * (1 - tau_in[[i]])
    }
    stack_SR <- list()
    for(i in 1:6){
      stack_SR[i] <- (image.TOAr[[i]] - path_refl[[i]]) / (tau_in[[i]] * tau_out[[i]])
    }
    image_SR <- do.call(stack, stack_SR)
  }
  image_SR <- saveLoadClean(imagestack = image_SR, 
                            stack.names = c("B", "G", "R", "NIR", "SWIR1", "SWIR2"), 
                            file = "image_SR", 
                            overwrite=TRUE)
  return(image_SR)
}  


#' Check needed SRTM grids from image extent
#' @param raw.image  image to calculate extent
#' @author Guillermo Federico Olmedo
#' @export
# Get links or optionally open web pages... 
# Check if the files are present on path o in a specific SRTM local repo
checkSRTMgrids <-function(raw.image){
  path = getwd()
  polyaoi <- SpatialPolygons(
    list(Polygons(list(Polygon(coords = matrix(
      c(xmin(raw.image), xmax(raw.image), xmax(raw.image),
        xmin(raw.image),xmin(raw.image), ymax(raw.image),
        ymax(raw.image), ymin(raw.image), ymin(raw.image),
        ymax(raw.image)), ncol=2, nrow= 5))), ID=1)))
  polyaoi@proj4string <- raw.image@crs
#   limits <- proj4::project(xy = matrix(polyaoi@bbox, ncol=2, nrow=2), proj = polyaoi@proj4string, 
#                            inverse = TRUE)
  limits <- extent(sp::spTransform(polyaoi, CRS("+proj=longlat +ellps=WGS84")))
  # I have to improve this. It should work ONLY for west and south coordinates.. maybe
  lat_needed <- seq(trunc(limits[3])-1, trunc(limits[4])-1, by=1)
  long_needed <- seq(trunc(limits[1])-1, trunc(limits[2])-1, by = 1)
  grids <- expand.grid(lat_needed, long_needed)
  result <- list()
  link <- "http://earthexplorer.usgs.gov/download/options/8360/SRTM1"
  for(i in 1:nrow(grids)){
    result[[i]] <- paste(link, ifelse(grids[i,1]>0,"N", "S"), abs(grids[i,1]),
                         ifelse(grids[i,2]>0,"E", "W"), "0", abs(grids[i,2]),"V3/", sep="")
  }
  print(paste("You need", nrow(grids), "1deg x 1deg SRTM grids"))
  print("You can get them here:")
  return(unlist(result))
}

#' Create a mosaic with SRTM grid from image extent
#' @param format  format of SRTM grid files
#' @param extent  minimal extent of mosaic
#' @author Guillermo Federico Olmedo
#' @export
# Should use checkSRTMgrids to get the files list and not use all from the folder...!
# Also look for files on path and local repo
prepareSRTMdata <- function(format="tif", extent){
  path = getwd()
  files <- list.files(path= path,  
                      pattern=paste("^[sn]\\d{2}_[we]\\d{3}_1arc_v3.", 
                              format, "$", sep=""), full.names = T) 
  if(length(files) < 1){stop(paste("You need to download SRTM grids and save them to working directory \n try: checkSRTMgrids()"))}
  stack1 <- list()
  for(i in 1:length(files)){
    stack1[[i]] <- raster(files[i])}
  stack1$fun <- mean
  if(length(files)>1){SRTMmosaic <- do.call(mosaic, stack1)}
  if(length(files)==1){SRTMmosaic <- stack1[[1]]}
  destino  <-  projectExtent(extent, extent@crs)
  mosaicp <- projectRaster(SRTMmosaic, destino)
  mosaicp <- saveLoadClean(imagestack = mosaicp, stack.names = "DEM", 
                           file = "DEM", overwrite=TRUE)
  return(mosaicp)
}

#' Calculates surface model used in METRIC
#' @description
#' DEM map is used to generate the surface representation of the image through of aspect and slope maps. This procedure helps to avoid differences in the surface temperature (and finally Evapotranspiration) caused by different incidence angles and/or elevations in mountainous areas.
#' @param DEM  raster with Digital elevation model 
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
METRICtopo <- function(DEM){
  aspect <- terrain(DEM, opt="aspect") 
  slope <- terrain(DEM, opt="slope") 
  aspect_metric <- aspect-pi  #METRIC expects aspect - 1 pi
  surface.model <- stack(DEM, slope, aspect_metric)
  surface.model <- saveLoadClean(imagestack = surface.model, 
                                 stack.names = c("DEM", "Slope", "Aspect"), 
                                 file = "surface.model", 
                                 overwrite=TRUE)
  return(surface.model)
}

#' Calculates solar angles
#' @description
#' Metadata, aspect and slope maps are combined to estimate solar angles for the entire image. 
#' @param surface.model   rasterStack with DEM, Slope and Aspect. See surface.model()
#' @param MTL             Landsat Metadata File
#' @param WeatherStation  Weather Station data
#' @details Narrowband transmittances,are calculated considering some radiation transfer models operated 
#' over a wide range of climates and locations across the world, this parameter vary with the cosine of 
#' the solar angle, atmospheric pressure and precipitable water vapor in the atmosphere, so the author must 
#' obtain accurate values of these three parameters.
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David
#' @author Fernando Fuentes Peñailillo
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
### Change to look in metadata for keyword instead of using line #
solarAngles <- function(surface.model, MTL, WeatherStation){
  path = getwd()
  if(class(WeatherStation)== "waterWeatherStation"){
    WeatherStation <- getDataWS(WeatherStation)
  }
  if(missing(MTL)){MTL <- list.files(path = path, pattern = "MTL.txt", full.names = T)}
  MTL <- readLines(MTL, warn=FALSE)
  Elev.line <- grep("SUN_ELEVATION",MTL,value=TRUE)
  sun.elevation <- (90 - as.numeric(regmatches(Elev.line, 
                                               regexec(text=Elev.line ,
                                                       pattern="([0-9]{1,5})([.]+)([0-9]+)"))[[1]][1]))*pi/180
  Azim.line <- grep("SUN_AZIMUTH",MTL,value=TRUE)
  sun.azimuth <- as.numeric(regmatches(Azim.line, 
                                       regexec(text=Azim.line ,
                                               pattern="([0-9]{1,5})([.]+)([0-9]+)"))[[1]][1])
  # latitude
  latitude <- surface.model[[1]]
  xy <- SpatialPoints(xyFromCell(latitude, cellFromRowCol(latitude, 1:nrow(latitude), 1)))
  xy@proj4string <- latitude@crs
  lat <- coordinates( spTransform(xy, CRS("+proj=longlat +datum=WGS84")))[,2] 
  values(latitude) <- rep(lat*pi/180,each=ncol(latitude))
  # declination
  time.line <- grep("SCENE_CENTER_TIME",MTL,value=TRUE)
  date.line <- grep("DATE_ACQUIRED",MTL,value=TRUE)
  sat.time <-regmatches(time.line,regexec(text=time.line,
                                          pattern="([0-9]{2})(:)([0-9]{2})(:)([0-9]{2})(.)([0-9]{2})"))[[1]][1]
  sat.date <-regmatches(date.line,regexec(text=date.line,
                                          pattern="([0-9]{4})(-)([0-9]{2})(-)([0-9]{2})"))[[1]][1]
  sat.datetime <- strptime(paste(sat.date, sat.time), 
                           format = "%Y-%m-%d %H:%M:%S", tz="GMT")
  DOY <-  sat.datetime$yday +1
  declination <- surface.model[[1]]
  values(declination) <- 0.409*sin((2*pi/365*DOY)-1.39)
  # hour angle
  hour.angle <- surface.model[[1]]
  time.line <- grep("SCENE_CENTER_TIME",MTL,value=TRUE)
  time.hours <-regmatches(time.line,regexec(text=time.line,
                                            pattern="([0-9]{2})(:)([0-9]{2})(:)([0-9]{2})(.)([0-9]{2})"))[[1]][1]
  time.hours <- strptime(time.hours, format = "%H:%M:%S", tz="GMT")
  Sc <- (0.1645*sin(4*pi*(DOY-81)/364)-0.1255*cos(2*pi*(DOY-81)/364)-0.025*sin(2*pi*(DOY-81)/364))*60
  time.hours <- (time.hours$hour*3600+time.hours$min*60+time.hours$sec+WeatherStation$long/15*3600)/60 + Sc
  values(hour.angle) <- abs((12-as.numeric(time.hours)/60)*15/180*pi)
  #hour.angle <- asin(-1*(cos(sun.elevation)*sin(sun.azimuth)/cos(declination))) first
  
  ## solar incidence angle, for horizontal surface
  incidence.hor <- acos(sin(declination) * sin(latitude) + cos(declination)*
                          cos(latitude)*cos(hour.angle))
  slope <- surface.model$Slope
  aspect <- surface.model$Aspect
  ##solar incidence angle, for sloping surface
  incidence.rel <- acos(sin(declination)*sin(latitude)*cos(slope) 
                        - sin(declination)*cos(latitude)*sin(slope)*cos(aspect)
                        + cos(declination)*cos(latitude)*cos(slope)*cos(hour.angle)
                        + cos(declination)*sin(latitude)*sin(slope)*cos(aspect)*cos(hour.angle)
                        + cos(declination)*sin(aspect)*sin(slope)*sin(hour.angle))
  ## End
  solarAngles <- stack(latitude, declination, hour.angle, incidence.hor, incidence.rel)
  solarAngles <- saveLoadClean(imagestack = solarAngles, 
                               stack.names = c("latitude", "declination", 
                                               "hour.angle", "incidence.hor", "incidence.rel"), 
                               file = "solarAngles", overwrite=TRUE)
  return(solarAngles)
}

#' Calculates Incoming Solar Radiation
#' @description 
#' This function calculates incoming solar radiation from surface model and solar angles.
#' @param surface.model   rasterStack with DEM, Slope and Aspect. See surface.model()
#' @param solar.angles    rasterStack with latitude, declination, hour.angle, incidence.hor and incidence.rel. See solarAngles()
#' @param WeatherStation  Weather Station data
#' @author Guillermo Federico Olmedo
#' @author Daniel de la Fuente Saiz
#' @author Fonseca-Luengo, David
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
incSWradiation <- function(surface.model, solar.angles, WeatherStation){
  if(class(WeatherStation)== "waterWeatherStation"){
    WeatherStation <- getDataWS(WeatherStation)
  }
  ea <- (WeatherStation$RH/100)*0.6108*exp((17.27*WeatherStation$temp)/(WeatherStation$temp+237.3))
  tau.sw <- SWtrasmisivity(Kt = 1, ea, surface.model$DEM, solar.angles$incidence.hor)
  DOY <- WeatherStation$DOY
  d <- sqrt(1/(1+0.033*cos(DOY * (2 * pi/365))))
  Rs.inc <- 1367 * cos(solar.angles$incidence.rel) * tau.sw / d^2
  Rs.inc <- saveLoadClean(imagestack = Rs.inc, 
                          file = "Rs.inc", overwrite=TRUE)
  return(Rs.inc)
}

#' Calculates Broadband Albedo from Landsat data
#' @description 
#' Broadband surface Albedo is estimated considering the integration of all 
#' narrowband at-surface reflectances following a weighting function with 
#' empirical coefficients (Tasumi et al., 2007).
#' @param image.SR   surface reflectance image with bands B, R, G, NIR, SWIR1, 
#' SWIR2
#' @param aoi        area of interest to crop images, if waterOptions("autoAoi")
#'  == TRUE will look for any object called aoi on .GlobalEnv
#' @param coeff      coefficient to transform narrow to broad band albedo. 
#' See Details.
#' @details 
#' There are differents models to convert narrowband data to broadband albedo. 
#' You can choose coeff="Tasumi" to use Tasumi et al (2008) coefficients, 
#' calculated for Landsat 7; coeff="Liang" to use Liang Landsat 7 coefficients 
#' or "Olmedo" to use Olmedo coefficients for Landsat 8.
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#'
#' M. Tasumi, Allen, R. G., and Trezza, R. 2007. "Estimation of at-surface reflection albedo from satellite for routine operational calculation of land surface energy balance". J. Hydrol. Eng. \cr
#'
#' Liang, S. (2000). Narrowband to broadband conversions of land surface albedo: I. Algorithms. Remote Sensing of Environment, 76(1), 213-238. \cr
#' @export
albedo <- function(image.SR, aoi, coeff="Tasumi"){
  if(coeff=="Tasumi"){wb <- c(0.254, 0.149, 0.147, 0.311, 0.103, 0.036)} 
  # Tasumi 2008
  if(coeff=="Olmedo") {wb <- c(0.246, 0.146, 0.191, 0.304, 0.105, 0.008)}
  # Calculated using SMARTS for Kimberly2-noc13 and Direct Normal Irradiance
  if(coeff=="Liang") {wb <- c(0.356, 0, 0.130, 0.373, 0.085, 0.072)} 
  # Liang 2001
  albedo <- calc(image.SR[[1]], fun=function(x){x *wb[1]})
  for(i in 2:6){
      albedo <- albedo + calc(image.SR[[i]], fun=function(x){x *wb[i]})
      #removeTmpFiles(h=0.0008) # delete last one... maybe 3 seconds
  }
  albedo <- aoiCrop(albedo, aoi) 
  if(coeff=="Liang"){
    albedo <- albedo - 0.0018
  }
  albedo <- saveLoadClean(imagestack = albedo, 
                          file = "albedo", overwrite=TRUE)
  return(albedo)
}

#' Estimate LAI from Landsat Data
#' @description
#' This function implements empirical models to estimate LAI (Leaf Area Index) for satellital images. Models were extracted from METRIC publications and other works developed on different crops.
#' @param method   Method used to estimate LAI from spectral data. 
#' @param image    image. top-of-atmosphere reflectance for method=="metric" | method=="metric2010" | method=="vineyard" | method=="MCB"; surface reflectance for method = "turner". Not needed if ESPA == TRUE
#' @param aoi      area of interest to crop images, if waterOptions("autoAoi") == TRUE will look for any object called aoi on .GlobalEnv
#' @param L        L factor used in method = "metric" or "metric2010" to estimate SAVI, defaults to 0.1
#' @details LAI is computed using the top-of atmosphere (at-satellite) reflectance value. 
#' LAI and other indices such NDVI, SAVI are used to predict characteristics of vegetation, 
#' depending on preferences of the user.
#' Available methods are: "metric", "metric2010", "MCB" and "vineyard". 
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David
#' @author Fernando Fuentes Peñailillo
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' Carrasco-Benavides, M., Ortega-Farias, S., Lagos, L., Kleissl, J., Morales-Salinas, L., & Kilic, A. (2014). Parameterization of the Satellite-Based Model (METRIC) for the Estimation of Instantaneous Surface Energy Balance Components over a Drip-Irrigated Vineyard. Remote Sensing, 6(11), 11342-11371. http://doi.org/10.3390/rs61111342\cr
#' Johnson, L. F. (2003). Temporal Stability of the NDVI-LAI Relationship in a Napa Valley Vineyard, 96-101. http://doi.org/10.1111/j.1755-0238.2003.tb00258.x
#' @export
## Cite Pocas work for LAI from METRIC2010
LAI <- function(method="metric2010", image, aoi, L=0.1){
  if(method=="metric" | method=="metric2010" | method=="vineyard" | method=="MCB"){
      toa.4.5 <- stack(image[[3]], image[[4]])}
  if(method=="turner"){
      sr.4.5 <- stack(image[[3]], image[[4]])}
  if(method=="turner"){
    sr.4.5 <- sr.4.5 * 0.0001}
  if(method=="metric2010"){
    SAVI_ID <- (1 + L)*(toa.4.5[[2]] - toa.4.5[[1]])/(L + toa.4.5[[1]] + 
                                                        toa.4.5[[2]])
    SAVI_ID <- saveLoadClean(imagestack = SAVI_ID, stack.names = "SAVI_ID", 
                             file = "SAVI_ID", overwrite=TRUE)
    LAI <- raster(SAVI_ID)
    LAI[SAVI_ID <= 0] <- 0
    LAI[SAVI_ID > 0 & SAVI_ID <= 0.817] <- 11 * SAVI_ID[SAVI_ID > 0 & 
                                                          SAVI_ID <= 0.817]^3 # for SAVI <= 0.817
    LAI[SAVI_ID > 0.817] <- 6
  }
  if(method=="metric"){
    SAVI_ID <- (1 + L)*(toa.4.5[[2]] - toa.4.5[[1]])/(L + toa.4.5[[1]] + toa.4.5[[2]])
    LAI <- log((0.69-SAVI_ID)/0.59)/0.91 *-1
    LAI[SAVI_ID > 0.817] <- 6
  }
  if(method=="vineyard"){
    NDVI <- (toa.4.5[[2]] - toa.4.5[[1]])/(toa.4.5[[1]] + toa.4.5[[2]])
    LAI <- 4.9 * NDVI -0.46 # Johnson 2003
  }
  ## method carrasco
  if(method=="MCB"){
    NDVI <- (toa.4.5[[2]] - toa.4.5[[1]])/(toa.4.5[[1]] + toa.4.5[[2]])
    LAI <- 1.2 - 3.08*exp(-2013.35*NDVI^6.41) 
  }
  if(method=="turner"){
    NDVI <- (sr.4.5[[2]] - sr.4.5[[1]])/(sr.4.5[[1]] + sr.4.5[[2]])
    LAI <- 0.5724+0.0989*NDVI-0.0114*NDVI^2+0.0004*NDVI^3 # Turner 1999
  }
  LAI[LAI<0]  <- 0
  LAI <- saveLoadClean(imagestack = LAI, stack.names = "LAI", 
                       file = "LAI", overwrite=TRUE)
  return(LAI)
}

#' Calculates short wave transmisivity
#' @description
#' Short wave transmisivity is estimated for broad-band considering an extended equation developed by Allen (1996), based from Majumdar et al.(1972), using coefficients developed by ASCE-EWRI (2005).
#' @param Kt            unitless turbidity coefficient 0<Kt<=1.0, where Kt=1.0 for clean air and Kt=0.5 for extremely turbid, dusty, or polluted air
#' @param ea            near-surface vapor pressure (kPa)
#' @param dem           digital elevation model 
#' @param incidence.hor solar incidence angle, considering plain surface
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#'
#' Majumdar, N.; Mathur, B. & Kaushik, S. Prediction of direct solar radiation for low atmospheric turbidity Solar Energy, Elsevier, 1972, 13, 383-394 \cr
#'
#' ASCE-EWRI The ASCE Standardized Reference Evapotranspiration Equation Report of the ASCE-EWRI Task Committee on Standardization of Reference Evapotranspiration, 2005 \cr
#' @export
SWtrasmisivity <- function(Kt = 1, ea, dem, incidence.hor){
  P <- 101.3*((293-0.0065 * dem)/293)^5.26
  W <- 0.14 * ea * P + 2.1
  tauB <- 0.98 * exp(((-0.00149 * P )/ (Kt * 
                      cos(incidence.hor)))-0.075*((W / cos(incidence.hor))^0.4))
  tauD <- raster(tauB)
  tauD[tauB >= 0.15] <- 0.35 - 0.36 * tauB 
  tauD[tauB < 0.15] <- 0.18 + 0.82 * tauB  
  tau.sw <- tauB + tauD
  # Next one it's from METRIC 2007, previous from METRIC 2010
  #sw.t <- 0.35 + 0.627 * exp((-0.00149 * P / Kt * 
  #                              cos(incidence.hor))-0.075*(W / cos(incidence.hor))^0.4)
  return(tau.sw)
}

#' Estimates Land Surface Temperature from Landsat Data
#' @description
#' Surface temperature is estimated using a modified Plank equation considering empirical constants for Landsat images. In addition, this model implements a correction of thermal radiance following to Wukelic et al. (1989).
#' @param thermalband     Satellite thermal band
#' @param sat             "L7" for Landsat 7, "L8" for Landsat 8 or "auto" to guess from filenames 
#' @param LAI             raster layer with leaf area index. See LAI()
#' @param aoi             area of interest to crop images, if waterOptions("autoAoi") == TRUE will look for any object called aoi on .GlobalEnv
#' @param WeatherStation  Weather Station data
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#'
#' Wukelic G. E.; Gibbons D. E.; Martucci L. M. & Foote, H. P. Radiometric calibration of Landsat thematic mapper thermal band Remote Sensing of Environment, 1989, 28, (339-347) \cr
#' @export
## Add Sobrino and Qin improvements to LST in ETM+
## Add Rsky estimation from WeatherStation
surfaceTemperature <- function(thermalband, sat="auto", LAI, aoi, 
                               WeatherStation){
  path=getwd()
  if(class(WeatherStation)== "waterWeatherStation"){
    WeatherStation <- getDataWS(WeatherStation)
  }
  if(sat=="auto"){sat = getSat(path)}
  if(sat=="L8"){
    if(missing(thermalband)){
      bright.temp.b10 <- raster(list.files(path = path, 
                                         pattern = "_toa_band10.tif"))
      bright.temp.b10 <- aoiCrop(bright.temp.b10, aoi) 
      Ts <- bright.temp.b10 * 0.1
    } else {
      bright.temp.b10 <- aoiCrop(thermalband, aoi) 
      Ts <- bright.temp.b10 * 0.1
    }
  }
  if(sat=="L7"){
    epsilon_NB <- raster(LAI)
    epsilon_NB <- 0.97 + 0.0033 * LAI  
    epsilon_NB[LAI > 3] <- 0.98
    if(missing(thermalband)){
      B6 <- raster(list.files(path = path, 
                pattern = "^L[EC]\\d+\\w+\\d+_B6_VCID_1.TIF$", full.names = T))
    }
    if(!missing(thermalband)){B6 <- thermalband}
    L_t_6 <-  0.067 * B6 - 0.06709
    if(missing(aoi) & exists(x = "aoi", envir=.GlobalEnv)){
      aoi <- get(x = "aoi", envir=.GlobalEnv)
      L_t_6 <- crop(L_t_6,aoi)
    }
    L7_K1 <- 666.09 
    L7_K2 <- 1282.71 
    Rp <- 0.91             #Allen estimo en Idaho que el valor medio era 0.91
    #tau_NB <- 1       #Allen estimo en Idaho que el valor medio era 0.866
    ## There is a equation on excel to calculate from DEM
    tau_NB <- 0.75+2e-5*WeatherStation$elev
    #R_sky <- 1        #Allen estimo en Idaho que el valor medio era 1.32
    ## There is a equation for R_sky on Metric Manual:
    R_sky <- (1.807e-10)*(WeatherStation$temp+273.15)^4*(1-0.26*
                                            exp(-7.77e-4*WeatherStation$temp^2))
    Rc <- ((L_t_6 - Rp) / tau_NB) - (1-epsilon_NB)*R_sky
    Ts <- L7_K2 / log((epsilon_NB*L7_K1/Rc)+1)}
  Ts <- saveLoadClean(imagestack = Ts, 
                      file = "Ts", overwrite=TRUE)
  return(Ts)
}

#' Calculates Long wave outgoing radiation
#' @description
#' This function estimates the long wave outgoing radiation using the Stefan-Boltzmann equation.
#' @param LAI  raster layer with leaf area index. See LAI()
#' @param Ts   Land surface temperature. See surfaceTemperature()
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
outLWradiation <- function(LAI, Ts){
  surf.emissivity <- 0.95 + 0.01 * LAI
  surf.emissivity[LAI>3] <- 0.98
  Rl.out <- surf.emissivity * 5.67e-8 * Ts^4
  Rl.out <- saveLoadClean(imagestack = Rl.out, 
                          file = "Rl.out", overwrite=TRUE)
  return(Rl.out)
}

#' Calculates long wave incoming radiation
#' @description
#' This function estimates the long wave incoming radiation using the Stefan-Boltzmann equation. In addition, empirical equation of Bastiaanssen (1995) with coefficients developed by Allen (2000) are used to estimate the effective atmospheric emissivity.
#' @param WeatherStation   Weather Station data
#' @param DEM              digital elevation model in meters.
#' @param solar.angles     rasterStack with latitude, declination, hour.angle, incidence.hor and incidence.rel. See solarAngles()
#' @param Ts               Land surface temperature. See surfaceTemperature()
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#'
#' Bastiaanssen W. Regionalization of surface flux densities and moisture indicators in composite terrain: A remote sensing approach under clear skies in Mediterranean climates. Ph.D. dissertation, CIP Data Koninklijke Bibliotheek, Den Haag, The Netherlands, 1995, p273 \cr
#'
#' Allen R. RAPID long-wave radiation calculations and model comparisons Internal report, University of Idaho, Kimberly, Idaho, 2000 \cr
#' @export
# Add a function to get "cold" pixel temperature so in can be used in the next function
incLWradiation <- function(WeatherStation, DEM, solar.angles, Ts){
  if(class(WeatherStation)== "waterWeatherStation"){
    WeatherStation <- getDataWS(WeatherStation)
  }
  ea <- (WeatherStation$RH/100)*0.6108*exp((17.27*WeatherStation$temp)/
                                             (WeatherStation$temp+237.3))
  tau.sw <- SWtrasmisivity(Kt = 1, ea, DEM, solar.angles$incidence.hor)
  #temp <-  WeatherStation$temp+273.15 - (DEM - WeatherStation$elev) * 6.49 / 1000 
  ## Temperature in Kelvin
  #Mountain lapse effects from International Civil Aviation Organization
  ef.atm.emissivity  <- 0.85*(-1*log(tau.sw))^0.09
  Rl.in <- ef.atm.emissivity * 5.67e-8 * Ts^4
  #Rl.in <- ef.atm.emissivity * 5.67e-8 * temp^4
  Rl.in <- saveLoadClean(imagestack = Rl.in, 
                         file = "Rl.in", overwrite=TRUE)
  return(Rl.in)
}

#' Estimates net radiation
#' @description
#' This function estimates net radiation considering a surface radiation balance. Equations use the information from the image, in addition of measurements of actual vapor pressure and altitude.
#' @param LAI     raster layer with leaf area index. See LAI()
#' @param albedo  broadband surface albedo. See albedo()
#' @param Rs.inc  incoming short-wave radiation
#' @param Rl.inc  incomin long-wave radiation
#' @param Rl.out  outgoing long-wave radiation
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
netRadiation <- function(LAI, albedo, Rs.inc, Rl.inc, Rl.out){
  surf.emissivity <- 0.95 + 0.01 * LAI 
  surf.emissivity[LAI>3] <- 0.98
  Rn <- (1- albedo)*Rs.inc + Rl.inc - Rl.out - (1-surf.emissivity)*Rl.inc
  Rn <- saveLoadClean(imagestack = Rn, 
                      file = "Rn", overwrite=TRUE)
  return(Rn)
}


#' Estimates Soil Heat Flux
#' @description
#' This function implements models to estimate soil heat flux for different surfaces and considering different inputs.
#' @param image    surface reflectance image
#' @param Ts       Land surface temperature. See surfaceTemperature()
#' @param albedo   broadband surface albedo. See albedo()
#' @param LAI      raster layer with leaf area index. See LAI()
#' @param Rn       Net radiation. See netRadiation()
#' @param aoi      area of interest to crop images, if waterOptions("autoAoi") == TRUE will look for any object called aoi on .GlobalEnv
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
soilHeatFlux <- function(image, Ts, albedo, LAI, Rn, aoi){
  sr.4.5 <- stack(image[[3]], image[[4]])
  NDVI <- (sr.4.5[[2]] - sr.4.5[[1]])/(sr.4.5[[1]] + sr.4.5[[2]])
  G <- ((Ts - 273.15)*(0.0038+0.0074*albedo)*(1-0.98*NDVI^4))*Rn
  G <- saveLoadClean(imagestack = G, file = "G", overwrite=TRUE)
  G <- raster(NDVI)
  e <- 2.71828
  G <- (0.05 + 0.18 * e^(-0.521*LAI)) * Rn
  G[LAI < 0.5] <- ((1.8*(Ts[LAI < 0.5] - 273.16) / Rn[LAI < 0.5]) + 0.084) * 
    Rn[LAI < 0.5]
  return(G)
}




