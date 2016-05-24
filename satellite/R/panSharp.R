if ( !isGeneric("panSharp") ) {
  setGeneric("panSharp", function(x, ...)
    standardGeneric("panSharp"))
}

#' Pan sharpen low resolution satellite channels by using the high resolution 
#' panchromatic channel.
#'
#' @description The function PAN sharpens the low resolution channels with the 
#' panchromatic channel. This is done by multiplying the normlized XS channel 
#' with the PAN channel (see Details).
#'
#' @param x Satellite or \code{raster::Raster*} object.
#' @param pan A raster::RasterLayer object of the panchromatic channel
#' @param pan_lp A raster::RasterLayer object containing a lowpass filtering of
#' pan 
#' @param filter Type of filter to be used for smoothing the PAN raster; one of 
#'  mean (default), Gauss, median.
#' @param winsize Size of the filter window in x and y direction; defaults to 3.
#' @param subset Logical; if \code{TRUE}, all layers except for the cropped ones 
#' are being dropped; if \code{FALSE}, the cropped layers are being appended to
#' the Satellite object.
#'
#' @return If x is a Satellite object, a Satellite object (with added 
#' pansharpened layers); if x is a \code{raster::Raster*} object, a 
#' \code{raster::Raster*} with pansharpened layer(s).
#' 
#' @export panSharp
#' 
#' @name panSharp
#'
#' @details Pan sharpen low resolution satellite channels by using the high 
#' resolution panchromatic channel. This function uses the same algorithm as the 
#' OTB Toolbox where "The idea is to apply a low pass filter to the 
#' panchromatic band to give it a spectral content (in the Fourier domain) 
#' equivalent to the XS data. Then we normalize the XS data with this low-pass 
#' panchromatic and multiplythe result with the original panchromatic band." 
#' (see \url{https://www.orfeo-toolbox.org/SoftwareGuide/SoftwareGuidech13.html#x41-2140011}).
#' 
#' @references 
#' Al-amri, Salem Saleh, Namdeo V. Kalyankar, and Santosh D. Khamitkar. "A 
#' comparative study of removal noise from remote sensing image." 
#' \url{http://ijcsi.org/articles/A-Comparative-Study-of-Removal-Noise-from-Remote-Sensing-Image.php}
#' 
#' Bhattacharya, Amit K., P. K. Srivastava, and Anil Bhagat. "A modified texture 
#' filtering technique for satellite images."
#' Paper presented at the 22nd Asian Conference on Remote Sensing. Vol. 5. 2001.
#' \url{http://a-a-r-s.org/aars/proceeding/ACRS2001/Papers/DPA3-08.pdf}
#' 
#' Randen, Trygve, and John Hakon Husoy. "Filtering for texture classification: 
#' A comparative study." Pattern Analysis and Machine Intelligence, IEEE 
#' Transactions on 21.4 (1999): 291-310. \url{http://dx.doi.org/10.1109/34.761261}.
#' 
#' PAN sharpening articles \cr
#' - \url{http://remotesensing.spiedigitallibrary.org/article.aspx?articleid=1726558} \cr
#' - \url{http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=1368950&url=http\%3A\%2F\%2Fieeexplore.ieee.org\%2Fxpls\%2Fabs_all.jsp\%3Farnumber\%3D1368950}
#'
#' @examples 
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' sat <- satellite(files)
#' 
#' \dontrun{
#' ## using 'satellite' object
#' sat_ps <- panSharp(sat)
#' 
#' par(mfrow = c(1, 2))
#' plot(getSatDataLayer(sat_ps, "B002n"), main = "raw", legend = TRUE)
#' plot(getSatDataLayer(sat_ps, "B002n_PAN_sharpend"), 
#'      main = "pan-sharpened", legend = TRUE)
#' dev.off()
#' }      
#' 
#' ## using 'RasterLayer' object
#' rst_b001n <- getSatDataLayer(sat, "B001n")
#' rst_panch <- getSatDataLayer(sat, getSatBCDEFromType(sat, type = "PCM"))
#' rst_b001n_ps <- panSharp(rst_b001n, rst_panch)
#' 
#' par(mfrow = c(1, 2))
#' plot(rst_b001n, main = "raw", legend = FALSE)
#' plot(rst_b001n_ps, main = "pan-sharpened", legend = FALSE)
#' dev.off()      
NULL


# Function using satellite object ----------------------------------------------
#' 
#' @rdname panSharp
#'
setMethod("panSharp", 
          signature(x = "Satellite"), 
          function(x, filter = c("mean", "Gauss", "median"), 
                   winsize = 1, subset = FALSE){
            
            pan <- getSatDataLayer(x, getSatBCDEFromType(x, type = "PCM"))
            pan_lp <- satellite:::pan_lpf(pan = pan, filter = filter[1], 
                                          winsize = winsize)
            
            bcde_solar <- getSatBCDEFromSpectrum(x, spectrum = "solar")
            bcde_solar <- 
              bcde_solar[!bcde_solar %in% getSatBCDEFromType(x, type = "PCM")]
            
            for(act_bcde in bcde_solar){
              act_pan <- panSharp(x = getSatDataLayer(x, act_bcde),
                                  pan = pan,
                                  pan_lp = pan_lp)
              
              layer_bcde <- paste0(act_bcde, "_PAN_sharpend")
              
              meta_param <- data.frame(getSatSensorInfo(x),
                                       getSatBandInfo(x, act_bcde, 
                                                      return_calib = FALSE),
                                       CALIB = "PAN_sharpend",
                                       createRasterMetaData(pan))
              info <- sys.calls()[[1]]
              info <- paste0("Add layer from ", info[1], "(", 
                             toString(info[2:length(info)]), ")")
              x <- addSatDataLayer(x, bcde = layer_bcde, data = act_pan,
                                   meta_param = meta_param,
                                   info = info, in_bcde = act_bcde)
            }
            if(subset == TRUE){
              x <- subset(x, cid = "PAN_sharpend")
            }
            return(x)
          }
)


# Function using raster::RasterStack object ------------------------------------
#' 
#' @rdname panSharp
#'
setMethod("panSharp", 
          signature(x = "RasterStack"),
          function(x, pan, filter = c("mean", "Gauss", "median"), winsize = 1){
            pan_lf <- pan_lpf(pan = pan, filter = filter[1], winsize = winsize)
            for(l in seq(nlayers(x))){
              x[[l]] <- panSharp(x = x[[l]], pan, pan_lf)  
              return(x)
            }
          }
)


# Function using raster::RasterLayer object ------------------------------------
#' 
#' @rdname panSharp
#'
setMethod("panSharp", 
          signature(x = "RasterLayer"), 
          function(x, pan, pan_lp, filter = c("mean", "Gauss", "median"), 
                   winsize = 1){
            if(missing(pan_lp)){
              pan_lp <- pan_lpf(pan, filter = filter[1], winsize = winsize)
            }
            x_panres <- raster::resample(x, pan, method = "ngb")
            x_pan <- x_panres / pan_lp * pan
            return(x_pan)
          }         
)


# Helper function for the computation of the low pass filter -------------------
#
#pan PAN raster
#filter Filter to be used
#winsizse Window size of filter
# 
pan_lpf <- function(pan, filter, winsize){
  filter_size <- winsize * res(pan)[1]
  sigma <- 3
  switch(filter,
         mean = {
           ftype <- raster::focalWeight(pan, d = filter_size, 
                                        type = "rectangle")
           fun <- sum
         },
         Gauss = {
           ftype <- raster::focalWeight(pan, d=c(sigma, filter_size), 
                                        type = "Gauss")
           fun <- sum
         },
         median = {
           ftype <- raster::focalWeight(pan, d = filter_size, 
                                        type = "rectangle" )
           fun <- stats::median
         }
  )
  ref <- raster::focal(pan, w = ftype, fun = fun)
  return(ref)
}
