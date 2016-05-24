if ( !isGeneric("calcHistMatch") ) {
  setGeneric("calcHistMatch", function(x, ...)
    standardGeneric("calcHistMatch"))
}
#' Illumination correction across scenes using histogram matching
#'
#' @description
#' This function adjusts the illumination of individual bands across two scenes
#' using a histogram match.
#'
#' @param x Satellite or raster::Raster* object providing the source band(s) to 
#' be adjusted.
#' @param bcde Band code which should be alligned
#' @param target The target band as raster::RasterLayer.
#' @param minv Lower limit of the possible range for transformation (if not 
#' provided, defaults to the minimum of both layers).
#' @param maxv Upper limit of the possible range for transformation (if not 
#' provided, defaults to the maximum of both layers).
#' @param use_cpp Logical. If \code{TRUE}, C++ functionality (via \strong{Rcpp}) 
#' is enabled, which leads to a considerable reduction of both computation time
#' and memory usage.
#'
#' @name calcHistMatch
#' 
#' @export calcHistMatch
#' 
#' @details The function is based on a histogram matching technique described
#' by Morovic et al. (2002).
#' 
#' @references Morovic J, Shaw J, Sun P-L (2002) A fast, non-iterative and exact 
#' histogram matching algorithm. Pattern Recognition Letters 23/1-3: 127-135, 
#' doi:10.1016/S0167-8655(01)00107-6.
#'
#' @examples
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' sat <- satellite(files)
#' target <- getSatDataLayer(sat, "B003n")
#'  
#' \dontrun{
#' ## histogram matching
#' calcHistMatch(sat, target, bcde = "B002n")
#' }
NULL


# Function using satellite object ----------------------------------------------
#' 
#' @return Satellite object with added atmospheric corrected layers
#' 
#' @rdname calcHistMatch
#'
setMethod("calcHistMatch", 
          signature(x = "Satellite"), 
          function(x, target, bcde = NULL, minv = 0L, maxv = 1023L,
                   use_cpp = TRUE){
            
            if(is.null(bcde)){
              bcde <- c(as.character(getSatBCDESolar(x)), 
                        as.character(getSatBCDEThermal(x)))
            }
            
            for(act_bcde in bcde){
              if(class(target) == "Satellite"){
                act_target <- getSatDataLayer(target, act_bcde)
              } else {
                act_target <- target
              }
              hm <- calcHistMatch(x = getSatDataLayer(x, act_bcde),
                                  target = act_target,
                                  minv = minv,
                                  maxv = maxv,
                                  use_cpp = use_cpp)
              
              meta_param <- getSatMeta(x, act_bcde)
              meta_bcde <- paste0(act_bcde, "_histm")
              
              info <- sys.calls()[[1]]
              info <- paste0("Adjust layer ", act_bcde, 
                             " using histogram matching.")
              x <- addSatDataLayer(x, bcde = meta_bcde, data = hm, 
                                   meta_param = meta_param,
                                   info = info, in_bcde = act_bcde)
            }
            return(x)
          })


# Function using raster::RasterStack object ------------------------------------
#' 
#' @return raster::RasterStack object with atmospheric corrected layers
#' 
#' @rdname calcHistMatch
#'
setMethod("calcHistMatch", 
          signature(x = "RasterStack"), 
          function(x, target, minv = 0L, maxv = 1023L, use_cpp = TRUE){
            # If not supplied, 'model' defaults to DOS2
            model <- model[1]
            
            for(l in seq(nlayers(x))){
              x[[l]] <- calcHistMatch(x = x, target = target, 
                                      minv = minv, maxv = maxv, 
                                      use_cpp = use_cpp)
            }
            return(x)
          })


# Function using raster::RasterLayer object ------------------------------------
#' 
#' @return raster::RasterLayer object with atmospheric corrected layer
#' 
#' @rdname calcHistMatch
#'
setMethod("calcHistMatch", 
          signature(x = "RasterLayer"), 
          function(x, target, minv = 0L, maxv = 1023L, use_cpp = TRUE){
            
            x <- round((x - minValue(x)) * (maxv - minv) / 
                         (maxValue(x) - minValue(x)) + minv)
            
            target <- round((target - minValue(target)) * (maxv - minv) / 
                              (maxValue(target) - minValue(target)) + minv)
            
            breaks <- seq(minv, maxv + 1)
            hs <- hist(x, breaks = breaks, right = FALSE)
            ht <- hist(target, breaks = breaks, right = FALSE)
            
            ## enable c++ functionality
            if (use_cpp) {
              t <- insertMinReqRem(hs$counts, ht$counts)
              
              ## or stick to base-r version  
            } else {
              t <- matrix(data = 0, nrow = length(hs$counts), 
                          ncol = length(ht$counts))
              
              for(j in seq(length(ht$counts))){
                for(i in seq(length(hs$counts))){
                  pixelsreq <- ht$counts[j] - sum(t[1:i,j], na.rm = TRUE)
                  pixelsrem <- hs$counts[i] - sum(t[i,1:j], na.rm = TRUE)
                  t[i,j] <- min(pixelsreq, pixelsrem)
                }
              }
            }
            
            df <-getValues(x)
            for(i in seq(length(df))){
              i_t <- df[i] + 1
              cpf <- cumsum(t[i_t,])
              set.seed(i)
              p <- sample(seq(1, cpf[length(cpf)]), 1)
              j <- which(p <= cpf)[1]
              df[i] <- j #ht$mids[j] + 0.5
              t[i_t,j] <- t[i_t,j] - 1
            }
            x <- setValues(x, df - 1)
          })  
