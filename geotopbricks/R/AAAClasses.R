NULL
#'
#' A GeotopRasterBrick: an object to manage raster maps provied by GEOtop!!

#' @name zoo-class
#' @rdname zoo-class
#' @exportClass zoo 
#' 
#'   @keywords classes
#' @examples showClass("zoo")


##   @exportClass varest
## 
## #    
## #  \code{varest} S3 class (formal definition) see \code{\link{VAR}} 
## #  
## 
## #  	\describe{
## #      The details of the class are reported on \code{\link{VAR}} documentation in "vars" package
## #   }
## 
## #  
## #  
## #  @title varest-class 
## #  
## #  @note Formal definition with \code{\link{setOldClass}} for the S3 class \code{varest2} 
## #  subset
## #  

##   
##   @author Bernhard Pfaff 
##   
##   @docType class
##   @aliases varest
##   @name varest-class
##   @rdname varest-class
##   
##   @keywords classes
##   @exportClass varest
##  
##   @examples showClass("varest")
##   
##    
## 
setOldClass("zoo")

NULL
#' 
#' A GeotopRasterBrick: an object to manage raster maps provied by GEOtop!! 
#' 
#'
#'  \describe{
#'     \item{\code{ascpath}:}{A \code{"zoo"} S3 object containing the names of ascii maps provided by GEOtop 
#' 
#'  }
#'  \item{\code{index}:}{A \code{"POSIXt"} S3 object containing time or dates on which raster layers of \code{brick} are referred
#' 
#'  }
#'  \item{\code{layer}:}{character. Name of the vertical layer at which raster map are referred
#' 
#'  }
#' \item{\code{brick}:}{A \code{"RasterBrick-class"} S4 object containing the  Raster-Layer maps imported from GEOtop output files   
#' 
#'  }   
#'}   
#' 
#' #' @note A \code{GeotopRasterBrick} object can be created by \code{new("GeotopRasterBrick", ...)} 
#' 
#' 
#' 
#' @docType class 
#' @title GeotopRasterBrick-class 
#' 
#' @keywords classes
#' @seealso \code{\link{Raster-class}}
#' 
#' @author Emanuele Cordano
#' @aliases GeotopRasterBrick
#' @name GeotopRasterBrick-class
#' @rdname GeotopRasterBrick-class
#' @exportClass GeotopRasterBrick
#' 
#' 
#' @examples showClass("GeotopRasterBrick")
#' 
#' 
#'  
#' @exportClass GeotopRasterBrick 
#' 

require("raster")

setClass("GeotopRasterBrick",representation(ascpath="zoo",index="POSIXt",layer="character",brick="RasterBrick"))