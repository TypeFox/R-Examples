#' Calculate textures from grey-level co-occurrence matrices (GLCMs) in R
#'
#' Enables calculation of image textures derived from grey-level co-occurrence 
#' matrics (GLCMs) in R. The texture calculation is coded in C++ to optimize 
#' computation time.
#'
#' @name glcm-package
#' @docType package
#' @author Alex Zvoleff, \email{azvoleff@@conservation.org}
#' @keywords package
#' @useDynLib glcm
NULL
#' Landsat 5 Surface Reflectance Image from February 6, 1986 (path 15, row 53)
#' 
#' Portion of Landsat 5 Surface Reflectance image from the Landsat Climate Data 
#' Record archive. This subset of the image includes only bands 1-4.
#'
#' @docType data
#' @name L5TSR_1986
NULL
#' Randomly generated 100x100 test image
#'
#' Used in testing the output from the GLCM texture statistics C++ code.
#' 
#' @docType data
#' @name test_raster
#' @examples
#' # The image was generated with the following code:
#' require(raster)
#' set.seed(0)
#' test_matrix <- matrix(runif(100)*32, nrow=10)
#' test_raster <- raster(test_matrix, crs='+init=EPSG:4326')
#' test_raster <- cut(test_raster, seq(0, 32))
NULL
#' GLCM textures calculated in EXELIS ENVI (for testing purposes)
#'
#' This is the output from running a "co-occurrence measures" calculation to 
#' calculate GLCM textures in EXELIS ENVI from the \code{test_raster} included 
#' in the \code{glcm} package. The following settings were used: window size 
#' 3x3; co-occurrence shift 1 row (y in ENVI), 1 column (x in ENVI); greyscale 
#' textures to compute: mean, variance, homogeneity, contrast, dissimilarity, 
#' entropy, second moment, correlation.
#' 
#' @docType data
#' @name expected_textures_3x3_1x1
#' @seealso  \code{\link{expected_textures_5x7_2x3}} 
#' \code{\link{expected_textures_5x3_n1xn2}}
NULL
#' GLCM textures calculated in EXELIS ENVI (for testing purposes)
#'
#' This is the output from running a "co-occurrence measures" calculation to 
#' calculate GLCM textures in EXELIS ENVI from the \code{test_raster} included 
#' in the \code{glcm} package. The following settings were used: window size 
#' 5x7; co-occurrence shift 2 rows (y in ENVI), 3 columns (x in ENVI); 
#' greyscale textures to compute: mean, variance, homogeneity, contrast, 
#' dissimilarity, entropy, second moment, correlation.
#' 
#' @docType data
#' @name expected_textures_5x7_2x3
#' @seealso \code{\link{expected_textures_3x3_1x1}} 
#' \code{\link{expected_textures_5x3_n1xn2}}
NULL
#' GLCM textures calculated in EXELIS ENVI (for testing purposes)
#'
#' This is the output from running a "co-occurrence measures" calculation to 
#' calculate GLCM textures in EXELIS ENVI from the \code{test_raster} included 
#' in the \code{glcm} package. The following settings were used: window size 
#' 5x3; co-occurrence shift -1 row (y in ENVI), -2 columns (x in ENVI); 
#' greyscale quantization levels 32; textures to compute: mean, variance, 
#' homogeneity, contrast, dissimilarity, entropy, second moment, correlation.
#' 
#' @docType data
#' @name expected_textures_5x3_n1xn2
#' @seealso \code{\link{expected_textures_3x3_1x1}} 
#' \code{\link{expected_textures_5x7_2x3}}
NULL
