#################################################
# qtlc - library for densitometric quantification
#
# Developer and Maintainer: Ivan D. Pavicevic
#                           ivanp84@gmail.com
#
# S3 OO model
# Created: 2015-01-11
#
## Last modified: 2015-12-22
## Library: 2015-12-22
#################################################



#' Creates TLC S3 object
#'
#' Create matrix from TLC image
#' @param ttiff File name of the TIFF image with scanned TLC plate.
#' @param turnv Boolean value determines to turn vertically data in the matrix. TRUE generates turned image which is useful for Cartesian coordinates, because without turning the coordinate system begins in the left corner of the monithor and rises left and down.
#' @param ... Additional parameters for TIFF image manipulation.
#'
#' @return An object of class \code{qtlc}, that contains TLC matrix and descriptions. The object contains: \item{file_name}{File name of of the TIFF image from which the TLC matrix was created.} \item{mat}{TLC matrix (or matrices if intensities Red, Green and Blue channels are not combined.)} \item{spots}{Coordinates of marked spots (using function \code{spot2D}).}
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' @keywords device array
#' 
#' @examples
#' #Creates test matrix.
#' # RGB channels stay separated, or
#' # intensities are combined.
#' fname01 <- system.file("extdata", "testTIFF.tiff", package="qtlc")
#' testTLC <- createTLC(fname01, RGB=TRUE)
#' print(testTLC)
#'
#' @export
#' @importFrom tiff readTIFF
#' 
createTLC <- function(ttiff, turnv=TRUE, ...) {
	
	object <- list();
	class(object) <- "qtlc";
	object$file_name <- ttiff;
	object$mat <- picmatrixTIFF(ttiff, ...);
	if (turnv == TRUE) {
		if (!is.na(dim(object$mat)[3])) {
			for(i in 1:3) {object$mat[,,i] <- rotatev(object$mat[,,i]);}}
			else {object$mat <- rotatev(object$mat);}}
	
	object$spots = NULL;		
	return(object);
    }
