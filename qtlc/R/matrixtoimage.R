#' Converts matrix back to image
#'
#' Using S3 qtlc object, extracts the matrix and converts to image plot.
#' @param object S3 object of working TLC.
#' @param show Boolean, default TRUE. Shows the plot of the image.
#' @param bkg If \code{show} is TRUE, then defines background color. Default is "thistle".
#' @param axes Boolean, default FALSE. Shows x,y axes if TRUE.
#' @param xlab Label of the x-axis.
#' @param ylab Label of the y-axis.
#' @param ... Additional graphical parameters.
#'
#' @return Returns image as matrix suitable for \code{plot}, or other graphics functinos.
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' # Converts test image to matrix,
#' # then matrix back to image.
#' fname01 <- system.file("extdata", "testTIFF.tiff", package="qtlc")
#' testTLC <- createTLC(fname01, RGB=FALSE)
#' print(testTLC)
#' matrixtoimage(testTLC, bkg="white")
#'
#' @export
#' @importFrom graphics par plot rasterImage
#' 
matrixtoimage <- function(object, show=TRUE, bkg = "thistle", axes=FALSE, xlab="", ylab="", ...) {
    mat.sl <- 1 - object$mat;
    mat.sl <- rotatev(mat.sl);
    mat.sl <- t(mat.sl);

    if(show==TRUE && exists("rasterImage")) {		 
		graphics::par(bg=bkg);
		graphics::plot(0,0, xlim=c(0, ncol(mat.sl)), ylim=c(0, nrow(mat.sl)), type="n", xlab=xlab, ylab=ylab, axes=axes, asp=1, ...);
		graphics::rasterImage(mat.sl, 0, 0, ncol(mat.sl), nrow(mat.sl));
		return(TRUE);
	} else {
		return(mat.sl);
	}
}

