#' Shows 3D plot of the TLC matrix.
#'
#' The function uses TLC matrix width, height and intensity values to make 3D plot.
#' @param object S3 object of the working TLC.
#' @param spot If the specific spot should be represented in 3D, but not entire TLC matrix. (Spot number is given as value, and spots are counted from left to right.)
#' @param ogl If TLC 3D plot use OpenGL library for fast and interactive 3D plot. (This functionality is based on the \code{rgl} package.) Otherwise the function is based on the plot3D package.
#' @param RGB If RGB matrices are present in the object, choose between R, G, or B.
#' @param grey Boolean, if TRUE, then tlc is greyscaled. Default value is FALSE.
#' @param ... Additional graphics parameters.
#'
#' @return None.
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' # Tests 3D plot of the entire matrix
#' fname01 <- system.file("extdata", "test025to100sp.tiff", package="qtlc")
#' testTLC <- createTLC(fname01, RGB=FALSE)
#' 
#' # now we'll imitate interactive spot2D function,
#' # and create spots coordinates automatically,
#' # for interactive version run testTLC <- spot2D(testTLC)
#' testTLC$spots$x <- c(40.93354, 83.18687, 121.59899, 160.01111, 203.54485,
#'                      239.39616, 280.36909, 320.06161, 362.31494, 399.44666,
#'                      439.13919, 480.11211, 518.52423, 559.49716, 599.18969)
#' testTLC$spots$y <- c(198.3160, 198.3160, 199.2833, 198.3160, 198.3160,
#'                      198.3160, 198.3160, 198.3160, 197.3487, 198.3160,
#'                      199.2833, 198.3160, 199.2833, 199.2833, 199.2833)
#'
#' testTLC <- select2D(testTLC, 30, 30)
#' testTLC <- matrices2D(testTLC)
#' testTLC <- summat2D(testTLC)
#'
#' # 3D without OpenGL, shows only spot 13
#' showtlc3D(testTLC, spot=13, ogl=FALSE, grey=FALSE)
#' # without openGL and greyscaled
#' showtlc3D(testTLC, spot=13, ogl=FALSE, grey=TRUE)
#' #openGL showtime
#' showtlc3D(testTLC, spot=13, ogl=TRUE)
#'
#' @export
#' @importFrom plot3D persp3D
#' @importFrom grDevices gray.colors rainbow
#' @importFrom rgl open3d persp3d
#'
showtlc3D <- function(object, spot=NULL, ogl=FALSE, RGB=NULL, grey=FALSE, ...) {

    if (is.null(spot) && is.null(RGB) && is.na(dim(object$mat)[3])) {s3D(object$mat, ogl, grey, ...);}
    else if (is.null(RGB) && !is.na((dim(object$mat)[3] == 3))) {cat("\nR,G,B channel must be specified!\n"); return(FALSE);}
    else if(!is.null(RGB)) {
	if (RGB == "R") {s3D(object$mat[,,1], ogl, grey, ...);}
	else if (RGB == "G") {s3D(object$mat[,,2], ogl, grey, ...);}
	else if (RGB == "B") {s3D(object$mat[,,3], ogl, grey, ...);}
	else {cat("\nMethod cannot interprate the matrix! Please verify data.\n"); return(FALSE);}}
    
	if (!is.null(spot) && !is.null(object$spots) && !is.null(object$spot_matrices)) {
		s3D(object$spot_matrices[,,spot], ogl, grey, ...);
		return(TRUE);
		} else if (!is.null(spot) && is.null(object$spots) || is.null(object$spot_matrices)) {
			cat("\nErr: Spots must be selected. Matrices must be formed.\n");
			return(FALSE);
			}
	
	return(TRUE);
	}

