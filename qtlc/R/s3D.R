#' Internal function used by showtlc3D
#'
#' Internal function used by showtlc3D
#' @param mat Matrix with x,y,Intensity dimensions
#' @param ogl If TLC 3D plot use OpenGL library for fast and interactive 3D plot. (This functionality is based on the \code{rgl} package.) Otherwise the function is based on the plot3D package.
#' @param grey Boolean, if TRUE, then tlc is greyscaled. Default value is FALSE.
#' @param ... Additional graphics parameters.
#'
#' @return None.
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' #Internal function.
#' }
#'
#' @importFrom plot3D persp3D
#' @importFrom grDevices gray.colors rainbow
#' @importFrom rgl open3d persp3d
#'
s3D <- function(mat, ogl, grey, ...) {
	
    if (ogl == FALSE ) {
        if (grey == TRUE) { plot3D::persp3D(1:nrow(mat), 1:ncol(mat), mat, zlab="Intensity", col=rev(grDevices::gray.colors(100, start = 0.0, end = 1.0, gamma = 2.2, alpha = NULL), ...));}
        else {
		plot3D::persp3D(1:nrow(mat), 1:ncol(mat), mat, zlab="Intensity", ...);
		}}
	
	##OpenGL verzija poziva
	else if (ogl == TRUE) {
            if (grey == TRUE) {color <- rev(grDevices::gray.colors(25, start = 0.0, end = 1.0, gamma = 2.2, alpha = NULL));}
            else {color <- rev(grDevices::rainbow(25, start = 0/6, end = 4/6));}
	    zcolor <- cut(mat, 25);
	    rgl::open3d();		
	    rgl::persp3d(1:nrow(mat), 1:ncol(mat), mat, col=color[zcolor], box=F, axes=T, zlab="Intensity", xlab="x", ylab="y", ...);
		}
	}	

