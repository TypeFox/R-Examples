bilinearInterpolator <-
function(oldMatrix, pointsInNewX, pointsInNewY) {
# 
# A function to interpolate a matrix to a new resolution
# 
# oldMatrix :: the original matrix
# pointsInNewX :: the number of points in the new matrix (rows)
# pointsInNewY :: the number of points in the new matrix (columns)
# 
	return(fields::interp.surface.grid(list(x=seq(nrow(oldMatrix)), y=seq(ncol(oldMatrix)), z=oldMatrix),
            				   list(x=seq(1,nrow(oldMatrix), length=pointsInNewX), y=seq(1,ncol(oldMatrix), length=pointsInNewY)))$z)
}
