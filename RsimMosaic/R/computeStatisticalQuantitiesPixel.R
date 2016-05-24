computeStatisticalQuantitiesPixel <-
function(i, j, img, useGradients=FALSE) {
#
# A function to compute the relevant statistical quantities for the pixel. Optionally it can also 
# output the values of relevant, nearby pixels as RGB colors at the Le, Le, U, UR, R, LoR, Lo, 
# and LoLe, corners of the image
#
# i :: index in i
# j :: index in j
# img :: an image matrix as outputed by the readJPEG function (from the jpeg library)
# useGradients :: a flag indicating if the values of the median RGB colors should be calculated
# 
	if(!useGradients){	
		pixelRGBandNeigArray <- c(
			img[i,j,1], img[i,j,2],img[i,j,3]
			)
	} else {
		pixelRGBandNeigArray <- c(
			img[i,j,1],     img[i,j,2],     img[i,j,3],
			img[i-1,j,1],   img[i-1,j,2],   img[i-1,j,3],
			img[i-1,j-1,1], img[i-1,j-1,2], img[i-1,j-1,3],
			img[i,j-1,1],   img[i,j-1,2],   img[i,j-1,3],
			img[i+1,j-1,1], img[i+1,j-1,2], img[i+1,j-1,3],
			img[i+1,j,1],   img[i+1,j,2],   img[i+1,j,3],
			img[i+1,j+1,1], img[i+1,j+1,2], img[i+1,j+1,3],
			img[i,j+1,1],   img[i,j+1,2],   img[i,j+1,3],
			img[i-1,j+1,1], img[i-1,j+1,2], img[i-1,j+1,3]
			)
	}
 	return(pixelRGBandNeigArray)
}
