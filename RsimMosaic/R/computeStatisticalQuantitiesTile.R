computeStatisticalQuantitiesTile <-
function(img, useGradients=FALSE) {
#
# A function to compute the relevant statistical quantities (at this moment only the median value)
# of the RGB colors for the entire image. Optionally it can also compute the median values of the 
# RGB colors at the Le, Le, U, UR, R, LoR, Lo, and LoLe, corners of the image
# 
#
# img :: an image matrix as outputed by the readJPEG function (from the jpeg library)
# useGradients :: a flag indicating if the values of the median RGB colors should be calculated
# 
		xHalfSize <- dim(img[,,1])[1] / 2
		yHalfSize <- dim(img[,,1])[2] / 2

		if(!useGradients) {
			rgbStatWGrad <- c(median(as.vector(img[,,1])),
		 	median(as.vector(img[,,2])),
		 	median(as.vector(img[,,3]))
		 	)
		} else {
	 	 	rgbStatWGrad <- c(median(as.vector(img[,,1])),
		 	median(as.vector(img[,,2])),
		 	median(as.vector(img[,,3])),

			median(as.vector(img[(1:xHalfSize),,1])),
 		 	median(as.vector(img[(1:xHalfSize),,2])),
 		 	median(as.vector(img[(1:xHalfSize),,3])),

			median(as.vector(img[(1:xHalfSize),(1:yHalfSize),1])),
 	 		median(as.vector(img[(1:xHalfSize),(1:yHalfSize),2])),
 	 		median(as.vector(img[(1:xHalfSize),(1:yHalfSize),3])),

			median(as.vector(img[,(1:yHalfSize),1])),
 	 		median(as.vector(img[,(1:yHalfSize),2])),
 	 		median(as.vector(img[,(1:yHalfSize),3])),

			median(as.vector(img[(xHalfSize:(2*xHalfSize-1)),(1:yHalfSize),1])),
 	 		median(as.vector(img[(xHalfSize:(2*xHalfSize-1)),(1:yHalfSize),2])),
 	 		median(as.vector(img[(xHalfSize:(2*xHalfSize-1)),(1:yHalfSize),3])),

			median(as.vector(img[(xHalfSize:(2*xHalfSize-1)),,1])),
 	 		median(as.vector(img[(xHalfSize:(2*xHalfSize-1)),,2])),
 	 		median(as.vector(img[(xHalfSize:(2*xHalfSize-1)),,3])),

			median(as.vector(img[(xHalfSize:(2*xHalfSize-1)),(yHalfSize:(2*yHalfSize-1)),1])),
 		 	median(as.vector(img[(xHalfSize:(2*xHalfSize-1)),(yHalfSize:(2*yHalfSize-1)),2])),
 		 	median(as.vector(img[(xHalfSize:(2*xHalfSize-1)),(yHalfSize:(2*yHalfSize-1)),3])),

			median(as.vector(img[,(yHalfSize:(2*yHalfSize-1)),1])),
 	 		median(as.vector(img[,(yHalfSize:(2*yHalfSize-1)),2])),
 	 		median(as.vector(img[,(yHalfSize:(2*yHalfSize-1)),3])),

			median(as.vector(img[(1:xHalfSize),(yHalfSize:(2*yHalfSize-1)),1])),
 	 		median(as.vector(img[(1:xHalfSize),(yHalfSize:(2*yHalfSize-1)),2])),
 	 		median(as.vector(img[(1:xHalfSize),(yHalfSize:(2*yHalfSize-1)),3]))
 	 		)
		}

 	 	return(rgbStatWGrad)
}
