NULL
#'
#' Graphic Representation of a Color bar, function written by John Colby
#' 
#' @param lut see reference \url{http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette} 
#' @param min see reference \url{http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette} 
#' @param max see reference \url{http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette} 
#' @param nticks see reference \url{http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette}
#' @param ticks see reference \url{http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette} 
#' @param title see reference \url{http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette} 
#' @param width,height width and height of the device 
#' @param digits specified number of significant digits
#' @param pdf character value for pdf output file.  Default is \code{NULL} and no pdf file is created. 
#' @param ncolmax maximum number of colors. Default is 100. 
#' 
#' @author John Colby \url{http://stackoverflow.com/users/412342/john-colby}
#' @export
#' @name color.bar
#' @note This function is taken from \url{http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette}. Plese visit the URL for major details and give your feedback if possible. 
#' @references \url{http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette}
#' @examples color.bar(colorRampPalette(c("light green", "yellow", "orange", "red"))(100), -1)
#' 
#' 
#' 

# OK color.bar.raster() <- function TO DO 
# Function to plot color bar 
#  http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='',width=1.75, height=5,ncolmax=100,digits=4,pdf=NULL) {
	
	
	if (length(lut)>ncolmax) {
		start_lut <- lut[1]
		end_lut <- lut[length(lut)]
		
		index <- as.integer((1:ncolmax)*length(lut)/ncolmax)
		lut <- lut[index] 
		
		lut[1] <- start_lut 
		lut[length(lut)] <- end_lut		
		
	}
	
	
	
	scale = (length(lut)-1)/(max-min)
	

	ticks <- signif(ticks,digits=digits)
	
	if (is.null(pdf)) 
		{dev.new(width=width, height=height)
	} else {
		
		pdf(pdf,width=width, height=height)
		
	}
	plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
	axis(2, ticks, las=1)
	for (i in 1:(length(lut)-1)) {
		y = (i-1)/scale + min
		y1 <- y+1/scale
#		y <- signif(y,digits=digits)
#		y1 <- signif(y1,digits=digits)
		rect(0,y,10,y1, col=lut[i], border=NA)
	}
	
	if (!is.null(pdf)) dev.off()
	
	
}
