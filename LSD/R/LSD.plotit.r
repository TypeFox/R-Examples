

### plotit ###


#' @export
#' @name plotit
#' @aliases LSD.plotit
#' @title Plotting wrapper function to plot plots in printable quality and all kinds of formats
#' @description Plotting wrapper function to save plots in R as "pdf", "ps", "jpeg", "png", "bmp" or "tiff".
#' @param filename name of the plot to be saved with the format type suffix.
#' @param sw scaling factor of weight.
#' @param sh scaling factor of height.
#' @param sres scaling factor of the resolution.
#' @param plotsfkt list of plots to be plotted.
#' @param ww width of window.
#' @param wh height of window.
#' @param pointsize the default pointsize of plotted text, interpreted as big points (1/72 inch) for plots to be saved.
#' @param dev.pointsize pointsize of plotted text, interpreted as big points (1/72 inch) for display in R.
#' @param paper needed only if filformat = "pdf" or "ps".
#' @param quality needed only if filformat = "jpeg".
#' @param units needed only if filformat = "jpeg", "png", "bmp" or "tiff".
#' @param bg backgroundcolor.
#' @param fileformat save the plot as "pdf", "ps", "jpeg", "png", "bmp" or "tiff".
#' @param saveit should plot be saved.
#' @param notinR should plot be not plotted in R.
#' @param addformat should plot be saved additionally in another format ("pdf", "ps", "jpeg", "png", "bmp" or "tiff").
#' @author Bjoern Schwalb
#' @seealso \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples data(homer)
#'
#' plotsfkt = function(){
#' 		colpal = c("white","black","yellow","wheat3")
#' 		align(homer,colpal = colpal,main = "D'OH!",asp = 1,axes = FALSE)
#' }
#'
#' plotit("homer",sw = 2,sh = 2,sres = 2,plotsfkt,saveit = TRUE,fileformat = "jpeg")
#'
#' plotit("homer",sw = 2,sh = 2,sres = 2,plotsfkt,saveit = TRUE,fileformat = "png")
#'
#' plotit("homer",sw = 2,sh = 2,sres = 2,plotsfkt,saveit = TRUE,fileformat = "bmp")
#'
#' plotit("homer",sw = 2,sh = 2,sres = 2,plotsfkt,saveit = TRUE,fileformat = "tiff")
#'
#' plotit("homer",sw = 2,sh = 2,sres = 2,plotsfkt,saveit = TRUE,fileformat = "ps")
#'
#' plotit("homer",sw = 2,sh = 2,sres = 2,plotsfkt,saveit = TRUE,fileformat = "pdf")
#' @keywords plot, jpeg, png, bmp, tiff, ps, pdf


plotit = function(filename,sw = 1,sh = 1,sres = 1,plotsfkt,ww = 7,wh = 7,pointsize = 12,dev.pointsize = 8,paper = "special",quality = 100,units = "px",bg = "white",fileformat = "jpeg",saveit = FALSE,notinR = FALSE,addformat = NULL)
{
	
	# switch between different file formats for plot saving (if saveit = TRUE) #
	
	pwidth = sw*480
	pheight = sh*480
	pres = sres*72
	if (saveit){
		switch(fileformat,"jpeg" = jpeg(filename = paste(filename,".",fileformat,sep=""), width = pwidth, height = pheight,units = units, pointsize = pointsize, quality = quality, bg = bg,res = pres),
				"png" = png(filename = paste(filename,".",fileformat,sep=""), width = pwidth, height = pheight,units = units, pointsize = pointsize, bg = bg,res = pres),
				"bmp" = bmp(filename = paste(filename,".",fileformat,sep=""), width = pwidth, height = pheight,units = units, pointsize = pointsize, bg = bg,res = pres),
				"tiff" = tiff(filename = paste(filename,".",fileformat,sep=""), width = pwidth, height = pheight,units = units, pointsize = pointsize, bg = bg,res = pres),
				"ps" = postscript(file = paste(filename,".",fileformat,sep=""), width = ww, height = wh, pointsize = pointsize, paper = paper),
				"pdf" = pdf(file = paste(filename,".",fileformat,sep=""), width = ww, height = wh, pointsize = pointsize, paper = paper))
		plotsfkt()
		dev.off()
	}
	
	# switch between different file formats for additional plot saving (if saveit = TRUE) #
	
	if (saveit){
		if (!is.null(addformat)){
			switch(addformat,"jpeg" = jpeg(filename = paste(filename,".",addformat,sep=""), width = pwidth, height = pheight,units = units, pointsize = pointsize, quality = quality, bg = bg,res = pres),
					"png" = png(filename = paste(filename,".",addformat,sep=""), width = pwidth, height = pheight,units = units, pointsize = pointsize, bg = bg,res = pres),
					"bmp" = bmp(filename = paste(filename,".",addformat,sep=""), width = pwidth, height = pheight,units = units, pointsize = pointsize, bg = bg,res = pres),
					"tiff" = tiff(filename = paste(filename,".",addformat,sep=""), width = pwidth, height = pheight,units = units, pointsize = pointsize, bg = bg,res = pres),
					"ps" = postscript(file = paste(filename,".",addformat,sep=""), width = ww, height = wh, pointsize = pointsize, paper = paper),
					"pdf" = pdf(file = paste(filename,".",addformat,sep=""), width = ww, height = wh, pointsize = pointsize, paper = paper))
			plotsfkt()
			dev.off()
		}
	}
	
	# should plotting be repressed in R (if notinR = TRUE) #
	
	if (!notinR){
		dev.new(width = ww,height = wh,pointsize = dev.pointsize)
		plotsfkt()
	}
}


### alias ###


LSD.plotit = plotit



