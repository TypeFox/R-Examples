#' Quick Map
#' 
#' \code{quick.map} creates and displays an image, identifying the threshold as
#' the background color, and create the gradient legend in the map.
#' 
#' An image is created of the map requested. A gradient legend
#' (\code{\link{legend.gradient}}) will be added if \code{pnts} (the position
#' of the legend) is specified.
#' 
#' @param sdm.asc an object of class 'asc' (adehabitat package), 'RasterLayer'
#' (raster package) or 'SpatialGridDataFrame' (sp package)
#' @param threshold to indicate the threshold limit of \code{sdm.asc}
#' @param bkgd.col to specify the background color
#' @param cols a set of 2 or more colors to be used in the image and the
#' gradient legend
#' @param zlim to specify the upper an lower limits, which are going to be the
#' labels of the gradient legend
#' @param pnts location information for adding the \code{legend.gradient}
#' @param ... other graphical parameters defined by image() or plot()
#' @return Nothing is returned, an image is created.
#' @author Lorena Falconi \email{lorefalconi@@gmail.com}
#' @examples
#' 
#' 
#' #create a matrix
#' tmat = { matrix(c( 0,0,0,1,0,0,1,1,0,1,
#'                    0,0,1,0,1,0,0,0,0,0,
#'                    0,1,NA,1,0,1,0,0,0,1,
#'                    1,0,1,1,1,0,1,0,0,1,
#'                    0,1,0,1,0,1,0,0,0,1,
#'                    0,0,1,0,1,0,0,1,1,0,
#'                    1,0,0,1,0,0,1,0,0,0,
#'                    0,1,0,0,0,1,0,0,0,1,
#'                    0,0,1,1,1,0,0,1,1,1,
#'                    1,1,1,0,0,0,0,1,1,1),nr=10,byrow=TRUE) }
#' 
#' #do the connected component labeling
#' tasc = ConnCompLabel(tmat)
#' 
#' #put in the gradient scale
#' pnts = cbind(x =c(1.1,1.2,1.2,1.1), y =c(0.9,0.9,0.7,0.7))
#' 
#' # Set the map and gradient leyend colors 
#' tasc.col=colorRampPalette(c("yellow","orange", "red"))(5)
#' 
#' #Create an image with the gradient legend
#' quick.map(tasc,0.09,bkgd.col = 'darkgrey', cols=tasc.col,
#'     axes=FALSE, xlim=c(0.0,1.35))
#' 
#' #########################
#' # Create an image with two colors: below the threshold and 
#' # above the threshold
#' 
#' # The next version of SDM Tools will let you create the legend.gradient  
#' # at a specific side of your image, and the user would not need to set 
#' # the coordinates.
#' 
#' # To create the legend.gradient at the bottom left of your image without  
#' # setting up the coordinates at the image you can do this:
#' 
#' xlim = c(-0.5,1)
#' ylim = c(0,1)
#' wid = diff(xlim)*0.05
#' ht = diff(ylim)*0.1
#' xvals = c(xlim[1]+wid,xlim[1]+2*wid,xlim[1]+2*wid,xlim[1]+wid)
#' yvals = c(ylim[1]+ht,ylim[1]+ht,ylim[1]+2*ht,ylim[1]+2*ht)
#' 
#' #Create the points for the legend.gradient
#' pnts=(cbind(xvals,yvals))
#' 
#' # Set the images colors: above the threshold is black and 
#' # below the threshold is darkgrey.
#' quick.map(tasc,0.09,bkgd.col = 'darkgrey', cols="black",
#'     axes=FALSE, xlim=c(-0.8, 1))
#' 
#' 
#' @export 
quick.map= function(sdm.asc, threshold, bkgd.col = 'grey',cols=heat.colors(100), zlim=NULL, pnts=NULL, ...){
	#check if raster from sp or raster package and convert if necessary
	if (any(class(sdm.asc) %in% 'RasterLayer')) sdm.asc = asc.from.raster(sdm.asc)
	if (any(class(sdm.asc) == 'SpatialGridDataFrame')) sdm.asc = asc.from.sp(sdm.asc)

	#get range of data
	trange=range(sdm.asc,na.rm=T)
	#make sure threshold is within range
	if (!(trange[1]< threshold & trange[2]>threshold)) stop('Cannot be 0')
	#rework sdm.asc
	sdm.asc[which(!is.na(sdm.asc) & sdm.asc<=threshold)] = 0
	tvals = seq(threshold,trange[2],(trange[2]-threshold)/(length(cols)+2))
	for (i in 1:length(cols)) sdm.asc[which(!is.na(sdm.asc) & sdm.asc<tvals[i] & sdm.asc>=tvals[i+1])] = 0
	#create the image
	if (is.null(zlim)) { lim = range(sdm.asc,na.rm=T) } else { lim = zlim }
	image(sdm.asc, col=c(bkgd.col,cols), xlab="", ylab="",zlim=lim,...)
	#add the scale legend
	if (!is.null(pnts)) legend.gradient(pnts,cols=c(bkgd.col,cols),limits=lim)
}
 
