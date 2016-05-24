#' Biplot Comparison of Matrices
#' 
#' \code{compare.matrix} compares the values within two matrices (e.g., ESRI
#' ArcInfo ASCII raster files) and produces a biplot that shows the frequency
#' of each data combination shared between the matrices. The plot is overlayed
#' with contour lines that demarcate parts of the the plot that share the same
#' frequency of data combinations. \cr \cr \bold{NOTE:} it is assumed the
#' matrices are of the same extent, cell size and scaled to be the same units.
#' 
#' 
#' @param x a matrix of data; the matrix can be a raster of class 'asc'
#' (adehabitat package), 'RasterLayer' (raster package) or
#' 'SpatialGridDataFrame' (sp package)
#' @param y a matrix of data of the same extent, cell size and class as 'x'
#' @param nbins number of equally spaced bins used to partition range of values
#' in 'x' & 'y'
#' @param ... other graphical parameters defined by image(), contour(), or
#' plot()
#' @return Nothing is returned but images are created.
#' @author Luke Shoo \email{luke.shoo@@jcu.edu.au}
#' @examples
#' 
#' 
#' #create some simple objects of class 'asc'
#' tasc = as.asc(matrix(rep(x=1:10, times=1000),nr=100)); print(tasc)
#' #modify the asc objects so that they are slightly different
#' tasc1 = tasc + runif(n = 10000, min = -1, max = 1)
#' tasc2 = tasc + rnorm(n = 10000, mean = 1, sd = 1)
#' 
#' #create some images
#' #basic plot showing the density of data combinations shared 
#' #by the two matrices
#' compare.matrix(tasc1,tasc2,20)
#' 
#' #same as previous but with data partioned amoung more bins
#' compare.matrix(tasc1,tasc2,50)
#' 
#' #same as previous but altering the number of contour levels 
#' #and adding more graphical functions
#' compare.matrix(tasc1,tasc2,50,nlevels=5, xlab='asc1',ylab='asc2',
#'   main='Comparison between asc and asc2', bg="grey")
#' 
#' 
#' @export compare.matrix
compare.matrix <-
function(x,y,nbins,...){
	#check if raster from sp or raster package and convert if necessary
	if (any(class(x) %in% 'RasterLayer')) { x = asc.from.raster(x); y = asc.from.raster(y) }
	if (any(class(x) == 'SpatialGridDataFrame')) { x = asc.from.sp(x); y = asc.from.sp(y) }
	#confirm same extents
	if(length(which(dim(x)==dim(y)))!=2) stop('matrix objects must be of the same extent')
	#select only positions where data is finite
	pos=which(is.finite(x))
	#create a dataframe where each asc file is a column of finite data
	tdata=data.frame(x=x[pos], y=y[pos])
	#get the range of the data and add a small number to maximum so it will be contained within a bin
	tt = range(c(tdata$x,tdata$y),na.rm=T); tt[2] = tt[2] + 1e-7
	#create the bin classes
	bins = seq(tt[1],tt[2],(tt[2]-tt[1])/nbins)
	mids = bins[1:(length(bins)-1)] + 0.5 * mean(diff(bins))
	tdata$xbin = cut(tdata$x,breaks=bins,labels=mids,include.lowest=T,right=T)
	tdata$ybin = cut(tdata$y,breaks=bins,labels=mids,include.lowest=T,right=T)
	#get frequencies (counts) of data combinations
	aggdata=aggregate(x=tdata$x, by=list(x=tdata$xbin, y=tdata$ybin), FUN=length)
	names(aggdata)[3] = 'counts'
	#create a data frame of all possible combinations of data within range
	tasc.p=expand.grid(x=mids, y=mids)  
	#merge tasc.p with aggdata
	mdata=merge(x=tasc.p, y=aggdata, all=TRUE)
	names(mdata)[3] = 'out'
	#convert to a matrix
	mdata <- matrix(mdata$out, nrow = length(mids), ncol=length(mids), byrow=TRUE)
	#create an image
	suppressWarnings(image(x=mids, y=mids, z=mdata, col=c(heat.colors(10)[10:1]),...))
	#overlay contours
	contour(x=mids, y=mids, z=mdata,  col="black", lty="solid", add=TRUE,...)
}

