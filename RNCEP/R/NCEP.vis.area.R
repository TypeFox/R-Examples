NCEP.vis.area <- function(wx.data, layer=1, show.pts=TRUE, draw.contours=TRUE, cols=heat.colors(64), transparency=.5, axis.args=NULL,
	map.args=NULL, grid.args=NULL, title.args=NULL, interp.loess.args=NULL, image.plot.args=NULL, contour.args=NULL, points.args=NULL){

## Load the needed libraries ##
#importFrom(maps,map)
#importFrom(fields,image.plot)
#importFrom(tgp,interp.loess)
#importFrom(graphics,contour)
#require(maps);require(fields);require(tgp)

## Make sure the weather data are in an array and not just a matrix ##
if (length(dimnames(wx.data)) == 2) { wx.data <- array(wx.data,dim=c(dim(wx.data),1), dimnames=c(dimnames(wx.data),1)) }

## Determine which layer to use, if not specified numerically ##
if(is.numeric(layer) == FALSE){
	if(any(gsub(x=dimnames(wx.data)[[3]], pattern = "_", replacement = " ") == gsub(x=layer, pattern = "-", replacement = " ")) == FALSE){ stop("None of the layers in the weather dataset matches the layer specified") }

	layer <- which(gsub(x=dimnames(wx.data)[[3]], pattern = "_", replacement = " ") == gsub(x=layer, pattern = "-", replacement = " "))
}

## Specify the values of longitude ##
if(mean(as.numeric(dimnames(wx.data)[[2]])) >= 180){
    longitudes <- as.numeric(dimnames(wx.data)[[2]])-360 
	} else {
	longitudes <- as.numeric(dimnames(wx.data)[[2]])
	}

## First create a map that spans the range of the data ##
orig.map.args <- map.args
if(is.null(map.args)){
	map.args <- list(database = "world", xlim = range(longitudes), ylim=range(as.numeric(dimnames(wx.data)[[1]])), fill=TRUE, col='gray', mar=c(4.1, 9, par("mar")[3], 2.1)) } else {
	map.args <- c(map.args, 
		if(is.null(map.args$col)) {list(col='gray')}, 
		if(is.null(map.args$database)) {list(database="world")}, 
		if(is.null(map.args$xlim)) {list(xlim=range(longitudes))}, 
		if(is.null(map.args$ylim)) {list(ylim=range(as.numeric(dimnames(wx.data)[[1]])))},
		if(is.null(map.args$fill)) {list(fill=TRUE)},
		if(is.null(map.args$mar)) {list(mar=c(4.1, 9, par("mar")[3], 2.1))}) }
## Make sure that the area specified contains some land mass ##
trying.map <- try(do.call('map', map.args), silent=TRUE)

## If no land masses are present in the map, must produce a regular plot instead ##
if(length(trying.map) != 4){
	dev.off()
	if(is.null(orig.map.args$mar)) {par(mar=c(5, 4, 4, 6) + 0.1)} else {par(mar=orig.map.args$mar)}
	if(!is.null(orig.map.args$bg)) {par(bg=orig.map.args$bg)}
	alt.map.args <- list(x=1, y=1, type='n', xlim=map.args$xlim, ylim=map.args$ylim, axes=FALSE, xlab='', ylab='')	
do.call('plot', alt.map.args)
} 

## Add axis labels ##
axs1 <- do.call('axis', c(axis.args, list(side=1)))
axs2 <- do.call('axis', c(axis.args, list(side=2)))

## Add the grid lines ##
if(is.null(grid.args)){
	grid.args <- list(lty=2, col='black') } else {
	grid.args <- c(grid.args, 
		if(is.null(grid.args$lty)) {list(lty=2)},
		if(is.null(grid.args$col)) {list(col='black')}) }
do.call('abline', c(grid.args, if(is.null(grid.args$h)) {list(h=axs2)}, if(is.null(grid.args$v)) {list(v=axs1)}))

## Add a title to the figure ##
if(is.null(title.args)){
	title.args <- list(main=dimnames(wx.data)[[3]][layer], xlab="Longitude", ylab="Latitude", cex.lab=1.25) } else {
	title.args <- c(title.args, if(is.null(title.args$main)) {list(main=dimnames(wx.data)[[3]][layer])},
			if(is.null(title.args$xlab)) {list(xlab="Longitude")},
			if(is.null(title.args$ylab)) {list(ylab="Latitude")},
			if(is.null(title.args$cex.lab)) {list(cex.lab=1.25)}) }
do.call('title', title.args)


## Overlay the filled contour describing the weather data ##
## First interpolate the original grid a bit ##
if(is.null(interp.loess.args)){
	interp.loess.args <- list(x=rep(longitudes, each=length(dimnames(wx.data)[[1]])),
			y=rep(as.numeric(dimnames(wx.data)[[1]]), length(longitudes)), 
			z=as.vector(wx.data[,,layer]), span=0.6)} else {
	interp.loess.args <- c(interp.loess.args, 
			if(is.null(interp.loess.args$x)) {list(x=rep(longitudes, each=length(dimnames(wx.data)[[1]])))},
			if(is.null(interp.loess.args$y)) {list(y=rep(as.numeric(dimnames(wx.data)[[1]]), length(longitudes)))},
			if(is.null(interp.loess.args$z)) {list(z=as.vector(wx.data[,,layer]))},
			if(is.null(interp.loess.args$span)) {list(span=0.6)}) }
tst <- do.call('interp.loess', interp.loess.args)

## Calculate transparency in hex notation ##
transparency <- substr(rgb(1,1,1,transparency), start=8, stop=9)

## Specify the colors to use and apply transparency ##
cols <- paste(substr(cols, start=1, stop=7), transparency, sep='')

## Then add the interpolated grid to the map
if(is.null(image.plot.args)){
	image.plot.args <- list(tst, add=TRUE, col=cols) } else {
	image.plot.args <- c(image.plot.args,
			if(is.null(image.plot.args$x)) {list(x=tst$x)},
			if(is.null(image.plot.args$y)) {list(y=tst$y)},
			if(is.null(image.plot.args$z)) {list(z=tst$z)},
			if(is.null(image.plot.args$add)) {list(add=TRUE)},
			if(is.null(image.plot.args$col)) {list(col=cols)}) }
	image.plot.args$add <- TRUE
	image.plot.args$col <- cols
	do.call('image.plot', image.plot.args)
par(xpd=TRUE)

## If desired, add contour lines ##
if(draw.contours == TRUE){
if(is.null(contour.args)){
	contour.args <- list(x=tst$x, y=tst$y, z=tst$z, add=TRUE, labcex=1) } else {
	contour.args <- c(contour.args, 
			if(is.null(contour.args$x)) {list(x=tst$x)},
			if(is.null(contour.args$y)) {list(y=tst$y)},
			if(is.null(contour.args$z)) {list(z=tst$z)}, 
			if(is.null(contour.args$add)) {list(add=TRUE)},
			if(is.null(contour.args$labcex)) {list(labcex=1)}) }
	contour.args$add=TRUE
do.call('contour', contour.args)
}

## Make sure that there is a box drawn around the plot ##
rect(xleft=par()$usr[1], ybottom=par()$usr[3], xright=par()$usr[2], ytop=par()$usr[4])

## If desired, add points to indicate grid points in the dataset ##
if(show.pts==TRUE) {
if(is.null(points.args)){
	points.args <- list(x=rep(longitudes, length(dimnames(wx.data)[[1]])), 
			y=rep(as.numeric(dimnames(wx.data)[[1]]), each=length(dimnames(wx.data)[[2]]))) } else {
	points.args <- c(points.args,
			if(is.null(points.args$x)) {list(x=rep(longitudes, length(dimnames(wx.data)[[1]])))},
			if(is.null(points.args$y)) {list(y=rep(as.numeric(dimnames(wx.data)[[1]]), each=length(dimnames(wx.data)[[2]])))}) }
do.call('points', points.args)
}			

par(xpd=FALSE)
## End function ##
}

