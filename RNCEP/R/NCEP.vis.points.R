NCEP.vis.points <- function(wx, lats, lons, cols=heat.colors(64), transparency=.5, connect=TRUE,
	axis.args=NULL, points.args=NULL, map.args=NULL, grid.args=NULL, title.args=NULL, image.plot.args=NULL, lines.args=NULL){	

## Make sure that the lengths of the input variables match ##
if(length(wx) != length(lats) | length(wx) != length(lons)) {stop("length of 'wx' different from lengths of 'lats' and/or 'lons'")}

## Load the needed libraries ##
#importFrom(maps,map)
#importFrom(fields,image.plot)
#require(maps); require(tgp)

## Determine the range to show on the map ##
lon.lim <- c(min(range(lons, na.rm=TRUE), na.rm=TRUE) - (diff(range(lons, na.rm=TRUE), na.rm=TRUE)*.1), 
	max(range(lons, na.rm=TRUE), na.rm=TRUE) + (diff(range(lons, na.rm=TRUE), na.rm=TRUE)*.1))
lat.lim <- c(min(range(lats, na.rm=TRUE), na.rm=TRUE) - (diff(range(lats, na.rm=TRUE), na.rm=TRUE)*.1), 
	max(range(lats, na.rm=TRUE), na.rm=TRUE) + (diff(range(lats, na.rm=TRUE), na.rm=TRUE)*.1))
## Determine their spatial range
rlon <- diff(lon.lim)
rlat <- diff(lat.lim)	
## Make sure that the range of lats and longs is large enough ##
	if(rlon < (rlat*.6)){
		xtra <- diff(c(rlon, (rlat*.6)))/2
		lon.lim <- c(lon.lim[1] - xtra, lon.lim[2] + xtra)
		}
	if(rlat < (rlon*.6)){
		xtra <- diff(c(rlat, (rlon*.6)))/2
		lat.lim <- c(lat.lim[1] - xtra, lat.lim[2] + xtra)
		}

## Make sure that the longitudes are specified correctly ##
if(mean(lon.lim) >= 180){
    lon.lim <- lon.lim-360
	lons <- lons-360
	}
if(mean(lon.lim) <= -180){
	lon.lim <- lon.lim+360
	lons <- lons+360
	}
	
## Specify the full range of colors to be considered ##
all.cols <- cols
	
## First create a map that spans the range of the data ##
orig.map.args <- map.args
if(is.null(map.args)){
	map.args <- list(database = "world", xlim=lon.lim, ylim=lat.lim, fill=TRUE, col='gray', mar=c(4.1, 9, par("mar")[3], 0.1)) } else {
	map.args <- c(map.args, 
		if(is.null(map.args$col)) {list(col='gray')}, 
		if(is.null(map.args$database)) {list(database="world")}, 
		if(is.null(map.args$xlim)) {list(xlim=lon.lim)}, 
		if(is.null(map.args$ylim)) {list(ylim=lat.lim)},
		if(is.null(map.args$fill)) {list(fill=TRUE)},
		if(is.null(map.args$mar)) {list(mar=c(4.1, 9, par("mar")[3], 0.1))}) }
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

## Make sure that there is a border drawn around the plot ##
rect(xleft=par()$usr[1], ybottom=par()$usr[3], xright=par()$usr[2], ytop=par()$usr[4])

## Add a title to the plot ##
if(is.null(title.args)){
	title.args <- list(main='Weather interpolated to points', xlab="Longitudes", ylab="Latitudes", cex.lab=1.25) } else {
	title.args <- c(title.args, if(is.null(title.args$main)) {list(main='Weather interpolated to points')},
			if(is.null(title.args$xlab)) {list(xlab="Longitudes")},
			if(is.null(title.args$ylab)) {list(ylab="Latitudes")},
			if(is.null(title.args$cex.lab)) {list(cex.lab=1.25)}) }
do.call('title', title.args)

## Optionally, draw a line connecting the points ##
if(connect == TRUE){
	if(is.null(lines.args)) { 
		lines.args <- list(x=lons, y=lats, lwd=2) } else {
		lines.args <- c(lines.args, if(is.null(lines.args$x)) {list(x=lons)},
				if(is.null(lines.args$y)) {list(y=lats)},
				if(is.null(lines.args$lwd)) {list(lwd=2)}) }
do.call('lines', lines.args)
}

## Calculate transparency in hex notation ##
transparency <- substr(rgb(1,1,1,transparency), start=8, stop=9)

## Specify the colors to use and apply transparency ##
cols <- paste(substr(cols, start=1, stop=7), transparency, sep='')

## Assign colors to the points based on the value of the weather variable at that point ##
brks <- seq(floor(min(wx, na.rm=TRUE)), ceiling(max(wx, na.rm=TRUE)), length.out=length(cols))
int <- findInterval(wx, brks, all.inside=TRUE)
cols <- cols[int]

## Then add the points to the map
if(is.null(points.args)){
	points.args <- list(x=lons, y=lats, col='black', cex=2, bg=cols, pch=21) } else {
	points.args <- c(points.args,
			if(is.null(points.args$x)) {list(x=lons)},
			if(is.null(points.args$y)) {list(y=lats)},
			if(is.null(points.args$col)) {list(col='black')},
			if(is.null(points.args$bg)) {list(bg=cols)},
			if(is.null(points.args$cex)) {list(cex=2)},
			if(is.null(points.args$pch)) {list(pch=21)}) }
		points.args$bg <- cols
do.call('points', points.args)		

## Add a colorbar style legend ##
if(is.null(image.plot.args)){
	image.plot.args <- list(x=lons, y=lats, z=wx, legend.only=TRUE, col=all.cols) } else {
	image.plot.args <- c(image.plot.args,
			if(is.null(image.plot.args$x)) {list(x=lons)},
			if(is.null(image.plot.args$y)) {list(y=lats)},
			if(is.null(image.plot.args$z)) {list(z=wx)},
			if(is.null(image.plot.args$col)) {list(col=all.cols)}) }
		image.plot.args$legend.only=TRUE
	do.call('image.plot', image.plot.args)

## End function ##
}
