# Kernel Densities from point object
#kde.points <- function(pts,h,n=200,lims=NULL) {
#	xy = coordinates(pts)
#	p4s = CRS(proj4string(pts))
# 	if (missing(h)) h = c(bandwidth.nrd(xy[,1]),bandwidth.nrd(xy[,2]))
# 	if (is.null(lims)) 
# 		{ lims=c(range(xy[,1]),range(xy[,2])) }
# 	else 
# 		{ lims = t(bbox(lims))}
# 	kd = kde2d(xy[,1],xy[,2],h,n,lims)
# 	temp = SpatialPoints(expand.grid(kd$x,kd$y))
# 	temp = SpatialPixelsDataFrame(temp,data.frame(kde=array(kd$z,length(kd$z))))
# 	proj4string(temp) = p4s
# 	temp}
# 
# 
# 
# 
# 
# level.plot <- function(grd,shades,index=1,add=FALSE) {
# 	x = unique(coordinates(grd)[,1])
# 	y = unique(coordinates(grd)[,2])
# 	z = data.frame(grd)[,index]
# 	zz = z[!is.na(z)]
# 	if(missing(shades)) shades = auto.shading(zz,cutter=range.cuts)
# 	levels = c(min(z,na.rm=TRUE),shades$breaks,max(z,na.rm=TRUE))
# 	z = matrix(z,length(x),length(y))
# 	storage.mode(z) = "double"
# 	if (!add) {
# 		xx = c(min(x),max(x))
# 		yy = c(min(y),max(y))
# 		plot(xx,yy,type='n',xaxt='n',yaxt='n',xaxs='i',yaxs='i',bty='n',xlab='',ylab='',asp=1) }
# 	cols = shades$cols
# 	.Internal(filledcontour(as.double(x),as.double(y),z,as.double(levels),col = cols))}


