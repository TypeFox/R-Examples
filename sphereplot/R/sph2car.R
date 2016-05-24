sph2car <-
function(long,lat,radius=1,deg=TRUE){
    if(is.matrix(long) || is.data.frame(long)){
        if(ncol(long) == 1){long = long[,1]}
		else	if(ncol(long) == 2){lat = long[, 2];long = long[, 1]}
		else	if(ncol(long) == 3){radius = long[, 3];lat = long[, 2];long = long[, 1]}
		}
	if(missing(long) | missing(lat)){stop('Missing full spherical 3D input data.')}
    if(deg){
    	long = long * pi/180
		lat = lat * pi/180
	}
return=cbind(x=radius*cos(long)*cos(lat),y=radius*sin(long)*cos(lat),z=radius*sin(lat))
}
