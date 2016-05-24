car2sph <-
function(x,y,z,deg=TRUE){
	if(is.matrix(x) || is.data.frame(x)){
		if(ncol(x) == 1){x=x[,1]}
        else	if(ncol(x) == 2){y = x[, 2];x = x[, 1]}
        else	if(ncol(x) == 3){z = x[, 3];y = x[, 2];x = x[, 1]}
    }
    if(missing(x) | missing(y) | missing(z)){stop('Missing full cartesian 3D input data.')}
	radius = sqrt(x^2 + y^2 + z^2)
    long = atan2(y, x)
    lat = asin(z/radius)
	if(deg){
		long=long*180/pi
		lat=lat*180/pi
	}
	lat[radius==0]=0
	return=cbind(long=long,lat=lat,radius=radius)
	}
