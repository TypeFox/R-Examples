#DEPRECATED FUNCTION (works but not used anymore)
interpolateWithinGrid2D <-
function(grid_string, x, y, default_z){
	# create grid
	parts 		= strsplit(grid_string,split=",",fixed=TRUE);
	x_values 	= as.numeric(strsplit(parts[[1]][1],split=" ",fixed=TRUE)[[1]]);
	y_values 	= as.numeric(strsplit(parts[[1]][2],split=" ",fixed=TRUE)[[1]]);
	z_values 	= as.numeric(strsplit(parts[[1]][3],split=" ",fixed=TRUE)[[1]]);
	NX 			= length(x_values);
	NY 			= length(y_values);
	NZ			= length(z_values);
	if(length(z_values)!=NX*NY){
		cat(sprintf("ERROR parsing grid string: NX=%d, NY=%d but NZ=%d != NX*NY=%d\n",NX,NY,NZ,NX*NY));
		return(0);
	}
	EPSILON		= 1e-5; # relative grid tolerance, mainly to circumvent rounding errors
	
	# find enclosing grid box
	if((!isInRange(x_values[1],x_values[NX],x,EPSILON)) || (!isInRange(y_values[1],y_values[NY],y,EPSILON))) return(default_z);
	x_du = x + abs(x)*EPSILON;
	x_dd = x - abs(x)*EPSILON;
	y_du = y + abs(y)*EPSILON;
	y_dd = y - abs(y)*EPSILON;
	for(xi in 2:NX){
		if(x_values[xi-1]<=x_du && x_values[xi]>=x_dd){ xi2 = xi; break; }
	}
	for(yi in 2:NY){
		if(y_values[yi-1]<=y_du && y_values[yi]>=y_dd){ yi2 = yi; break; }
	}
	x1 = x_values[xi2-1];
	x2 = x_values[xi2];
	y1 = y_values[yi2-1];
	y2 = y_values[yi2];
	z11 = z_values[(xi2-2)*NY + (yi2-1)];
	z12 = z_values[(xi2-2)*NY + yi2];
	z21 = z_values[(xi2-1)*NY + (yi2-1)];
	z22 = z_values[(xi2-1)*NY + yi2];
	
	# interpolate bilinearly
	tz1 = z11*(x2-x)/(x2-x1) + z21*(x-x1)/(x2-x1);
	tz2 = z12*(x2-x)/(x2-x1) + z22*(x-x1)/(x2-x1);
	z = tz1*(y2-y)/(y2-y1) + tz2*(y-y1)/(y2-y1);
	return(z);
}
