isInRange <-
function(xmin,xmax,x,epsilon){
	return ((x+abs(x)*epsilon>=xmin) && (x-abs(x)*epsilon<=xmax));
}
