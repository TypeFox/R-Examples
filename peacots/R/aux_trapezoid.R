aux_trapezoid <-
function(x,y){
	N = length(x);
	return(sum((x[2:N]-x[1:N-1])*(y[2:N]+y[1:N-1])*0.5));
}
