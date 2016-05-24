stdev.f <-
function(x)
{
	## computes the standard deviation of x using the Fortran interface
	return(.Fortran("rs_stdev",s=as.double(0),
	                as.double(x),as.double(mean.f(x)),
	                as.integer(length(x)),package="stima")$s)
}
