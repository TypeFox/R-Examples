rsq.f <-
function(x)
{
	## computes the mean of vector x using the Fortran interface
	n1<-dim(x)[1]
	n2<-dim(x)[2]
	return(.Fortran("rs_rsq",r=as.double(0),as.matrix(x),
	       as.integer(n1),as.integer(n2),package="stima")$r)
}
