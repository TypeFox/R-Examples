resid.f <-
function(x)
{
	## computes the residual of vector x using the Fortran interface
	n1<-dim(x)[1]
	n2<-dim(x)[2]
	return(.Fortran("rs_resid",r=as.double(numeric(n1)),as.matrix(x),
	       as.integer(n1),as.integer(n2),package="stima")$r)
}
