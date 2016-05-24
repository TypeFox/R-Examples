cor.f <-
function(x,y)
{
	n<-length(x)
	corr<-0
	if (n!=length(y))
		print("Length Error")
	else
		corr<-.Fortran("rs_cor",corr=as.double(0),
		               as.double(x),as.double(y),as.integer(n),package="stima")$corr
	return(corr)
}
