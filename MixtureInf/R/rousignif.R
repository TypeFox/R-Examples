rousignif <-
function(x)
#This function return the output combining round and signif functions.
#x:	a single numeric or a numeric vector or matrix.
{
	a=NULL
	if (is.matrix(x))
	{
		a=c(nrow(x),ncol(x))
		b1=rownames(x)
		b2=colnames(x)
		x=as.vector(x)
	}
	x[abs(x)>=0.001]=round(x[abs(x)>0.001],3)
	x[abs(x)<=0.001]=signif(x[abs(x)<=0.001],3)
	if (is.null(a)==F)
	{
		x=matrix(x,a)
		rownames(x)=b1
		colnames(x)=b2
	}
	x
}
