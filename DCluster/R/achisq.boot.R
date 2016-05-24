achisq.boot<-function(data, i, ...)
{
	data$Observed<-data$Observed[i]
	return( achisq.stat(data, ...)$T )
}
