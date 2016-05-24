stone.boot<-function(data, i, ...)
{
	data$Observed<-data$Observed[i]

	return( stone.stat(data, ...)[1] )
}
