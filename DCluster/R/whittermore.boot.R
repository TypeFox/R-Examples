whittermore.boot<-function(data, i, ...)
{
	data$Observed<-data$Observed[i]

	return( whittermore.stat(data, ...) ) 
}
