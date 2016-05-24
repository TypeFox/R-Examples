pottwhitt.boot<-function(data, i)
{
	data$Observed<-data$Observed[i]

	return( pottwhitt.stat(data)$T )

}
