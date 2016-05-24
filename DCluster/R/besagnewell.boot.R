besagnewell.boot<-function(data, i, ...)
{
	data$Observed<-data$Observed[i]
	return( besagnewell.stat(data, ...) )
}
