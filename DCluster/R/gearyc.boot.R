gearyc.boot<-function(data, i, ...)
{
	data$Observed<-data$Observed[i]
	return( gearyc.stat(data, ...) )

}
