tango.boot<-function(data, i, ...)
{
	data$Observed<-data$Observed[i]

	tango.stat(data, ...)
}
