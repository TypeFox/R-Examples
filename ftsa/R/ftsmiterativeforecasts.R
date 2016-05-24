ftsmiterativeforecasts <-
function(object, components, iteration = 20)
{
	xname = object$xname
	yname = object$yname
	foreca = matrix(, length(object$x), iteration)
	data = object
	for(i in 1:iteration)
	{
		fore = forecast(ftsm(data, order = components), h = 1)$mean$y
		foreca[,i] = fore
		newdata = cbind(data$y, fore)
		colnames(newdata) = as.numeric(colnames(data$y))[1]:(max(as.numeric(colnames(data$y))) + 1)
		data = fts(data$x, newdata)
	}
	colnames(foreca) = (max(as.numeric(colnames(object$y)))+1):(max(as.numeric(colnames(object$y)))+iteration)
	forecafts = fts(object$x, foreca, xname = xname, yname = yname)
	return(forecafts)
}

