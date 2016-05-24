forecastfplsr <-
function(object, components = 2, h = 20)
{
	xname = object$xname
	yname = object$yname
	foreca = matrix(, length(object$x), h)
	data = object
	for(i in 1:h)
	{
		fore = fplsr(data, order = components)$Ypred
		foreca[,i] = fore$y
		newdata = cbind(data$y, fore$y)
		colnames(newdata) = as.numeric(colnames(data$y))[1]:(max(as.numeric(colnames(data$y))) + 1)
		data = fts(data$x, newdata)
	}
	colnames(foreca) = (max(as.numeric(colnames(object$y)))+1):(max(as.numeric(colnames(object$y)))+h)
	forecafts = fts(object$x, foreca, xname = xname, yname = yname)
	return(forecafts)
}

