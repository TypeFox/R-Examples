`modelimage` <-
function(model, color.palette = topo.colors, log = TRUE, plottitle = "RDE Model Selection", ...)
{
	# find out whether there's only one parameter
	# to name x-axis
	if(length(model$params) == 1)
	{
		xlab <- names(model$params)[1]
		x <- model$params[[1]]
	}
	else
	{
		xlab <- "kernel parameters"
		x <- 1:nrow(model$errs)
	}
	
	# scale logarithmic?
	if(log)
	{
		x <- log10(x)
		xlab <- paste("log10(", xlab, ")", sep = "")
	}
	
	ylab <- "dimension"
	y <- 1:ncol(model$errs)
	
	# draw filled contour
	filled.contour(x = x, y = y, z = model$errs,
		plot.title = title(main = plottitle, xlab = xlab, ylab = ylab),
		color.palette = color.palette,
		plot.axes = 	{
					axis(1); axis(2); 
					lines(x = x, y = model$rds, lwd = 2);
				}, ...)
}

