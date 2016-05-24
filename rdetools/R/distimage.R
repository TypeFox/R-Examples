`distimage` <-
function(model, color.palette = terrain.colors, log = TRUE, plottitle = "Distance of Ys", ...)
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
		x <- 1:nrow(model$Yd)
	}
	
	# scale logarithmic?
	if(log)
	{
		x <- log10(x)
		xlab <- paste("log10(", xlab, ")", sep = "")
	}
	
	ylab <- "dimension"
	y <- 1:ncol(model$Yd)
	
	# draw filled contour
	filled.contour(x = x, y = y, z = model$Yd, 
		plot.title = title(main = plottitle, xlab = xlab, ylab = ylab),
		color.palette = color.palette,
		plot.axes = 	{
					axis(1); axis(2); 
					lines(x = x, y = model$rds, lwd = 2);
				}, ...)
}

