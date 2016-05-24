# blend & blended.landscape methods

blend = function(x, warn = FALSE){
	UseMethod("blend")
}

blend.character = function(x, warn = FALSE){
	# import, then hand off to blend.list
	blend(import(x), warn = warn)	
}


blend.list = function(x, warn = FALSE){
	# sanitize, then hand off to blend.landscape.list
	blend(sanitize(x), warn = warn)
}


blend.landscape.list = function(x, warn = FALSE){
	# lapply blendOnce across landscapes and then bundle it up if needed
	unbundled = lapply(
		1:length(x$landscape.names),
		function(i) blendOnce(
				x$native.landscapes[[i]], 
				x$exotic.landscapes[[i]], 
				x$landscape.names[[i]],
				warn = warn
			)
	)
	if(length(unbundled) > 1){
		bundle(unbundled)
	}else{
		unbundled[[1]]
	}
}





##########################

# blended.landscape methods

print.blended.landscape = function(x, ...){
	print(x$results.table, ...)
}


plot.blended.landscape = function(
	x,
	main = paste(x$name, "homogenization vs. exotic occupancy"),
	pch = 16,
	cex = .6,
	cex.axis = 1.3,
	cex.lab = 1.3,
	cex.main = 1.6,
	lwd = 2,
	xaxs = "i",
	xlim = c(0, 1),
	ylab = expression(paste(Delta, " mean similarity")),
	...
){
	baseline.scipen = getOption("scipen")
	on.exit(options(scipen = baseline.scipen))
	options(scipen = 6)
	
	with(
		x$species.delta.table,
		{
			`exotic occupancy` = occupancy
			`delta J` = delta.J.Bars
			plot(				
				y = `delta J`,
				x = `exotic occupancy`,
				main = main,
				pch = pch,
				cex = cex,
				xaxs = xaxs,
				xlim = xlim,
				ylab = ylab,
				lwd = lwd,
				cex.axis = cex.axis,
				cex.lab = cex.lab,
				cex.main = cex.main,
				...
			)
		}
	)
	with(x, lines(scoop))
	abline(h = 0, lty = 2, lwd = 1.5)
	with(
		x, 
		{
			abline(v = p.Star, lty = 3, lwd = 1.5)
			abline(v = p.Star/2, lty = 3, lwd = 1.5)
			
		}
	)
}