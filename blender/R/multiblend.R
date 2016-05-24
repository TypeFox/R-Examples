# Functions for blending multiple landscapes & dealing with bundles

import = function(directory){
	# find the relevant CSV files and then put them in one big list
	files = list.files(
		directory, 
		pattern = " table\\.csv$", 
		full.names = TRUE
	)
	landscapes = structure(
		lapply(
			files, 
			function(x) read.csv(x, row.names = 1)
		),
		names = gsub("\\.csv$", "", files)
	)
}




sanitize = function(landscapes){
	# check a raw list of occupancy matrices for errors, then repackage
	# them for blend.landscape.list
	natives = landscapes[grepl("native", names(landscapes))]
	exotics = landscapes[grepl("exotic", names(landscapes))]
	
	native.names = gsub(" native.*$", "", names(natives))
	exotic.names = gsub(" exotic.*$", "", names(exotics))
	
	if(
		any(
			native.names[native.names %in% exotic.names] !=
				exotic.names[exotic.names %in% native.names]
		)
		){
		stop("blender is having trouble lining up your native and exotic landscapes.  See ?blend")	
	}
	
	if(length(native.names) == 0){
		stop("blender can't identify any complete landscapes.  See ?blend")
	}
	
	structure(
		list(
			native.landscapes = natives,
			exotic.landscapes = exotics,
			landscape.names = native.names[native.names %in% exotic.names]
		),
		class = "landscape.list"
	)
}




bundle = function(prebundle.list){
	# write up a summary data frame and append it to a list of ecoblender
	# results to make a bundle
	summary = as.data.frame(
		do.call(
			rbind, 
			lapply(prebundle.list, function(i) i$results.table)
		)
	)
	
	structure(
		c(prebundle.list, summary = list(summary)),
		class = "blended.landscape.bundle",
		names = c(row.names(summary), "summary")
	
	)
	
}





################
# bundle methods
################


print.blended.landscape.bundle = function(x, ...){
	print(x$summary)
}

plot.blended.landscape.bundle = function(x, ...){
	
	
	start.par = par("mfrow")
	on.exit(par(mfrow = start.par))
	par(mfrow = c(2, 2))
	
	subplot = function(x, y, ...){
		plot(y ~ x, ...)
		abline(0, 1)
		legend(
			"bottomright", 
			paste(
				"R2 =", 
				round(1 - var(x - y, na.rm = TRUE)/var(y, na.rm = TRUE), 3)
			), 
			bty = "n"
		)
	}
	
	with(
		x$summary,
		{
			subplot(
				x = J.Bar,
				y = J.Star,
				main = "Predicted vs. observed\nmean similarity",
				xlab = expression(bar(J)),
				ylab = "J*",
				...
			)
			
			
			subplot(
				x = delta.J.Bar,
				y = delta.J.Star,
				main = "Predicted vs. observed\nchange in mean similarity",
				xlab = expression(paste(Delta, bar(J))),
				ylab = expression(paste(Delta, " J*")),
				...
			)
			
			subplot(
				x = threshold,
				y = p.Star,
				main = "Predicted vs. observed\nhomogenization threshold",
				xlab = "threshold",
				ylab = "p*",
				...
			)
			
			subplot(
				x = nadir,
				y = p.Star / 2,
				main = "Predicted vs. observed\nsimilarity nadir",
				xlab = "nadir",
				ylab = "p*/2",
				...
			)
			
		}
	)
}