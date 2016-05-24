plot.fluxxes <-
function(x, dims, subs = NULL, folder = getwd(), ask = TRUE, ...){
	# ambient levels
	al <- list(CH4 = 1870, N2O = 323, CO2 = 388.5)
	# extract data and table
	xres <- x$flux.res[[1]]
	xt <- x$flux.table
	# prepare plot sizes
	d <- dims*200*5/360
	# check on subs
	if(is.null(subs)){
		# control plotting behavior
   		if (ask) {
        	# start graphics device
        	dev.new(width=d[2], height=d[1])
        	devAskNewPage(TRUE)
    	}
    	if(!ask) dev.new(width=d[2], height=d[1])
		# check consistency
		if(length(xres)!=prod(dims)){warning("Number of subfigures and number of positions on plot did not match.")}
		# prepare plotting
		par(mfrow = dims)
		for(i in c(1:length(xres))){
			plot.fluxx(xres[i], ...)
		}
	}
	else{
		if(length(subs)==1){
			spots <- xt[,subs]
		}
		else{
			spots <- as.factor(as.character(apply(xt[,subs], 1, function(x) paste(x, collapse="."))))
		}	
		for(j in levels(spots)){
			# subset fluxes and names
			tmp.flux <- xres[spots==j]
			# check consistency
			if(length(tmp.flux)!=prod(dims)){warning("Number of subfigures and number of positions on plot did not match.")}
			# start pdf device
			pdf(file=paste(folder, "/", j, ".pdf", sep=""), width=d[2], height=d[1])
			# prepare plotting
			par(mfrow = dims)		
			#plot
			for(i in c(1:length(tmp.flux))){
				plot.fluxx(tmp.flux[i], ...)
			}
			dev.off()
		}
	}	
}

