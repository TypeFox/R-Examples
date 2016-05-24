plot.fluss <-
function(x, subs, dims, folder = getwd(), xlims = NULL, ...){
	# extract data and table
	xr <- x$flux.res
	xt <- x$flux.table
	# prepare plot size
	d <- dims*200*5/360
	# prepare xlims when xlims is NULL
	if(is.null(xlims)){
		xlims <- range(as.vector(sapply(xr, function(x) range(x$fl.dat$orig.dat$time))))
	}
	# check on subs
	if(is.null(subs)){
		# start graphics device
		dev.new(width=d[2], height=d[1])
		# prepare plotting
		par(mfrow = dims)
		for(i in c(1:length(xr))){
			zero.line <- switch(xr[[i]]$fluss$ghg, CH4 = 1870, N2O = 323, CO2 = 388.5)
			plot.flux(xr[[i]], zero.line = zero.line, main=rownames(xt)[i], xlims = xlims, cex.axis=1.2, cex.lab=1.4, ...)
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
			tmp.flux <- xr[spots==j]
			tmp.nms <- rownames(xt)[spots==j]
			# start pdf device
			pdf(file=paste(folder, j, ".pdf", sep=""), width=d[2], height=d[1])
			# prepare plotting
			par(mfrow = dims)		
			#plot
			for(i in c(1:length(tmp.nms))){
				zero.line <- switch(xr[[i]]$fluss$ghg, CH4 = 1870, N2O = 323, CO2 = 388.5)
				plot.flux(tmp.flux[[i]], zero.line = zero.line, main=tmp.nms[i], xlims = xlims, cex.axis=1.2, cex.lab=1.4, ...)
			}
		dev.off()
		}
	}	
}

