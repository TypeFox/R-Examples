plot.fluxes <-
function(x, dims, ghg = "all", subs = NULL, folder = getwd(), xlims = NULL, ask = TRUE, ...){
	# ambient levels
	al <- list(CH4 = 1870, N2O = 323, CO2 = 388.5)
	# extract data and table
	xres <- x$flux.res
	xt <- x$flux.table
	# prepare plot sizes
	d <- dims*200*5/360
	# check on ghgs
	show <- !is.na(xres)
	if(ghg!="all"){
		tmp <- rep(FALSE, 3)
		names(tmp) <- c("CO2", "CH4", "N2O")
		sel <- match(ghg, c("CO2", "CH4", "N2O")[show])
		if(sum(is.na(sel))!=0){warning(cat("Check your ghg argument.", "\n", "Fluxes can only be plotted when they have been estimated."))}
		tmp[sel[!is.na(sel)]] <- TRUE
		show <- tmp
	}
	# check on subs
	if(is.null(subs)){
		# control plotting behavior
   		if (ask) {
        	# start graphics device
        	dev.new(width=d[2], height=d[1])
        	if(sum(show)>1){
        		devAskNewPage(TRUE)
        	}    	
    	}
    	if(show["CO2"]){
    		if(!ask) dev.new(width=d[2], height=d[1])
			# select data
			xr <- xres$CO2
			# check consistency
			if(length(xr)!=prod(dims)){warning("Number of subfigures and number of positions on plot did not match.")}
			# prepare xlims when xlims is NULL
			if(is.null(xlims)){
				xlims <- range(as.vector(sapply(xr, function(x) range(x$fl.dat$orig.dat$time))))
			}
			# prepare plotting
			par(mfrow = dims)
			for(i in c(1:length(xr))){
				zero.line <- al$CO2
				plot.flux(xr[[i]], zero.line = zero.line, main=rownames(xt)[i], xlims = xlims, cex.axis=1.2, cex.lab=1.4, ...)
			}
		}
    	if(show["CH4"]){
			if(!ask) dev.new(width=d[2], height=d[1])
			# select data
			xr <- xres$CH4
			# check consistency
			if(length(xr)!=prod(dims)){warning("Number of subfigures and number of positions on plot did not match.")}
			# prepare xlims when xlims is NULL
			if(is.null(xlims)){
				xlims <- range(as.vector(sapply(xr, function(x) range(x$fl.dat$orig.dat$time))))
			}
			# prepare plotting
			par(mfrow = dims)
			for(i in c(1:length(xr))){
				zero.line <- al$CH4
				plot.flux(xr[[i]], zero.line = zero.line, main=rownames(xt)[i], xlims = xlims, cex.axis=1.2, cex.lab=1.4, ...)
			}
		}
    	if(show["N2O"]){
			if(!ask) dev.new(width=d[2], height=d[1])
			# select data
			xr <- xres$N2O
			# check consistency
			if(length(xr)!=prod(dims)){warning("Number of subfigures and number of positions on plot did not match.")}
			# prepare xlims when xlims is NULL
			if(is.null(xlims)){
				xlims <- range(as.vector(sapply(xr, function(x) range(x$fl.dat$orig.dat$time))))
			}
			# prepare plotting
			par(mfrow = dims)
			for(i in c(1:length(xr))){
				zero.line <- al$N2O
				plot.flux(xr[[i]], zero.line = zero.line, main=rownames(xt)[i], xlims = xlims, cex.axis=1.2, cex.lab=1.4, ...)
			}
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
			if(show["CO2"]){
				# select data
				xr <- xres$CO2
				# subset fluxes and names
				tmp.flux <- xr[spots==j]
				# check consistency
				if(length(tmp.flux)!=prod(dims)){warning("Number of subfigures and number of positions on plot did not match.")}
				# prepare names
				tmp.nms <- rownames(xt)[spots==j]
				# start pdf device
				pdf(file=paste(folder, "/CO2.", j, ".pdf", sep=""), width=d[2], height=d[1])
				# prepare xlims when xlims is NULL
				if(is.null(xlims)){
					xlims <- range(as.vector(sapply(xr, function(x) range(x$fl.dat$orig.dat$time))))
				}
				# prepare plotting
				par(mfrow = dims)		
				#plot
				for(i in c(1:length(tmp.nms))){
					zero.line <- al$CO2
					plot.flux(tmp.flux[[i]], zero.line = zero.line, main=tmp.nms[i], xlims = xlims, cex.axis=1.2, cex.lab=1.4, ...)
				}
				dev.off()
			}
			if(show["CH4"]){
				# select data
				xr <- xres$CH4
				# subset fluxes and names
				tmp.flux <- xr[spots==j]
				# check consistency
				if(length(tmp.flux)!=prod(dims)){warning("Number of subfigures and number of positions on plot did not match.")}
				# prepare names
				tmp.nms <- rownames(xt)[spots==j]
				# start pdf device
				pdf(file=paste(folder, "/CH4.", j, ".pdf", sep=""), width=d[2], height=d[1])
				# prepare xlims when xlims is NULL
				if(is.null(xlims)){
					xlims <- range(as.vector(sapply(xr, function(x) range(x$fl.dat$orig.dat$time))))
				}
				# prepare plotting
				par(mfrow = dims)		
				#plot
				for(i in c(1:length(tmp.nms))){
					zero.line <- al$CH4
					plot.flux(tmp.flux[[i]], zero.line = zero.line, main=tmp.nms[i], xlims = xlims, cex.axis=1.2, cex.lab=1.4, ...)
				}
				dev.off()
			}
			if(show["N2O"]){
				# select data
				xr <- xres$N2O
				# subset fluxes and names
				tmp.flux <- xr[spots==j]
				# check consistency
				if(length(tmp.flux)!=prod(dims)){warning("Number of subfigures and number of positions on plot did not match.")}
				# prepare names
				tmp.nms <- rownames(xt)[spots==j]
				# start pdf device
				pdf(file=paste(folder, "/N2O.", j, ".pdf", sep=""), width=d[2], height=d[1])
				# prepare xlims when xlims is NULL
				if(is.null(xlims)){
					xlims <- range(as.vector(sapply(xr, function(x) range(x$fl.dat$orig.dat$time))))
				}
				# prepare plotting
				par(mfrow = dims)		
				#plot
				for(i in c(1:length(tmp.nms))){
					zero.line <- al$N2O
					plot.flux(tmp.flux[[i]], zero.line = zero.line, main=tmp.nms[i], xlims = xlims, cex.axis=1.2, cex.lab=1.4, ...)
				}
				dev.off()
			}
		}
	}	
}

