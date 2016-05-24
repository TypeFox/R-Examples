plotSinglescan.heat <-
function(data.obj, singlescan.obj, standardized = TRUE, show.marker.labels = FALSE, show.chr.boundaries = TRUE, label.chr = TRUE, threshold.above = NULL, color = "blue", light.dark = "light", scale.fun = NULL){
	
	
	D1.results <- singlescan.obj$singlescan.results
	marker.names <- data.obj$marker.names
	ind.markers <- singlescan.obj$geno.for.pairscan
		
	if(is.null(D1.results)){
		stop("singlescan() must be run before plotting the results")
		}

	chr <- unique(data.obj$chromosome)
	traits <- names(D1.results)
	chr.names <- chr
	chr.names[which(chr.names == 0)] <- "C"

	results.to.plot <- NULL
	for(r in 1:length(D1.results)){
		if(standardized){
			results.to.plot <- cbind(results.to.plot, D1.results[[r]][,"t.stat"])
			}else{
			results.to.plot <- cbind(results.to.plot, D1.results[[r]][,"slope"])
			}
		}
	colnames(results.to.plot) <- names(D1.results)

	if(show.marker.labels){
		show.chr.boundaries <- FALSE
		}

	#get coordinates of the chromosome boundaries
	if(show.chr.boundaries){
		chromosomes <- data.obj$chromosome
		u_chr <- unique(chromosomes[which(!is.na(chromosomes))])
		chr.boundaries <- apply(matrix(u_chr, ncol = 1), 1, function(x) max(which(chromosomes == x))) + 0.5
		chr.boundaries <- c(0, chr.boundaries)
		}

	
	if(!is.null(threshold.above)){
		#set all values above the threshold to the threshold value
		large.val.locale <- which(results.to.plot > threshold.above)
		results.to.plot[large.val.locale] <- threshold.above
		}
	
		if(!is.null(scale.fun)){
			cl <- call(scale.fun, results.to.plot)
			results.to.plot <- eval(cl)
			}
	
	
		col <- get.col(color, light.dark = light.dark)
		my.palette <- colorRampPalette(col)
		num.traits <- dim(results.to.plot)[2]
		num.markers <- dim(results.to.plot)[1] 
		image(x = 1:num.markers, y = 1:num.traits, z = results.to.plot, col = my.palette(50), axes = FALSE, xlab = "Chromosome", ylab = "")
		#place a black border around the plot so we can see the edges
		polygon(x = c(0.5, num.markers+0.5, num.markers+0.5, 0.5), y = c(0.5, 0.5, num.traits+0.5, num.traits+0.5))
		
		#add the trait labels
		par(xpd = TRUE)

		x <- 0-num.markers*0.05
		y <- seq(1, num.traits, 1)
		text(x = rep(x, length(y)), y = y, labels = traits)
		par(xpd = FALSE)

		#add the chromosome labels
		chr.cols <- rep(c("darkgray", "white"), ceiling(length(chr)/2))

		poly.perc = 0.3; label.size <- 0.9
		plot.dim <- par("usr")
		poly.height <- plot.dim[1]*poly.perc
		poly.max = 0.48; poly.min <- poly.max-poly.height; poly.mid <- mean(c(poly.max, poly.min))
	
	
		if(show.chr.boundaries){
			par(xpd = TRUE)
			for(i in 1:(length(chr.boundaries)-1)){
	
				polygon(x = c(chr.boundaries[i], chr.boundaries[i+1], chr.boundaries[i+1], chr.boundaries[i]), y = c(poly.min, poly.min, poly.max, poly.max), col = chr.cols[i])
				segments(x0 = chr.boundaries[i+1], x1 = chr.boundaries[i+1], y0 = 0.5, y1 = num.traits+0.5, col = "darkgray")
				if(label.chr){
					text(x = mean(c(chr.boundaries[i], chr.boundaries[i+1])), y = poly.mid, cex = label.size, labels = chr.names[i])
					}
				}
				
			par(xpd = FALSE)
			}
	

		par(xpd = TRUE)
		if(show.marker.labels){
			text(x = 1:num.markers, y = rep(poly.max, num.markers), labels = marker.names, adj = 1, srt = 90, cex = 0.5)
			}
		par(xpd = FALSE)

	

}
