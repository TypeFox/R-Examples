plotSinglescan <-
function(data.obj, singlescan.obj, chr = NULL, traits = NULL, show.alpha.values = NULL, standardized = TRUE, show.marker.labels = FALSE, mark.covar = FALSE, mark.chr = TRUE, plot.type = "h", overlay = FALSE, trait.colors = NULL, show.rejected.markers = FALSE, show.selected.markers = FALSE, relative.spacing = FALSE, gap.grades = 10, min.gap = 1, max.gap = 10, min.gap.genome = 1000, max.gap.genome = 1000000, fixed.scale = FALSE, ymax = NULL, cex.axis = 0.5, cex.points = 1, cex.alpha = 0.5, cex.labels = 2, cex.legend = 0.7){
	
	if(singlescan.obj$alpha.thresh[1] == "no permutations"){
		show.alpha.values <- FALSE
		}
	
	if(is.null(show.alpha.values)){alpha.to.show = singlescan.obj$alpha; show.alpha.values = TRUE}
	if(show.alpha.values){alpha.to.show = singlescan.obj$alpha}
	if(!show.alpha.values){alpha.to.show <- NULL}
	if(is.numeric(show.alpha.values)){alpha.to.show <- show.alpha.values}
	
	all.alpha <- singlescan.obj$alpha
	
	if(show.rejected.markers && show.selected.markers){
		stop("show.rejected.markers and show.rejected.markers cannot both be TRUE.")
		}
	
	D1.results <- singlescan.obj$singlescan.results
	D1.num <- rownames(D1.results[[1]])
	D1.order <- order(as.numeric(D1.num))
	
	marker.names <- data.obj$marker.names
	marker.num <- data.obj$marker.num
	
	if(show.rejected.markers || show.selected.markers){
	ind.markers <- data.obj$geno.for.pairscan
	if(is.null(ind.markers)){
		stop("select.markers.for.pairscan() be run before showing rejected or selected markers for pairscan.")
		}
	ind.marker.num <- get.marker.num(data.obj, colnames(ind.markers))
	}

	

	num.geno <- dim(D1.results[[1]])[1]
	num.pheno <- length(D1.results)
	
	
	if(is.null(D1.results)){
		stop("singlescan() must be run before plotting the results")
		}



	#pull out positional information about the covariates
	covar.info <- get.covar(data.obj)
	# str(covar.info)
	
	covar.names <- covar.info$covar.names
	covar.chr <- get.marker.chr(data.obj, covar.names)
	covar.loc <- get.marker.location(data.obj, covar.names)
	covar.table <- covar.info$covar.table	
	
	#integrate the covariate names with the marker names
	marker.chr <- get.marker.chr(data.obj, marker.names)
	marker.loc <- get.marker.location(data.obj, marker.names)

	marker.table <- cbind(c(marker.chr, covar.chr), c(marker.loc, covar.loc), c(marker.names, covar.names))
	
	sorted.marker.table <- sortByThenBy(marker.table, col.type = c("n", "n"))
	#make sure phenotypic covariates are at the end
	pheno.covar.locale <- which(sorted.marker.table[,1] == 0)
	if(length(pheno.covar.locale) > 0){
		pheno.covar.table <- sorted.marker.table[pheno.covar.locale,]
		sorted.marker.table <- rbind(sorted.marker.table[-pheno.covar.locale,], pheno.covar.table)
		}

	all.names <- sorted.marker.table[,3]
	all.num  <- get.marker.num(data.obj, all.names)

	
	if(is.null(chr)){
		chr <- unique(data.obj$chromosome)
		if(!is.null(covar.names)){
			chr <- unique(c(chr, covar.chr))
			}
		}

	all.chr <- as.numeric(sorted.marker.table[,1])
	all.loc <- as.numeric(sorted.marker.table[,2])

	if(is.null(traits)){
		traits <- names(D1.results)
		}
		
	covars <- colnames(covar.table)
	if(mark.covar && is.null(covars)){
		stop("pheno2covar() or marker2covar() needs to be run before covariates can be marked.")
		}

	col.mat <- matrix("black", nrow = num.geno, ncol = num.pheno)

	if(!overlay){
		if(mark.covar){
			col.mat[which(rownames(D1.results[[1]]) %in% covars),] <- "red"
			}
		}else{
		if(is.null(trait.colors)){
			trait.colors <- c("black", "blue", "purple", "darkgreen")
			}
		if(length(trait.colors) < length(traits)){
		 	trait.colors <- rep(trait.colors, length(traits)/4)
		 	trait.colors <- trait.colors[1:length(traits)]
			}
		for(i in 1:length(traits)){
			col.mat[,i] <- trait.colors[i]
			}
		if(mark.covar){
			col.mat[which(rownames(D1.results[[1]]) %in% covars),] <- "red"
			}
		}

	min.chr.length <- 0.1
	u_chr <- unique(all.chr)
	chr.lengths <- rep(NA, length(u_chr))
	for(ch in 1:length(u_chr)){
		chr.locale <- which(all.chr == u_chr[ch])
		chr.location <- all.loc[chr.locale]
		chr.lengths[ch] <- max(chr.location) - min(chr.location)
		}
	chr.lengths <- chr.lengths/max(chr.lengths)	
	chr.lengths[which(chr.lengths < min.chr.length)] <- min.chr.length
	


	if(relative.spacing){
		gaps.genome <- segment.region(min.gap.genome, max.gap.genome, gap.grades)
		gaps.plotting <- segment.region(min.gap, max.gap, gap.grades)
		
		relative.locale <- rep(0, length(all.chr))
		start.ind <- 1
		for(ch in 1:length(u_chr)){
			chr.locale <- which(all.chr == chr[ch])
			chr.location <- all.loc[chr.locale]
				
			if(length(chr.location) > 1){
				#find the marker distances from the start of the chromosome (0)
				neighbor.dists <- apply(consec.pairs(c(0,chr.location)), 1, function(x) x[2]-x[1])				

				for(i in 1:length(chr.location)){
					if(ch == 1 && i == 1){					
						last.pos <- 0
						}else{
						last.pos <- 	relative.locale[start.ind-1]
						}
					gap.locale <- which(gaps.genome <= neighbor.dists[i])
					if(length(gap.locale) == 0){
						relative.locale[start.ind] <- last.pos + min.gap
						}else{
						relative.locale[start.ind] <- last.pos + gaps.plotting[max(gap.locale)]
						}
					start.ind = start.ind + 1
					}			
				}else{
					#if there's only one marker in the chromosome, add the maximum
					#gap between the last marker from the last chromosome and this
					#marker
					relative.locale[start.ind] <- relative.locale[start.ind-1] + max.gap
					start.ind = start.ind + 1
					}
				}#end looping through chromosomes
		}else{
		relative.locale <- 1:length(all.chr)
		}
			
	
	
	#pull out the results we are plotting and sort them 
	#so the markers that were used as covariates get 
	#put in the right place. We can use the marker numbers
	#to order them
	result.order <- order(as.numeric(rownames(D1.results[[1]])))
	ordered.results <- lapply(D1.results, function(x) x[result.order,])
	ordered.col <- col.mat[result.order,]
	ordered.chr <- all.chr[result.order]
	ordered.names <- all.names[result.order]
	ordered.num <- all.num[result.order]
	
	p.covar.names <- ordered.names[which(ordered.chr == 0)]
	
	results.rows <- which(ordered.chr %in% chr)
	results.el <- which(names(ordered.results) %in% traits)
	results.to.plot <- NULL
	for(r in results.el){
		if(standardized){
			results.to.plot <- cbind(results.to.plot, ordered.results[[r]][results.rows,"t.stat"])
			}else{
			results.to.plot <- cbind(results.to.plot, ordered.results[[r]][results.rows,"slope"])
			}
		}


	final.cols <- ordered.col[results.rows, results.el]
	colnames(results.to.plot) <- names(D1.results)[results.el]

	

		if(show.rejected.markers || show.selected.markers){
			if(show.selected.markers){
				ind.locale <- match(c(ind.marker.num, colnames(covar.table)), rownames(results.to.plot))
				}else{
				not.selected <- setdiff(rownames(results.to.plot), c(ind.marker.num, colnames(covar.table)))
				ind.locale <- match(not.selected, rownames(results.to.plot))
				}
			}

		
		if(is.numeric(alpha.to.show)){
				alpha.to.include <- which(all.alpha %in% sort(alpha.to.show))
				if(length(alpha.to.include) < length(show.alpha.values)){
					cant.find <- setdiff(show.alpha.values, all.alpha)
					warning("The following alpha values were not calculated by singlescan():\n")
					cat(cant.find, sep = "\n")
					}
				}else{
				alpha.to.include <- NULL
				}
				
			if(length(alpha.to.include) > 0){
				alpha.thresholds <- rep(NA, length(alpha.to.include))
				for(a in 1:length(alpha.to.include)){
					alpha.thresholds[a] <- singlescan.obj$alpha.thresh[[alpha.to.include[a]]]
					}
				}else{
				alpha.thresholds = NULL	
				}
		

	if(!overlay){
		layout.mat <- get.layout.mat(length(results.el), "upright")
		}else{
		layout.mat <- matrix(1, 1, 1)
		}

	layout(layout.mat)	
	for(p in 1:ncol(results.to.plot)){

		#figure out the axix label
		if(!overlay){
			if(standardized){
				y.label <- paste(colnames(results.to.plot)[p], "[|Eff|/se]", sep = " ")
				}else{
				y.label <- paste(colnames(results.to.plot)[p], "Eff", sep = " ")	
				}
			}else{
			if(standardized){
				y.label <- "[|Eff|/se]"
				}else{
				y.label <- "Eff"
				}				
			}

		pheno.res <- results.to.plot[,p]

		if(standardized){
			if(!overlay){
				all.vals <- c(pheno.res, unlist(alpha.thresholds), 0)
				if(fixed.scale){
					all.vals <- c(results.to.plot, unlist(alpha.thresholds), 0)	
					}
				}else{
				all.vals <- c(results.to.plot, unlist(alpha.thresholds), 0)
					}
			}else{
				if(!overlay){
					all.vals <- pheno.res
					}else{
					all.vals <- results.to.plot	
					}
			}

		#create the window
		if(p == 1 || !overlay){
			par(mar = c(3, 4, 3, 2) + 0.1)
			plot.new()
			if(is.null(ymax)){
				ymax <- max(all.vals[is.finite(all.vals)])
				}else{
				ymax <- ymax	
				}
			plot.window(xlim = c(0, max(relative.locale, na.rm = TRUE)), ylim = c(min(all.vals, na.rm = TRUE), ymax))
		
			#shade the chromosome regions
			# if(!show.marker.labels){
				markers.used.locale <- which(ordered.num %in% rownames(results.to.plot))
				chr.id <- ordered.chr[markers.used.locale]
				poly.relative.locale <- c(0, relative.locale)
				poly.chr.id <- c(min(chr.id[which(chr.id != 0)]), chr.id)
				par(xpd = TRUE)
				for(ch in 1:length(chr)){
					x.min <- min(poly.relative.locale[which(poly.chr.id == chr[ch])])
					x.max <- max(poly.relative.locale[which(poly.chr.id == chr[ch])])
					if(ch %% 2 == 1 && mark.chr){ #shade chromosome regions if mark.chr is TRUE
						polygon(x = c(x.min, x.min, x.max, x.max), y = c(min(all.vals, na.rm = TRUE), ymax, ymax, min(all.vals, na.rm = TRUE)), col = "lightgray", border = NA)
						}

					if(!show.marker.labels){ #plot the chromosome labels if we are not plotting marker labels
						if(chr[ch] == 0){
							covar.x <- seq(x.min, x.max, length(covar.names))
							covar.y <- rep(min(all.vals, na.rm = TRUE)-((max(all.vals, na.rm = TRUE)-min(all.vals))*0.05), length(covar.names))
							text(x = covar.x, y = covar.y, labels = p.covar.names, cex = cex.axis, adj = 1, font = 2, srt = 90)
							}else{
							text(x = mean(c(x.min, x.max)), y = min(all.vals, na.rm = TRUE)-((max(all.vals, na.rm = TRUE)-min(all.vals, na.rm = TRUE))*0.05), labels = chr[ch], cex = cex.axis, font = 2)
							}
						}
					}
				par(xpd = FALSE)
				# }
		
				abline(h = 0)
				
				if(show.marker.labels){
					if(standardized){
						if(mark.chr){y.val <- max(all.vals, na.rm = TRUE)*-0.08}else{y.val <- max(all.vals, na.rm = TRUE)*-0.05}
						}else{
						if(mark.chr){y.val <- min(all.vals, na.rm = TRUE) - (max(all.vals, na.rm = TRUE)*0.08)}else{y.val <- min(all.vals, na.rm = TRUE) - (max(all.vals, na.rm = TRUE)*0.05)}
						
						}
				    par(xpd = TRUE)
					text(x = relative.locale, y = y.val, labels = ordered.names, srt = 90, cex = cex.axis, adj = 1)
					par(xpd = FALSE)
					}

			if(show.selected.markers){
				ind.position <- relative.locale[ind.locale]
				}
			if(show.rejected.markers){
				ind.position <- relative.locale[ind.locale]
				}

			if(standardized){

				#add the lines for the significance thresholds
				if(length(alpha.thresholds) > 0){
					for(a in 1:length(alpha.thresholds)){				
						segments(x0 = min(poly.relative.locale), y0 = alpha.thresholds[a], x1 = max(poly.relative.locale), lty = a, lwd = cex.points)
						}
	
					#add labels for significance lines
					par(xpd = TRUE)
					for(a in 1:length(singlescan.obj$alpha)){
						text(x = max(relative.locale)*1.02, y = alpha.thresholds[a], labels = paste("p =", singlescan.obj$alpha[a]), cex = cex.alpha, adj = 0)
						}
					}
									
				par(xpd = FALSE)
				
								
				if(!overlay){
					mtext(colnames(results.to.plot)[p], cex = cex.labels)
					mtext(y.label, side = 2, line = 2.5)
					}else{
					mtext(y.label, side = 2, line = 2.5)	
					}
		
				axis(2, cex.axis = cex.axis)
				abline(h = 0)
				
				}
			}
			
			if(show.selected.markers || show.rejected.markers){

				if(standardized){
					par(xpd = TRUE)
					points(ind.position, rep((0-max(pheno.res, na.rm = TRUE)*0.02), length(ind.locale)), col = "red", pch = "*", cex = cex.points)
					par(xpd = FALSE)
					}else{
					points(ind.position, (pheno.res[ind.locale]+(max(pheno.res, na.rm = TRUE)*0.02)), col = "red", pch = "*", cex = cex.points)		
					}
				}

		
			if(p == 1){
				par(xpd = TRUE)
				if(mark.covar){
					if(plot.type == "p" || plot.type == "b"){
						legend(x = (0-(max(relative.locale)*0.04)), y = max(all.vals, na.rm = TRUE)*1.15, legend = "covariate", pch = 16, col = "red", cex = cex.legend)
						}
					if(plot.type == "h"){
						legend(x = (0-(max(relative.locale, na.rm = TRUE)*0.04)), y = max(all.vals, na.rm = TRUE)*1.15, legend = "covariate", lty = 1, col = "red", cex = cex.legend)
						}
					}
				if(show.selected.markers){
					legend(x = (max(relative.locale, na.rm = TRUE)*0.94), y = max(all.vals, na.rm = TRUE)*1.15, legend = "selected", pch = "*", col = "red", cex = cex.legend)
					}
				if(show.rejected.markers){
					legend(x = (max(relative.locale, na.rm = TRUE)*0.94), y = max(all.vals, na.rm = TRUE)*1.15, legend = "rejected", pch = "*", col = "red", cex = cex.legend)
					}
				par(xpd = FALSE)
				}
			
			
		marker.chr <- ordered.chr[which(names(pheno.res) %in% ordered.num)]
		if(standardized){
			#plot the effect sizes
			for(ch in 1:length(chr)){
				chr.locale <- which(marker.chr == chr[ch])
				if(ch != 0){
					points(x = relative.locale[chr.locale], y = pheno.res[chr.locale], type = plot.type, col = ordered.col[chr.locale,p], pch = 16, lwd = cex.points, cex = cex.points)
					}else{
					points(x = relative.locale[chr.locale], y = pheno.res[chr.locale], type = "h", col = ordered.col[chr.locale,p], pch = 16, cex = cex.points)	
					}
				}
			
			}else{
			
			for(ch in 1:length(chr)){
				chr.locale <- which(marker.chr == chr[ch])
				if(ch != 0){
					points(x = relative.locale[chr.locale], y = pheno.res[chr.locale], type = plot.type, col = ordered.col[chr.locale,p], pch = 16, cex = cex.points, lwd = cex.points)
					}else{
					points(x = relative.locale[chr.locale], y = pheno.res[chr.locale], type = "h", col = ordered.col[chr.locale,p], pch = 16, lwd = cex.points)	
					}
				}
			}	
				
			}
			
		
		if(overlay){
			par(xpd = TRUE)
			if(show.rejected.markers || show.selected.markers){
				legend(x = (max(relative.locale, na.rm = TRUE)*0.82), y = max(all.vals, na.rm = TRUE)*1.15, legend = traits, pch = 16, col = trait.colors[1:length(traits)], cex = cex.legend)
				}else{
				legend(x = (max(relative.locale, na.rm = TRUE)*0.94), y = max(all.vals, na.rm = TRUE)*1.15, legend = traits, pch = 16, col = trait.colors[1:length(traits)], cex = cex.legend)	
				}
			par(xpd = FALSE)
		}
		
			

}
