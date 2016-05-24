plotNetwork <-
function(data.obj, collapsed.net = TRUE, trait = NULL, phenotype.labels = NULL, main.lwd = 4, inter.lwd = 3, label.cex = 1.5, percent.bend = 15, chr.gap = 1, label.gap = 5, positive.col = "brown", negative.col = "blue"){

	positive.col <- get.col(positive.col)[2]
	negative.col <- get.col(negative.col)[2]

	# require(shape)
	circle.dens = 0.0005
	center.x = 1; center.y = 1; radius = 2

	# chr.rel.length <- c(98.5, 103.9, 82.7, 88.6, 90.2, 79, 89.1, 76.2, 75.1, 77.9, 88, 63.9, 67.3, 66.4, 59, 57.8, 61.3, 59.4, 56.9)
	chr.rel.length <- c(1)
	
	
	used.markers <- sort(unique(as.numeric(c(data.obj$var.to.var.influences[,1], data.obj$var.to.var.influences[,2]))))
	used.marker.names <- get.marker.name(data.obj, used.markers)
	used.marker.chr <- get.marker.chr(data.obj, used.markers)
	u_chr <- unique(used.marker.chr)
	all.pos <- as.numeric(get.marker.location(data.obj, used.markers))
	
	#for each chromosome distribute the marker locations 
	#between 1 and the max for easier plotting
	for(ch in 1:length(u_chr)){
		marker.chr.locale <- which(used.marker.chr == u_chr[ch])
		marker.pos <- as.numeric(all.pos[marker.chr.locale])
		rec.marker.pos <- marker.pos - (min(marker.pos)-1)
		norm.marker.pos <- rec.marker.pos/max(rec.marker.pos)
		all.pos[marker.chr.locale]  <- norm.marker.pos
		}
	
	
	all.pos[which(all.pos == 0)] <- 1
	chr <- sort(unique(as.numeric(used.marker.chr)))
	just.chr <- chr[which(chr != 0)]
	covar.info <- get.covar(data.obj)
	pheno.covar <- which(covar.info$covar.type == "p")
	
	if(length(pheno.covar) > 0){
		covar.names <- covar.info$covar.names[pheno.covar]
		num.true.chr = length(just.chr)
		names(just.chr) <- rep("chr", length(just.chr))
		chr.names <- c(covar.names, just.chr)
		names(covar.names) <- rep("covar", length(covar.names))
		chr <- c(covar.names, just.chr)
		}else{
		num.true.chr = length(chr)
		chr.names <- chr
		names(chr) <- rep("chr", length(chr))
		}
		
	
	if(length(chr.rel.length) != num.true.chr){
		# warning("The relative lengths of the chromosomes will not be plotted.")
		rel.length <- rep(1, num.true.chr)
		}else{
		rel.length <- chr.rel.length/max(chr.rel.length)
		}

	if(length(pheno.covar) > 0){
		rel.length <- c(rep(0.2, length(covar.names)), rel.length)
		}
	

	#============================================================================================
	# internal functions
	#============================================================================================
	get.circle <- function(radius, center.x = 1, center.y = 1, dens = circle.dens){
		t <- seq(0,2*pi, dens)	
		x <- radius*cos(t)+center.x
		y <- radius*sin(t)+center.y
		result <- list(x,y); names(result) <- c("x", "y")
		return(result)
		}
	
	
	get.block.coord <- function(radius.coord, start, pts.per.chr, block.rel.locale){
		coord.x <- radius.coord$x[start:(start+pts.per.chr[i]-1)] #get the relative x and y coordinates for the block
		coord.y <- radius.coord$y[start:(start+pts.per.chr[i]-1)]
		x.start <- round(length(coord.x)*as.numeric(block.rel.locale[2]))
		if(x.start == 0){x.start = 1}
		x.end <- round(length(coord.x)*as.numeric(block.rel.locale[3]))
		if(x.end > length(coord.x)){x.end = length(coord.x)}
		if(x.start == x.end){
			if(x.start == 1 && x.end < length(coord.x)){
				x.end = x.end + 1
				}
			if(x.end >= length(coord.x) && x.start > 1){
				x.start <- x.start - 1
				}
			}	
		x.coord <- coord.x[x.start:x.end]
		y.start <- round(length(coord.x)*as.numeric(block.rel.locale[2]))
		if(y.start == 0){y.start = 1}
		y.end <- round(length(coord.y)*as.numeric(block.rel.locale[3]))
		if(y.end == 0){y.end = 1}
		if(y.end > length(coord.y)){y.end = length(coord.y)}
		if(y.start == y.end){
			if(y.start == 1 && y.end < length(coord.y)){
				x.end = x.end + 1
				}
			if(y.end >= length(coord.y) && y.start > 1){
				y.start <- y.start - 1
				}
			}	
		y.coord <- coord.y[y.start:y.end]
		return(cbind(rep(names(chr.blocks.locale)[ch], length(x.coord)), x.coord, y.coord))
		}

	#assign a chromosome and relative position to each block
	get.chr.pos <- function(block){
		marker.locale <- which(used.markers %in% get.marker.num(data.obj, block))
		chr <- unique(used.marker.chr[marker.locale])
		if(length(chr) > 1){
			chr.char <- paste(chr, collapse = ", ")
			stop(paste("There is linkage between markers on chromosomes ", chr.char,". Please try a high r2.thresh.", sep = ""))
			}
		min.pos <- min(all.pos[marker.locale])
		max.pos <- max(all.pos[marker.locale])
		if(min.pos == max.pos){
			max.pos <- max.pos + 1
			}
		total.length <- max(all.pos[used.marker.chr == chr])
		if(chr == 0){chr <- block}
		return(c(chr, min.pos/total.length, max.pos/total.length))
		}
	
	get.block.effect <- function(block){
		effects <- rep(NA, length(var.to.pheno)); names(effects) <- names(var.to.pheno)
		block.char <- as.logical(is.na(suppressWarnings(as.numeric(block[1]))))
		if(block.char && !var.char){
			block <- get.marker.num(data.obj, block)
			}
		marker.locale <- lapply(var.to.pheno, function(x) match(block, x[,1]))
		for(i in 1:length(marker.locale)){
			effects[i] <- mean(as.numeric(var.to.pheno[[i]][marker.locale[[i]], "|t.stat|"]), na.rm = TRUE)
			}
		return(effects)	
		}
		
	#this function returns the points on a line between two points
	get.line <- function(x0, y0, x1, y1, dens = circle.dens){
		x.pts <- segment.region(x0, x1, 1/dens, alignment = "ends")
		slope <- (y1-y0)/(x1-x0)
		y.pts <- slope*(x.pts - x0) + y0
		result <- list(x.pts, y.pts); names(result) <- c("x", "y")
		return(result)
		}
	#============================================================================================
	# end internal functions
	#============================================================================================

	#get coordinates for the multiple concentric circles we will use
	chr.radius <- get.circle(radius)
	inner.bar.radius = get.circle(radius*0.98)
	inter.radius = get.circle(radius*0.96)
	bend.to.radius = get.circle(radius/2)

	#divide into chromosomes
	num.chr = length(chr)
	gap = round((length(chr.radius$x)*chr.gap)/100) #number of values to skip for gap between chromosomes


	label.gap <- round((length(chr.radius$x)*label.gap)/100)
	full.length <- length(chr.radius$x) - (gap*num.chr) - label.gap
	full.chr <- rel.length*(1/sum(rel.length))
	pts.per.chr <- floor(full.length*full.chr)
	
	#A block with one marker needs to be able to get at least one point
	#on the map
	min.block.val <- 1/pts.per.chr
	
	if(collapsed.net){
		adj.mat <- data.obj$collapsed.net
		}else{
		adj.mat <- data.obj$full.net	
		}

		if(is.null(adj.mat)){
			stop("get.network() must be run before plotting the collapsed network.")
			}
		
		if(collapsed.net){
			blocks <- data.obj$linkage.blocks.collapsed
			}else{
			blocks <- data.obj$linkage.blocks.full	
			}
			
		all.markers <- as.vector(unlist(blocks))
				

		chr.pos <- t(sapply(blocks, get.chr.pos))
		colnames(chr.pos) <- c("chromosome", "min.position", "max.position")
		#center the mark positions
		chr.pos[which(chr.pos[,1] == 0),2:3] <- 0.5
		
		
		#get the average effect size for the var to 
		#pheno effects for each block
		var.to.pheno <- data.obj$max.var.to.pheno.influence
		var.char <- as.logical(is.na(suppressWarnings(as.numeric(var.to.pheno[[1]][1,1]))))
	
		block.effects <- t(sapply(blocks, get.block.effect))
		rel.block.effects <- block.effects/max(block.effects, na.rm = TRUE)
		
		#and each phenotype
		if(is.null(trait)){
			pheno <- names(data.obj$max.var.to.pheno.influence)				
			}else{
			pheno <- trait
			trait.locale <- which(trait %in% names(data.obj$max.var.to.pheno.influence))
			if(length(trait.locale) < length(trait)){
				not.found <- which(!trait %in% names(data.obj$max.var.to.pheno.influence))
				message("I couldn't find the following traits:")
				cat(trait[not.found], sep = "\n")
				return()
				}
			}
			
		if(is.null(phenotype.labels)){
			pheno.names <- pheno
			}else{
			pheno.names <- phenotype.labels
			if(length(phenotype.labels) != length(pheno)){
				stop("I'm detecting the wrong number of phenotype labels.")
				}	
			}
			
		#get the circle coordinates for each trait
		trait.circ <- vector(mode = "list", length = length(pheno))
		names(trait.circ) <- pheno
		start.radius <- 1.05
		gap.radius <- 0.05
		for(i in 1:length(pheno)){
			effect.radius = get.circle(radius*start.radius)
			trait.circ[[i]] <- effect.radius
			start.radius <- start.radius + gap.radius
			}
		
		
		label.radius = get.circle((radius*start.radius)+0.1)
		
		#if we need to filter chr.pos and adj.mat to include
		#only the phenotypes we are including
		all.block.pheno <- c(rownames(chr.pos), pheno)
		adj.mat <- adj.mat[,colnames(adj.mat) %in% all.block.pheno, drop = FALSE]
	
		main.effect.mat <- adj.mat[,which(colnames(adj.mat) %in% pheno), drop = FALSE]
	
	chr.coord.table <- NULL #the coordinates of the chromosomes
	block.coord.table <- NULL #the coordinates of the blocks for plotting interaction polygons
	inner.bar.coord.table <- NULL #the coordinates of the blocks for plotting target bars
	
	plot.new()
	#give the right margin a bit more room to write covariate names
	plot.window(xlim = c(min(label.radius$x), max(label.radius$x)*1.25), ylim = c(min(label.radius$y), max(label.radius$y)))
	par(mar = c(2,2,2,0))
	plot.dim <- par("usr")
	
	#add a legend
	x.range <- plot.dim[2] - plot.dim[1]
	y.range <- plot.dim[4] - plot.dim[3]
	legend.x.coord <- plot.dim[2] - (x.range*0.25)
	legend.y.coord <- plot.dim[3] + (y.range*0.1)
	legend(x = legend.x.coord, y = legend.y.coord, legend = c("Enhancing", "Suppressing"), col = c(positive.col, negative.col), lty = 1, lwd = 3)
	
	#segment the y coordinates of the gap region
	#in the outermost circle to evenly space the
	#label sticks and gaps
	pheno.label.starts <- round(segment.region(1, label.gap, num.points = length(pheno), alignment = "center"))

	#put the labels between the outermost trait circle and the edge of the plot
	label.x <- max(trait.circ[[length(trait.circ)]]$x) + (plot.dim[2] - max(trait.circ[[length(trait.circ)]]$x))*0.25
	arrow.start.x <- max(trait.circ[[length(trait.circ)]]$x) + (plot.dim[2] - max(trait.circ[[length(trait.circ)]]$x))*0.2
	
	for(tr in length(trait.circ):1){

		#add light gray bars to help the eye track the phenotypes
		#they should be staggered so the label sticks don't overlap
		circ.x <- trait.circ[[tr]]$x[pheno.label.starts[tr]:length(trait.circ[[tr]]$x)]
		circ.y <- trait.circ[[tr]]$y[pheno.label.starts[tr]:length(trait.circ[[tr]]$y)]
		points(circ.x, circ.y, type = "l", col = "lightgray", lwd = main.lwd)

		#and add phenotype labels
		arrow.end.x <- trait.circ[[tr]]$x[pheno.label.starts[tr]]
		arrow.y <- trait.circ[[tr]]$y[pheno.label.starts[tr]]
		segments(x0 = arrow.start.x, x1 = arrow.end.x, y0 = arrow.y, lwd = main.lwd, col = "lightgray")
		text(label.x, arrow.y, labels = pheno.names[tr], adj = 0)

		}
	
	
	start = label.gap
	for(i in 1:length(chr)){
		chr.x.coord <- chr.radius$x[start:(start+pts.per.chr[i]-1)]
		chr.y.coord <- chr.radius$y[start:(start+pts.per.chr[i]-1)]
		points(chr.x.coord, chr.y.coord, type = "l", lwd = main.lwd)
		chr.coord.table <- rbind(chr.coord.table, cbind(rep(chr[i], length(chr.x.coord)), chr.x.coord, chr.y.coord))
		
		if(names(chr)[i] == "covar"){
			text.adj = 0
			}else{
			text.adj = 0.5
			}
			
		text(mean(label.radius$x[start:(start+pts.per.chr[i]-1)]), mean(label.radius$y[start:(start+pts.per.chr[i]-1)]), chr.names[i], adj = text.adj, cex = label.cex)
		
		# get the number of blocks on this chromosome
		chr.blocks.locale <- which(chr.pos[,1] == chr[i])
		if(length(chr.blocks.locale) > 0){
			for(ch in 1:length(chr.blocks.locale)){
				main.effects <- main.effect.mat[chr.blocks.locale[ch],,drop=FALSE]
				block.rel.locale <- chr.pos[chr.blocks.locale[ch],,drop=FALSE]
				
				for(ph in 1:length(pheno)){
						
					if(main.effects[ph] != 0){ #if there are significant effects of this block, add them to the circle
						if(main.effects[ph] < 0){trait.col = negative.col}else{trait.col = positive.col}
						
						block.coord <- get.block.coord(radius.coord = trait.circ[[ph]], start, pts.per.chr, block.rel.locale)
	
						if(!collapsed.net || names(chr)[i] == "covar" || dim(block.coord)[1] == 1){
							pt.locale <- which(trait.circ[[ph]]$x == block.coord[1,2])
							pt.x1 <- trait.circ[[ph]]$x[pt.locale-1]; pt.x2 <- trait.circ[[ph]]$x[pt.locale]
							pt.y1 <- trait.circ[[ph]]$y[pt.locale-1]; pt.y2 <- trait.circ[[ph]]$y[pt.locale]
							points(c(pt.x1, pt.x2), c(pt.y1, pt.y2), type = "l", lwd = main.lwd, col = trait.col)
							# points(as.numeric(block.coord[,2]), as.numeric(block.coord[,3]), type = "p", pch = 16, col = trait.col)
							}else{
							points(as.numeric(block.coord[,2]), as.numeric(block.coord[,3]), type = "l", lwd = main.lwd, col = trait.col)
							}
						}
					
					#collect positions of the blocks for polygons and inner target bars on slightly smaller circles
					block.coord.table <- rbind(block.coord.table, get.block.coord(radius.coord = inter.radius, start, pts.per.chr, block.rel.locale))
					inner.bar.coord.table <- rbind(inner.bar.coord.table, get.block.coord(radius.coord = inner.bar.radius, start, pts.per.chr, block.rel.locale))				
	
					} #end looping through phenotypes
				} #end looping through blocks
			} #end case for if there are blocks on this chromosome
		start = start + pts.per.chr[i] + gap
		}
		
		
		#add the interactions
		just.m <- adj.mat[,-which(colnames(adj.mat) %in% pheno), drop = FALSE]
		#adjust the values, so we can scale the line width by effect size
		# just.m <- just.m/max(abs(just.m))
		for(i in 1:dim(just.m)[1]){
			sig.locale <- which(just.m[i,] != 0)
			
			if(length(sig.locale) > 0){
				
				for(s in 1:length(sig.locale)){
					if(just.m[i,sig.locale[s]] > 0){edge.col <- positive.col}else{edge.col = negative.col}
				
					#find the start and stop positions
					start.block <- rownames(just.m)[i]
					start.block.coord <- block.coord.table[which(block.coord.table[,1] == start.block),,drop=FALSE]
				
					end.block <- colnames(just.m)[sig.locale[s]]
					end.block.coord <- block.coord.table[which(block.coord.table[,1] == end.block),,drop=FALSE]
	
	
					#draw a polygon to connect the start block and stop position
					start.inter.x <- as.numeric(start.block.coord[,2])
					start.inter.y <- as.numeric(start.block.coord[,3])

					end.inter.x <- as.numeric(end.block.coord[,2])
					end.inter.y <- as.numeric(end.block.coord[,3])

					#sort the corners of the polygon according to their distance from the center of the circle
					start.mat <- matrix(c(start.inter.x[1], start.inter.y[1], start.inter.x[length(start.inter.x)], start.inter.y[length(start.inter.y)]), ncol = 2, byrow = TRUE)
					end.mat <- matrix(c(end.inter.x[1], end.inter.y[1], end.inter.x[length(end.inter.x)], end.inter.y[length(end.inter.y)]), ncol = 2, byrow = TRUE)
					rownames(start.mat) <- c("dist.start.min", "dist.start.max") 
					rownames(end.mat) <- c("dist.end.min", "dist.end.max") 
					start.dist <- order(apply(start.mat, 1, function(x) dist(matrix(c(x, center.x, center.y), ncol = 2, byrow = TRUE))))
					end.dist <- order(apply(end.mat, 1, function(x) dist(matrix(c(x, center.x, center.y), ncol = 2, byrow = TRUE))))
					
					#add a bar to indicate the block at the source end of the interaction
					start.bar.coord <- inner.bar.coord.table[which(inner.bar.coord.table[,1] == start.block),,drop=FALSE]
					points(as.numeric(start.bar.coord[,2]), as.numeric(start.bar.coord[,3]), type = "l", lwd = main.lwd, col = "darkgray")
					
					#add a bar to indicate the block at the target end of the interaction
					end.bar.coord <- inner.bar.coord.table[which(inner.bar.coord.table[,1] == end.block),,drop=FALSE]
					points(as.numeric(end.bar.coord[,2]), as.numeric(end.bar.coord[,3]), type = "l", lwd = main.lwd, col = "darkgray")
					

					#find the midpoint of the line between source and target
					mid.point.x <- mean(c(mean(start.inter.x), mean(end.inter.x)))
					mid.point.y <- mean(c(mean(start.inter.y), mean(end.inter.y)))
					
					
					#get the points on the line between the midpoint and the center of the circle
					pts.to.center <- get.line(mid.point.x, mid.point.y, center.x, center.y, dens = circle.dens)
					#find the point that makes this line the correct percentage long
					shifted.x <- pts.to.center$x[round(((percent.bend/100))*length(pts.to.center$x))]
					shifted.y <- pts.to.center$y[round(((percent.bend/100))*length(pts.to.center$y))]
					if(length(shifted.x) == 0){shifted.x = mid.point.x}
					if(length(shifted.y) == 0){shifted.y = mid.point.y}

					
					#make a spline curve based on the start and stop of the line plus this shifted midpoint
					inter.curve <- xspline(c(mean(start.inter.x), shifted.x, mean(end.inter.x)), y = c(mean(start.inter.y), shifted.y, mean(end.inter.y)), shape = 1, draw = FALSE)
					points(inter.curve$x, inter.curve$y, col = edge.col, type = "l", lwd = inter.lwd)
					# points(inter.curve$x, inter.curve$y, col = edge.col, type = "l", lwd = sqrt(abs(just.m[i,sig.locale[s]])))

					#add an arrowhead pointed at the target
					arrow.rad <- atan2((mean(end.inter.y)-shifted.y), (mean(end.inter.x)-shifted.x))
					arrow.deg <- arrow.rad*(180/pi)
					Arrowhead(x0 = mean(end.inter.x), y0 = mean(end.inter.y), arr.col = edge.col, arr.adj = 1, lcol = edge.col, angle = arrow.deg, arr.lwd = inter.lwd)
										
				}
			} #end case for when there are significicant interactions in this row
		}#end looping through markers
		
}
