plotAncClim <- 
function(x, clades = NULL, col, density = TRUE, 
     tipmode = 1, nchar = 3, cex, tipspace, cladespace = 1, 
     lwd, ylab = ""){
	
	# get data:
	# ---------
	tr <- x$tree
	tips <- tr$tip.label
	nbtips <- length(tips)
	clim <- x$means
	
	if (density)
		cd <- x$central.density
		
	if (missing (cex)) cex <- 1
	if (missing (lwd)) lwd <- 1
		
	# if no clades are given use terminal sisters 
	# for coloration instead. Species in a grade
	# will receive their own color	
	if (is.null(clades)) {
		ts <- terminal.sisters(tr)
		nts <- tr$tip.label[!tr$tip.label %in% ts]
		for (i in seq(along = nts))
			ts <- rbind(ts, rep(nts[i], 2))
		lts <- vector(mode = "list", length = dim(ts)[1])
		for (i in seq(along = lts))
			lts[[i]] <- ts[i, ]
		clades <- lts
	}
	
	# abbreviate taxon labels:
	# -------------------------
	tips <- tr$tip
	abb <- function(tips)
		paste(head(unlist(strsplit(tips, "")), nchar), collapse = "")
	tips <- sapply(tips, abb)
		
	# calculate x-coordinate for nodes bases on node ages
	# and 'max_age'
	# -------------
	nodeages <- c(rep(0, nbtips), -branching.times(tr))
	max_age <- -max(branching.times(tr))
	
	# coordinates for plotting: node ages versus clima values
	# --------------------------------------------------------
	xy <- cbind(nodeages, clim)
	
	# plot coordinate system
	# ----------------------
	if (!density) yrange <- range(xy[, 2])			
	else yrange <- c(min(cd[1, ]), max(cd[2, ]))
	if (tipmode == 3) yrange[2] <- yrange[2] + 0.05 * diff(yrange)
	if (tipmode == 0) tipspace <- 0
	if (missing(tipspace)){
		tipspace <- 1 - (4 / nbtips) 	
	}
	tipspace <- -tipspace * max_age
	plot(x = c(max_age, tipspace), y = yrange,
		type = "n",
		xlab = "Time (Ma)", ylab = ylab,
		bty = "c", xaxp = c(floor(max_age), 0, 5))
	
	# calculate space needed for tip labels:
	# -------------------------------------
	tipxpos <- max(strwidth(tips, units = "f", cex = cex)) * 1:nbtips
	tipxpos <- tipxpos * cladespace
	
	# calculate coordinates for tiplabels
	# -----------------------------------
	tex <- cbind(tipxpos, xy[1:nbtips, 2])
	rownames(tex) <-  tr$tip.label
	for (i in 2:length(clades)){
		id <- which(rownames(tex) %in% clades[[i]])
		extraspace <- (-tex[1, 1]  + tex[1, 1] * 1) * (i - 1)
	    tex[id, 1] <- tex[id, 1] + extraspace
	}
	rownames(tex) <-  tips
	maxtip <- max(tex[, 1])
	tex[, 1] <-  tex[, 1] / maxtip * tipspace
	totalspace <- max(tex[, 1])
			
	# edge colors:
	# ------------	
	n <- noi(tr, clades)
	n <- lapply(n, descendants, tree = tr, internal = TRUE)
	lincol <- rep("grey", dim(xy)[1])
	if (missing(col))
		col <- rainbow(length(n))
	for (i in seq(along = n)){
		lincol[n[[i]]] <- col[i]
	}
	
	# plot edges:
	# ----------------
	colind <- vector() # used to order tip colors!
	for (i in seq(along = tr$edge[, 1])){
		ind <- tr$edge[i, ]
		lines(xy[ind, 1], xy[ind, 2], lwd = lwd, col = lincol[ind[2]])
		colind <- c(colind, ind[2])
	}
	
	# plot zero line:
	# ----------------
	if (min(xy[, 2]) < 0 & max(xy[, 2] > 0))
	lines(x = c(min(xy[, 1]), 0), y = rep(0, 2), lwd = 0.8 * lwd, 
	    lty = "12", col = "gray50")
	
	# density:
	# ----------------
	if (density & tipmode != 0){
		for (i in seq(along = cd[1, ]))
		lines(rep(tex[i, 1], 2), cd[, i], 				
		    col = lincol[1:nbtips][i], lwd = lwd, lty = "12")
	}
	
	# plot tiplabels:
	# ----------------
	if (tipmode != 0)
	    if (tipmode == 1)
	        text(tex[, 1], tex[, 2], rownames(tex), cex = cex, 		
	            adj = 0.5, col = lincol[1:nbtips], font = 4)
	    
	    if (tipmode == 2){
	    	points(tex, col = lincol[1:nbtips], pch = 18)
	    	id <- mean(yrange)
	    	id <- tex[, 2] < id
	    	offset <- strheight(".", cex = cex, font = 4)
	    	text(tex[id, 1], cd[2,id] + offset, rownames(tex)[id], 
	    	    cex = cex, adj = c(0, 0.5), font = 4, srt = 90, 
	    	    col = lincol[1:nbtips][id])
	        text(tex[!id, 1], cd[1, !id] - offset, rownames(tex)[!id], 
	            cex = cex, adj = c(1, 0.5), srt = 90, 
	            col = lincol[1:nbtips][!id], font = 4)
	    }
	    if (tipmode == 3){
	        points(tex, col = lincol[1:nbtips], pch = 18)
	        text(tex[, 1], max(cd[2,]), rownames(tex), cex = cex,
	            adj = c(0, 0.5), col = lincol[1:nbtips], font = 4, 
	            srt = 45)
	    }
}
