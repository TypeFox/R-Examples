plotPNO <- function(x, subset = NULL, thinning = NULL, xlab = NULL, tail_threshold = 0, wm = FALSE, legend.pos = "topleft"){
	
	# calculate weighted means:
	# -------------------------
	wmean <- pno.weighted.mean(x, subset = subset)
		
	# subset matrix
	# ---------------------
	if (!is.null(subset))
		x <- x[, c(1, which(names(x) %in% subset))]
		
	for (i in 2:(dim(x)[2])){
		nf <- sum(x[, i])
		x[, i] <- x[, i] / nf
	}
	
	# delete 'zero tails'
	# ---------------------
	if (dim(x)[2] > 2)
	    zeros <- which(apply(x[, -1], 1, sum) <= tail_threshold)
	else 
	    zeros <- which(x[, -1] <= tail_threshold)
	if (length(zeros) > 0)
			x <- x[-zeros, ]
			
	# thin matrix
	# --------------------
	if (!is.null(thinning)) 
		x <- x[seq(1, dim(x)[1], length.out = thinning), ]
	
	
	col <- rainbow(dim(x)[2] - 1)
	max_val  <-  max(x[, -1])
	plot(x[, 1], x[, 2], type = "l", col = col[1], 
	    ylim = c(0, 	max_val),
		main = "Predicted niche occupancy",
		xlab = xlab,
		ylab = ""
		)
		
	if (dim(x)[2] > 2)
	    for (i in 3:(dim(x)[2]))
		    lines(x[, 1], x[, i], col = col[i - 1])
		
	# plot legend:
	# --------------------
	if (!is.null(legend.pos)){
		if (is.list(legend.pos))
		    legend(x = legend.pos$x, y = legend.pos$y,  			    legend = colnames(x)[2:dim(x)[2]], fill = col)
		else {
			if (legend.pos == "topleft")
			    lxy <- legend(x = min(x[, 1]), y = max(x[, -1]),  			        legend = colnames(x)[2:dim(x)[2]], fill = col)
			
			if (legend.pos == "bottomleft"){
				lxy <- legend(x = min(x[, 1]), y = min(x[, -1]),  			        legend = colnames(x)[2:dim(x)[2]], 
				    fill = col, plot = FALSE)$rect
				xx <- min(x[, 1])
				yy <- min(x[, -1]) + lxy$h
			    legend(x = xx, y = yy, fill = col, 					legend = colnames(x)[2:dim(x)[2]])
			}
			if (legend.pos == "topright"){
				lxy <- legend(x = min(x[, 1]), y = min(x[, -1]),  			        legend = colnames(x)[2:dim(x)[2]], 
				    fill = col, plot = FALSE)$rect
				xx <- max(x[, 1]) - lxy$w
				yy <- max(x[, -1]) 
			    legend(x = xx, y = yy, fill = col, 					legend = colnames(x)[2:dim(x)[2]])
			}
			if (legend.pos == "bottomright"){
				lxy <- legend(x = min(x[, 1]), y = min(x[, -1]),  			        legend = colnames(x)[2:dim(x)[2]], 
				    fill = col, plot = FALSE)$rect
				xx <- max(x[, 1]) - lxy$w
				yy <- min(x[, -1]) + lxy$h
			    legend(x = xx, y = yy, fill = col, 					legend = colnames(x)[2:dim(x)[2]])
			}
		}
	}
	
		
	# plot weighted means:
	# --------------------
	if (wm){
		for (i in seq(along = wmean))
			lines(rep(wmean[i], 2), range(x[, -1]), col = col[i],				lty = 3, lwd = 3)
	}
}