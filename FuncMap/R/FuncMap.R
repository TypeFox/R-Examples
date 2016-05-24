FuncMap <-
function(fwb = foodweb, pkg = "none",
	method = "abs",
	sm.title = FALSE, newpage = TRUE, debug = FALSE) {
	
	# argument foodweb must be of class foodweb
		
#########################################################

#	a couple of convenience functions before we begin:

	p2cX <- function(r, theta) x <- r*cos(theta*2*pi/360)
	p2cY <- function(r, theta) x <- r*sin(theta*2*pi/360)

#########################################################

	fm <- fwb$funmat
	
	# Figure out which functions go in which categories
	
	inc <- colSums(fm) # total incoming function calls
	out <- rowSums(fm) # total outgoing calls
	# note that inc and out are named vectors, very handy!
	
	tot <- inc + out # total action for each func
	for (n in 1:length(tot)) if (tot[n] == 0) tot[n] <- NA # for standalones
	mid <- inc/out

	# properly interpreted, use mid to assign function to an axis:
	# mid = 0 is a top level function (source)
	# mid = Inf is a bottom level function (sink)
	# mid = 0/0 gives NaN, these are functions that neither
	# call nor get called (remove them eventually)
	# any other value is both a source and sink (mid)

	w.nan <- which(is.nan(mid)) # these are the stand alones
	nf <- nrow(fm) # no. of functions
	nsa <- length(w.nan) # no. of stand alones
	
	# Express in polar coordinates for plotting the axes
	# These are absolute units (native in grid)
	# unless method = rank or norm
	
	axis <- rep("mid", length(inc)) # all axes are mid to begin
	theta <- rep(120, length(inc)) # initialize angle
	r <- rep(NA, length(inc)) # initialize radius
	nomid <- FALSE # logical for when there is no mid axis
	
	w <- which(mid == Inf) # set up sink axis
	theta[w] <- 240
	axis[w] <- "sink"
	r[w] <- inc[w] # set radius for sink axis
	
	if (method == "rank") r[w] <- rank(inc[w], ties.method = "first")
	if (method == "norm") {
		if (length(unique(r[w])) == 1) stop("Can't do a norm plot, only 1 value on sink axis")
		r[w] <- (r[w]-min(r[w]))/max(r[w]-min(r[w]))		}
	
	w <- which(mid == 0) # set up source axis
	theta[w] <- 0
	axis[w] <- "source"
	r[w] <- out[w] # set radius for source axis
	
	if (method == "rank") r[w] <- rank(out[w], ties.method = "first")
	if (method == "norm") {
		if (length(unique(r[w])) == 1) stop("Can't do a norm plot, only 1 value on source axis")
		r[w] <- (r[w]-min(r[w]))/max(r[w]-min(r[w]))
		}

	w <- which(axis == "mid") # these are left overs: not set to sink or source
	r[w] <- tot[w] # set radius for mid axis
	if (all(is.na(r[w]))) nomid <- TRUE # these pkgs have no managers
	
	if (method == "rank" & !(nomid)) {
		# no matter what, the length can't be truncated
		r[w] <- rank(tot[w], ties.method = "first", na.last = "keep")
		}

	if (method == "norm" & (!nomid)) {
		if (length(unique(na.omit(r[w])))== 1) stop("Can't do a norm plot, only 1 value on mid axis")
		mn <- min(r[w], na.rm = TRUE)
		r[w] <- (r[w]-mn)/max(r[w]-mn, na.rm = TRUE)
		}
		
	# Put into data frame
		
	pts <- data.frame(axis = axis, theta = theta, radius = r, 
		row.names = names(inc))
	
	# Set up center hole which improves the display
	
	ch <- floor(min(na.omit(pts$radius)))
	if (method == "norm") ch <- 0.1; pts$radius <- pts$radius + ch
	if (method == "abs") pts$radius <- pts$radius - ch
	if (method == "rank") pts$radius <- pts$radius - ch
	
	# Convert to Cartesian coord as these are needed for the plot functions
	
	pts$x <- p2cX(pts$r, pts$theta)
	pts$y <- p2cY(pts$r, pts$theta)
	
	pts2 <- na.omit(pts) # Cleaned up copy for various tasks
	
	# Set up grid graphics viewport
	
	md <- max(abs(range(pts2[,c("x", "y")])))*1.3 # max dimension
	vp <- viewport(x = 0.5, y = 0.5, width = 0.8, height = 0.8,
		xscale = c(-md, md), yscale = c(-md, md), name = "HivePlotTest")
	if (newpage) grid.newpage()
	pushViewport(vp)
	
	# Add title & scale information.  sm.title for comparison plots
		
	if (sm.title) grid.text(pkg, y = 1.0, gp = gpar(fontsize = 14))
	if (sm.title) {
		if (method == "abs") grid.text("call count", y = 0.93, gp = gpar(fontsize = 8))
		if (method == "norm") grid.text("normalized call count", y = 0.93, gp = gpar(fontsize = 8))
		if (method == "rank") grid.text("ranked call count", y = 0.93, gp = gpar(fontsize = 8))
		}
	
	if (!sm.title) {
		grid.text(paste("Hive Plot Function Map of", pkg, "Package", sep = " "),
			y = 1.1, gp = gpar(fontsize = 14))	
		grid.text(paste(nf, "functions total;", nsa, "are stand alone", sep = " "),
			y = 1.05, gp = gpar(fontsize = 10))
		
		if (method == "norm") grid.text("position along axis is number of calls normalized to 0:1 interval",
			y = 0, gp = gpar(fontsize = 10))
			
		if (method == "abs") grid.text("position along axis is count of total calls",
			y = 0, gp = gpar(fontsize = 10))
			
		if (method == "rank") grid.text("position along axis is total calls by rank",
			y = 0, gp = gpar(fontsize = 10))
		}
	
	# Add axes using pts2 & label them
	# Two similar code sections for nomid = T/F
	
	if (!nomid) {
			
		x0 <- p2cX(r = c(ch, ch, ch), theta = c(0, 120, 240)) # this creates
		y0 <- p2cY(r = c(ch, ch, ch), theta = c(0, 120, 240)) # the center gap
	
		r.source <- subset(pts2, axis == "source")
		rmax.source <- max(r.source$radius)
		r.mid <- subset(pts2, axis == "mid")
		rmax.mid <- max(r.mid$radius)		
		r.sink <- subset(pts2, axis == "sink")
		rmax.sink <- max(r.sink$radius)	
		x1 <- p2cX(c(rmax.source, rmax.mid, rmax.sink), c(0, 120, 240))
		y1 <- p2cY(c(rmax.source, rmax.mid, rmax.sink), c(0, 120, 240))
		r.max <- max(rmax.source, rmax.mid, rmax.sink)
	
		grid.segments(x0 = x0, y0 = y0, x1 = x1, y1 = y1,
			gp = gpar(col = c("green", "blue", "red"), lwd = 3),
			default.units = "native")

		if (!sm.title) { # label offset along axis/radius
			if (method == "norm") fix1 <- fix2 <- fix3 <- 1.3
			if (method == "rank" | method == "abs") {
				fix1 <- ceiling(rmax.source * 1.1)
				fix2 <- ceiling(rmax.sink * 1.1)
				fix3 <- ceiling(rmax.mid * 1.1)
				# really small packages need a kluge:
				if (r.max < 5) {fix1 <- fix1 - 0.5; fix2 <- fix2 - 0.5; fix3 <- fix3 - 0.5}
				}

			grid.text("source",
				x = p2cX(fix1, 0), y = p2cY(fix1, 0), just = "left",
				default.units = "native", gp = gpar(fontsize = 10))
			grid.text("sink",
				x = p2cX(fix2, 240), y = p2cY(fix2, 240),
				default.units = "native", gp = gpar(fontsize = 10))
			grid.text("manager",
				x = p2cX(fix3, 120), y = p2cY(fix3, 120),
				default.units = "native", gp = gpar(fontsize = 10))
			}
		}		

	if (nomid) {
			
		x0 <- p2cX(r = c(ch, ch), theta = c(0, 240)) # this creates
		y0 <- p2cY(r = c(ch, ch), theta = c(0, 240)) # the center gap
	
		r.source <- subset(pts2, axis == "source")
		rmax.source <- max(r.source$radius)
		r.sink <- subset(pts2, axis == "sink")
		rmax.sink <- max(r.sink$radius)	
		x1 <- p2cX(c(rmax.source, rmax.sink), c(0, 240))
		y1 <- p2cY(c(rmax.source, rmax.sink), c(0, 240))
		r.max <- max(rmax.source, rmax.sink)
		grid.segments(x0 = x0, y0 = y0, x1 = x1, y1 = y1,
			gp = gpar(col = c("green", "red"), lwd = 3),
			default.units = "native")

		if (!sm.title) { # label offset along axis/radius
			if (method == "norm") fix1 <- fix2 <- 1.3
			if (method == "rank" | method == "abs") {
				fix1 <- ceiling(rmax.source * 1.1)
				fix2 <- ceiling(rmax.sink * 1.1)
				# really small packages need a kluge:
				if (r.max < 5) {fix1 <- fix1 - 0.5; fix2 <- fix2 - 0.5}
				}

			grid.text("source",
				x = p2cX(fix1, 0), y = p2cY(fix1, 0), just = "left",
				default.units = "native", gp = gpar(fontsize = 10))
			grid.text("sink",
				x = p2cX(fix2, 240), y = p2cY(fix2, 240),
				default.units = "native", gp = gpar(fontsize = 10))
			}
		}		

	# Add debugging grid if requested
	
	if (debug) { 
		grid.rect(gp = gpar(lty = "dashed", col = "gray")) # reference/guide

		if (method == "abs") {		
				grid.circle(x = 0, y = 0, r = 1:r.max,
					gp = gpar(col = "grey"), default.units = "native")			}

		if (method == "norm") {
			rad = c(0.25, 0.5, 0.75, 1.0) + ch
			grid.circle(x = 0, y = 0, rad,
				gp = gpar(col = "grey"), default.units = "native")			}
			
		if (method == "rank") {
			grid.circle(x = 0, y = 0, 1:r.max,
				gp = gpar(col = "grey"), default.units = "native")			}
		}
	
	# Now, connect the functions with calls using splines to produce arcs
	
	# Step through the fm: where ever there is a 1, capture the row/col
	# to access the x,y coordinates, assigning colors & widths along the way
	# Collect it all in a df for one drawing call in a moment
	
	x.st <- y.st <- x.end <- y.end <- grp <- c()
	
	for (r in 1:nrow(fm)) {
		
		for (c in 1:ncol(fm)) {
			if (!is.na(fm[r,c]) && fm[r,c] == 1) {
				if (pts$y[c] == pts$y[r] & pts$x[c] == pts$x[r]) {
					next
					# These functions appear to connect to themselves
					# but it's an artifact of two functions being
					# assigned to the same axis and radius
					# because colSums and rowSums are the same
					}
				x.st <- c(x.st, pts$x[r])
				y.st <- c(y.st, pts$y[r])
				x.end <- c(x.end, pts$x[c])
				y.end <- c(y.end, pts$y[c])
				
				if (pts$axis[r] == "source") {
					if (pts$axis[c] == "sink") grp <- c(grp, 1)
					if (pts$axis[c] == "mid") grp <- c(grp, 2)
					} # grp defined below
					
				if (pts$axis[r] == "mid") grp <- c(grp, 3)
				}
			}	
		}

	# Put into df for later use

	cur <- data.frame(x.st = x.st, y.st = y.st, x.end = x.end, y.end = y.end,
		grp = grp)

	# now, count the occurences of unique rows so as to make arcs of varying lwd

	lwds <- rep(1, dim(cur)[1])
	cur <- aggregate(lwds, by = as.list(cur), FUN = sum)
	colnames(cur)[6] <- "lwd"

	# grid.curve does not accept a vector of curvature values, so separate
	# by grp for different colors and curvatures
	# This also allows one to skip certain groups when a pkg doesn't have them
	
	cur1 <- subset(cur, grp == 1) # from source to sink
	cur2 <- subset(cur, grp == 2) # from source to mid
	cur3 <- subset(cur, grp == 3) # from mid to sink

	if (!dim(cur1)[1] == 0) { # needed b/c some pkgs only call sink via managers
		grid.curve(cur1$x.st, cur1$y.st, cur1$x.end, cur1$y.end,
		default.units = "native",	ncp = 5, square = FALSE,
		gp = gpar(col = "green", lwd = cur1$lwd), curvature = -0.5)
		}
		
	if (!dim(cur2)[1] == 0) {
		grid.curve(cur2$x.st, cur2$y.st, cur2$x.end, cur2$y.end,
		default.units = "native",	ncp = 5, square = FALSE,
		gp = gpar(col = "green", lwd = cur2$lwd), curvature = 0.5)
		}
		
	if (!dim(cur3)[1] == 0) {
		grid.curve(cur3$x.st, cur3$y.st, cur3$x.end, cur3$y.end,
		default.units = "native",	ncp = 5, square = FALSE,
		gp = gpar(col = "blue", lwd = cur3$lwd), curvature = 0.5)
		}

	# Add points
	
	if (!method == "rank") grid.points(pts$x, pts$y, pch = 20, gp = gpar(cex = 0.5))
	
	# Add a center point representing the stand alone functions
	
	if (nsa > 0) grid.points(0, 0, pch = 20, gp = gpar(cex = 0.5))

	# clean up

	popViewport()
	
	ans <- list(points = pts, curves = cur)
	}

