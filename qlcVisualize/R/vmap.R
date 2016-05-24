# ====================
# make a dirichlet/voronoi tessellation of points in a window
# this is just a convenience wrapper around "dirichlet" from spatstat
# ====================

voronoi <- function(points, window) {

  p <- spatstat::ppp(points[,1],points[,2],window = window)
  v <- spatstat::dirichlet(p)

  if (!is.null(attr(p, "rejects"))) {
    rejected <- cbind(attr(p, "rejects")$x, attr(p, "rejects")$y)
    index <- apply(rejected, 1, function(x) {
                which(points[,1] == x[1] & points[,2] == x[2])
             })
    attr(v, "rejects") <- unlist(index)
  }

	return(v)
}

# ====================
# plotting of a voronoi-map (v-map)
# default plotting of tessellations in spatstat is not easy to use with colour filling
# ====================

vmap <- function(tessellation, col = NULL, add = FALSE, outer.border = "black", border = "grey", lwd = 1, ...) {

	if (!add) {
		plot(0,0
			, xlim = tessellation$window$xrange
			, ylim = tessellation$window$yrange
			, type = "n"
			, axes = FALSE
			, ann = FALSE
		)
	}

	tiles <- spatstat::tiles(tessellation)

	# repeat colors if necessary
	col <- rep(col, length.out = length(tiles))

	# plot all tiles individually, to allow for separate colors
	# vectorize the plotting using polygon()

	poly <- function(tile) {
	  parts <- sapply(tile$bdry, function(poly) {
	              coor <- cbind( x = poly$x, y = poly$y)
	              rbind(coor, c(NA, NA))
	               }, simplify = FALSE)
	  do.call(rbind, parts)
	}

	coor <- sapply(tiles, poly)
	coor <- do.call(rbind, coor)
	coor <- head(coor, -1)

	nr_polys <- sapply(tiles, function(tile){ length(tile$bdry) } )

	polygon(coor
	        , col = rep(col, times = nr_polys)
	        , border = border
	        , lwd = lwd
	        , ...
	        )

	# add outer border
	spatstat::plot.owin(tessellation$window
						, add = TRUE
						, border = outer.border
						, lwd = lwd
						)
}
