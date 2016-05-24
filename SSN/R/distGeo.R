distGeo <-
function(x.row, y.row, x.col, y.col, range = 1, minorp = 1, rotate = 0)
{
	if(range < 0)
		stop("range parameter less than 0 in distGeo")
	if(rotate < 0 || rotate > 180)
		stop("rotation parameter beyond 0 - 180 range in distGeo")
	if(minorp < 0 || minorp > 1)
		stop("minor range proportion beyond 0 - 1 range")
	# total number of row observations
	n.row <- length(x.row)
	# total number of column observations
	n.col <- length(x.col)
	# expand all x-coordinates for rows
	sxr <- outer(as.vector(x.row), rep(1, times = n.col))
	# expand all x-coordinates for columns
	sxc <- outer(rep(1, times = n.row), as.vector(x.col))
	# find difference in x-coordinates between all pairwise locations
	sxdif <- sxr - sxc
	# expand all x-coordinates for rows
	syr <- outer(as.vector(y.row), rep(1, times = n.col))
	# expand all x-coordinates for columns
	syc <- outer(rep(1, times = n.row), as.vector(y.col))
	# find difference in x-coordinates between all pairwise locations
	sydif <- syr - syc
	# rotate coordinates
	newx <- cos(rotate*.0174533)*sxdif - sin(rotate*.0174533)*sydif
	newy <- sin(rotate*.0174533)*sxdif + cos(rotate*.0174533)*sydif
	# scale coordinates by minor and major axes */
	newx <- newx/(range*minorp)
	newy <- newy/range
	# compute distance for the scaled and rotated coordinates */
	sqrt(newx^2 + newy^2)
}

