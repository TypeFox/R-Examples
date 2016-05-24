ternary.grid <-
function(grid.size)
{
	lx = 2 / sqrt(3)
	nx = trunc(lx * grid.size)
	ny = grid.size
	x = seq(0, lx, length = nx)
	y = seq(0, 1, length = ny)
	xy = as.matrix(expand.grid(x = x, y = y))
	xyi = as.matrix(expand.grid(x = 1:nx, y = 1:ny))
	abc = xy2abc(xy)
	i = rowSums(abc < 0) == FALSE
	abc = abc[i,]
	xy = xy[i,]	
	xyi = xyi[i,]
	res = data.frame(xy, xyi, abc)
	names(res) = c('x', 'y', 'xi', 'yi', 'a', 'b', 'c')
	res
}

