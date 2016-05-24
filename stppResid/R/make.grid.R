make.grid <- function(stw, grid = c(10, 10))
{
	if(!is.stwin(stw))
		stop("stw must be an object of class stwin")
	if((grid[1] < 1) || (grid[2] < 1))
		stop("grid dimensions must be positive")
	xc <- stw$xcoord
	yc <- stw$ycoord
	tc <- stw$tcoord	
	x.node <- seq(xc[1], xc[2], length.out = grid[1] + 1)
	y.node <- seq(yc[1], yc[2], length.out = grid[2] + 1)
	xmin <- x.node[-length(x.node)]
	xmax <- x.node[-1]
	ymin <- y.node[-length(y.node)]
	ymax <- y.node[-1]
	gr1 <- expand.grid(ymin, xmin)
	gr2 <- expand.grid(ymax, xmax)			
	grid.full <- data.frame(cbind(gr1[,2], gr2[,2], gr1[,1], gr2[,1]))
	names(grid.full) <- c("xmin", "xmax", "ymin", "ymax")
	y <- list(grid.full = grid.full)
	class(y) <- "stgrid"
	return(y)
}