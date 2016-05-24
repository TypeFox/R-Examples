brownian.bridge <-
function(x, y, time.lag, location.error, area.grid=NULL, cell.size=NULL, time.step=10, max.lag=NULL){

    # Input: 
    #   x, y = Vectors of coordinates of ordered animal locations, in UTMs.
    #   time.lag = Vector of time differences (in minutes) between successive 
    #           locations. Length(time.lag) = length(x) - 1. 
    #   location.error = Either the standarad deviation of normally distributed 
	#			location error or a vector of SDs (1 for each location).
    #   area.grid = Matrix or data frame of "x" and "y" coords for brownian bridge. 
    #           If missing, defaults to minimum/maximum x/y minus/plus 1 SD of the  
    #           range of x/y. If area.grid is provided, cell sizes must be square and uniform. 
    #   cell.size = Cell size for grid, if grid not provided.
    #   time.step = The Brownian bridge probability density function 
	#			must be integrated to find the fraction of time spent in each region. While the 
	#			probability density function cannot be integrated, it can be approximated by 
	#			discretizing time into arbitrarily small intervals of time.step. The default 
	#			is 10 units (same as time.lag). Longer time.step speeds up estimation, 
	#			but reduces precision.
	#	max.lag = maximum lag between successful locations to consider in the analysis
    # Output: 
    #   UD = List with estimated Brownian.Motion.Variance, "x", "y", and "z", 
    #       where x and y are grid center point coordinates and z is the estimated
    #       probability of use with sum(z) = 1.0.
    
    if(is.null(x) | is.null(y) | (length(x) != length(y))) {
        stop("data is missing or unequal number of x and y coordinates")
    }
    if(is.null(location.error)) stop("must specify 'location.error'")
    if(is.null(area.grid) & is.null(cell.size)) {
        stop("'area.grid' or 'cell.size' must be specified")
    }
    if(!is.null(area.grid) & is.null(cell.size)){
        cell.size <- abs(area.grid[1,1] - area.grid[2,1]) 
    }
    if(is.null(area.grid) & !is.null(cell.size)){
        range.x <- range(x)
        range.y <- range(y)
        min.grid.x <- round(range.x[1] - 1*sd(x))
        max.grid.x <- round(range.x[2] + 1*sd(x))
        min.grid.y <- round(range.y[1] - 1*sd(y))
        max.grid.y <- round(range.y[2] + 1*sd(y))
        x. <- seq(min.grid.x, max.grid.x, cell.size) 
        y. <- seq(min.grid.y, max.grid.y, cell.size)
        area.grid <- merge(x., y.)
    }
	if(is.null(max.lag)){
		max.lag = max(time.lag)+1
	}    
    if(length(location.error) == 1){
        location.error <- rep(location.error, length(x))
    }
    
    n.locs <- length(x)

    BMvar <- brownian.motion.variance(n.locs, time.lag, location.error, x, y, max.lag)
    BMvar <- rep(BMvar, times=length(x))

    # Use 10 units (generally minutes) as default.
    if(is.null(time.step)) time.step <- 10

    grid.size <- nrow(area.grid)
    probability <- rep(0, grid.size)
    T.Total <- sum(time.lag)

    bbmm <- vector("list", 4)
    names(bbmm) <- c("Brownian motion variance", "x", "y", "probability")
    class(bbmm) <- "bbmm"
    
	probability <- NULL
	int <- 0
	for(i in 1:(n.locs-1)){
		if(time.lag[i] <= max.lag){
			theta <- NULL
			tm <- 0
			while(tm <= time.lag[i]){
				alpha <- tm/time.lag[i]
				mu.x <- x[i] + alpha*(x[i+1] - x[i])
				mu.y <- y[i] + alpha*(y[i+1] - y[i])
				sigma.2 <- time.lag[i]*alpha*(1-alpha)*BMvar[i] + 
					((1-alpha)^2)*(location.error[i]^2) + 
					(alpha^2)*(location.error[i+1]^2)
				ZTZ <- (area.grid[,1] - mu.x)^2 + (area.grid[,2] - mu.y)^2
				theta <- (1/(2*pi*sigma.2))*exp(-ZTZ/(2*sigma.2)) 
				int <- int + theta
				tm <- tm + time.step
			}
		}
	}
	#Scaling probabilities so they sum to 1.0
	probability <- int/T.Total
	probability <- probability/sum(probability)
    
    bbmm[[4]] <- probability   
    bbmm[[1]] <- BMvar[1]
    bbmm[[2]] <- area.grid[,1]
    bbmm[[3]] <- area.grid[,2]

    return(bbmm)

}
