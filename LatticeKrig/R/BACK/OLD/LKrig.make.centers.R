LKrig.make.centers <- function(grid.info, nlevel, 
    NC.buffer, distance.type) {
    #
    # actual number of grid points is determined by the spacing delta
    # delta is used so that centers are equally
    # spaced in both axes and NC is the maximum number of grid points
    # along the larger range. So for a rectangular region
    # along the longer side there will be NC grid but for the shorter dimension
    # there will be less than NC.
    # Finally note that if NC.buffer is greater than zero the number of gird points in
    # both dimensions will be increased by this buffer around the edges:
    # A total of  NC + 2* NC.buffer grid points along the longer dimension
    #
    delta.level1 <- grid.info$delta
    delta.save <- mx <- my <- rep(NA, nlevel)
    #
    # build up series of nlevel multi-resolution grids
    # and accumulate these in a master list  called grid
    grid.all.levels <- NULL
    # loop through multiresolution levels decreasing delta by factor of 2
    # and compute number of grid points.
    # build in hook for buffer regions to differ in x and y currently they are the same except for
    # cylinder geometry
    NC.buffer.x <- ifelse(distance.type == "cylinder", 0, NC.buffer)
    NC.buffer.y <- NC.buffer
    for (j in 1:nlevel) {
        delta <- delta.level1/(2^(j - 1))
        delta.save[j] <- delta
        # the width in the spatial coordinates for NC.buffer grid points at this level.
        buffer.width.x <- NC.buffer.x * delta
        buffer.width.y <- NC.buffer.y * delta
        
        
        # rectangular case
        grid.list <- list(x = seq(grid.info$xmin - buffer.width.x, 
            grid.info$xmax + buffer.width.x, delta), y = seq(grid.info$ymin - 
            buffer.width.y, grid.info$ymax + buffer.width.y, 
            delta))
        if (distance.type == "cylinder") {
            # for a grid on a clyinder assume that xmin and xmax are the same point
            # and drop the last grid point at xmax
            grid.list$x <- grid.list$x[-length(grid.list$x)]
        }
        # save results
        mx[j] <- length(grid.list$x)
        my[j] <- length(grid.list$y)
        grid.all.levels <- c(grid.all.levels, list(grid.list))
    }
    
    # end multiresolution level loop
    # create a useful index that indicates where each level starts in a
    # stacked vector of thebasis function parameters.
    offset <- as.integer(c(0, cumsum(mx * my)))
    return(list(delta = delta, mx = mx, my = my, m = sum(mx * 
        my), NC.buffer.x = NC.buffer.x, NC.buffer.y = NC.buffer.y, 
        m.domain = (mx - 2 * NC.buffer.x) * (my - 2 * NC.buffer.y), 
        offset = offset, grid = grid.all.levels, delta.save = delta.save))
}
