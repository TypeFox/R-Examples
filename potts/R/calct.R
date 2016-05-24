##############################################################################
#
# This file contains various functions used in calculating and storing
# the t(x) statistics for potts models.
#
# calc_t_full does the full image
#
# calc_t_innergrid does just an inner grid of the image
#
# calc_t does both by calling calc_t_full or calc_t_innergrid as appropriate
#
# everthing assumes colors are represented as 1...ncolor
# 
##############################################################################
calc_t_full <- function(x,ncolor) {
  # clever way to count memberships
  t_stat <- apply(matrix(1:ncolor, ncol=1), 1, function(s) sum(x==s))
  t_stat <- c(t_stat, 0)
  t_stat[ncolor+1] <-
    sum(apply(x, 1, function(r) sum(r==c(r[-1],r[1])))) +
      sum(apply(x, 2, function(r) sum(r==c(r[-1],r[1]))))
  names(t_stat) <- c(1:ncolor, "*")
  t_stat
}

calc_t_innergrid <- function(x, ncolor, grid, i, j) {
  ####################################################################
  # first copy a small grid to work on, augmented by 1 on each side
  dg <- dim(grid)
  dx <- dim(x)
  
  # copy the rows one at a time, cause I got tired of trying to be clever.
  horiz                <- (j-1):(j+dg[2])
  horiz[horiz == 0]    <- dx[2]
  horiz[horiz > dx[2]] <- horiz[horiz > dx[2]] - dx[2]
  
  vert               <- (i-1):(i+dg[1])
  vert[vert == 0]    <- dx[1]
  vert[vert > dx[1]] <- vert[vert > dx[1]] - dx[1]

  x_grid <- x[vert,horiz]
  
  x_grid[2:(dg[1]+1),2:(dg[2]+1)] <- grid
  # x_grid done!
  ####################################################################
  t_stat <- apply(matrix(1:ncolor, ncol=1), 1, function(s) sum(grid==s))
  t_stat <- c(t_stat, 0)
  t_stat[ncolor+1] <-
    # left one == grid
    sum(grid == x_grid[2:(dg[1]+1),1:dg[2]]) +
      # up one == grid
      sum(grid == x_grid[1:dg[1], 2:(dg[2]+1)]) +
        # right side == grid right edge
        sum(grid[,dg[2]] == x_grid[2:(dg[1]+1),dg[2]+2]) +
          # bottom side == grid bottom edge
          sum(grid[dg[1],] == x_grid[dg[1]+2,2:(dg[2]+1)])
  names(t_stat) <- c(1:ncolor, "*")
  t_stat
}

calc_t <- function(x, ncolor, grid=NULL, i=NULL, j=NULL) {
  if (is.null(grid)) {
    calc_t_full(x,ncolor)
  } else {
    calc_t_innergrid(x, ncolor, grid, i, j)
  }
}

