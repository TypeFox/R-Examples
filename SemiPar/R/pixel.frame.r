########## S-function: pixel.frame ##########

# For obtaining a grid of frame values from
# a pixel grid.

# Last changed: 04/16/98

pixel.frame <- function(pixel.grid)
{
   width.val <- pixel.grid[2] - pixel.grid[1]
   len <- length(pixel.grid)
   return(c(pixel.grid-0.5*width.val,pixel.grid[len]+0.5*width.val))

}

######### End of S-function pixel.frame ########
