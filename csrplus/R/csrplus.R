#' Moving Window Approach to Quadrat Counts
#' 
#' Calculates quadrat counts to estimate the intensity of a spatial point process through the moving window approach proposed by Bailey and Gatrell (1995).  
#' Event counts are computed within a window of a set size over a fine lattice of points within the region of observation.
#' 
#' @usage mwin(xcoord, ycoord, boundaryx, boundaryy, gridx, gridy, windowsizel, windowsizew)
#' 
#' @param xcoord a vector containing the x-coordinates of the observed point process
#' @param ycoord a vector containing the y-coordinates of the observed point process
#' @param boundaryx a vector containing the x-coordinates of the boundary
#' @param boundaryy a vector containing the y-coordinates of the boundary
#' @param gridx integer for the number of column lattice points 
#' @param gridy integer for the number of row lattice points 
#' @param windowsizel integer for the length of the window
#' @param windowsizew integer for the width of the window
#' 
#' @return xgrid a vector of the x-coordinates of the lattice points within the boundary
#' @return ygrid a vector of the y-coordinates of the lattice points within the boundary
#' @return quadrat a vector containing the number of events within each window sampled along the lattice
#' 
#' @import sp
#' 
#' @details 
#' The function first constructs a rectangular space based on the maximum and minimum values of the boundary.  
#' It then places a lattice over the rectangle, with a number of points in the lattice defined by the user.
#' The "point.in.polygon" function determines whether or not a specific lattice point falls within the region of observation.
#' If so, the function counts the number of events in a window of a specified size centered on that point.  If not, the function moves on to the next lattice point.
#' The function returns the location of lattice points within the region of observation in addition to the quadrat count at each point.  
#' 
#' @source
#' Cressie, Noel. "Chapter 8: Spatial Point Patterns." Statistics for Spatial Data. Revised ed. New York: John Wiley and Sons, 1993. N. pag. Print.
#' 
#' Gatrell, Anthony C. "Chapter 3: Introductory Methods for Point Patterns." Interactive Spatial Data Analysis. By Trevor C. Bailey. N.p.: Routledge, 1995. N. pag. Print.
#' 
#' @examples
#' # To load data corresponding to the location of earthquakes in California:
#' data(quake)
#' 
#' # To load data corresponding to the boundary:
#' data(boundary)
#' 
#' # To compute quadrat counts with a 40 x 40 lattice and 1 x 1 unit window:
#' m <- mwin(quake[,3], quake[,2], boundary[,1], boundary[,2], 40, 40, 1, 1)
#' 
#' # To plot the results (with the shading corresponding to the quadrat count):
#' layout(matrix(c(1,2), nc=2), widths = c(4, 1))
#' palette(rev(heat.colors((max(as.numeric(m$quadrat))-min(as.numeric(m$quadrat))))))
#' plot(m$xgrid, m$ygrid, col=m$quadrat, pch=15, cex=.8, 
#'      xlab="X-Coordinates", ylab="Y Coordinates", main="Quadrat Count")
#' lines(boundary[,1], boundary[,2])
#' breaks <- seq(min(as.numeric(m$quadrat)), (max(as.numeric(m$quadrat))), by=1)
#' plot.new()
#' plot.window(xlim = c(0, 1),ylim = range(breaks),xaxs = "i", yaxs = "i")
#' rect(0, breaks[-length(breaks)],1, breaks[-1],
#' col = rev(heat.colors(length(breaks) - 1)))
#' axis(2)
#' 
#' @export

mwin <- function(xcoord, ycoord, boundaryx, boundaryy, gridx, gridy, windowsizel, windowsizew){
  xmin <- min(boundaryx)
  xmax <- max(boundaryx)
  ymin <- min(boundaryy)
  ymax <- max(boundaryy)
  # To determine one string of the horizontal lattice:
  horizontalgrid <- seq(xmin, xmax, by=(abs(xmax-xmin))/gridx)
  # To replicate the horizontal string over a rectangular area: 
  b <- mat.or.vec(0,1)
  for (i in 1:length(horizontalgrid))
  {
    a <- rep(horizontalgrid[i], gridy)
    b <- cbind(b, a)
  }
  # To determine one string of the vertical lattice:
  verticalgrid <- seq(ymin, ymax, by=(abs(ymax-ymin))/gridy)
  # To specify a complete grid over a rectangular area:
  grid <- mat.or.vec(0,2)
  for (i in 1:dim(b)[[2]])
  {
    column <- cbind(b[,i], verticalgrid)
    grid <- rbind(column, grid)
  }
  # The above commands placed the lattice of points was over the rectangle, using
  # the minimum and maximum values of the boundary as vertices.  However, as the
  # boundary is not guaranteed to be a rectangle, it is necessary to use the 
  # "point.in.polygon" function to isolate the points in the area of observation. 
  p <- point.in.polygon(grid[,1], grid[,2], boundaryx, boundaryy)
  a <- length(which(p == 1)) # "1" corresponds to a point in the boundary
  grid <- cbind(grid, p)
  # To count the number of events in the moving window:
  quadrat <- matrix(nrow=a, ncol=1) 
  j <- 1
  for (i in 1:dim(grid)[[1]])
  {
    if(grid[i,3]==1)
    {
      a <- which(xcoord < (grid[i,1] + (windowsizel/2)) & xcoord > 
                   (grid[i,1] - (windowsizel/2)) & ycoord < (grid[i,2] + (windowsizew/2))
                 & ycoord > (grid[i,2] - (windowsizew/2)))
      quadrat[j,1] <- length(a)
      j <- j+1
    }
  }
  b <- which(grid[,3] == 1)
  # To define the lattice points within the boundary:
  gridnew <- cbind(grid[b,1], grid[b,2])
  # To connect the lattice points with the moving window counts:
  gridnew2 <- cbind(gridnew, quadrat)
  results <- data.frame(xgrid = gridnew2[,1], ygrid = gridnew2[,2], quadrat =
                          quadrat)
  return (results)
}

#' Nearest Neighbor Test Statistic:  Pielou
#' 
#' Calculates the nearest neighbor test statistic for a spatial point process based on that proposed by Pielou (1959).  
#' The test statistic is based on sample point to event distances.  It is equal to sum from i=1 to n of [(pi*lambda*(x_i)^2)/(n)], where x_i is equal to the point-to-event distance, lambda is equal to the intensity of the point process, and n is equal to the number of sample point-to-event distances.  The test statistic follows a normal distribution with mean 1 and variance 1/n.
#' 
#' @param xcoord a vector containing the x-coordinates of the observed point process
#' @param ycoord a vector containing the y-coordinates of the observed point process
#' @param points integer for the number of events for which to calculate point-to-event distances
#' @param boundaryx a vector containing the x-coordinates of the boundary
#' @param boundaryy a vector containing the y-coordinates of the boundary
#' @param lambda integer for the estimated intensity of the point process (likely the area of the region of observation divided by the number of events)
#' 
#' @return pstat a integer for the calcualted test statistic
#' 
#' @import sp
#' 
#' @details 
#' 
#' The function begins by proposing a sample point within the region of observation.  It then calculates and stores the distance from this point to the nearest event.  
#' The process is repeated for the number of sample points specified by the user.
#' Based on the stored distances, and the input for lambda, the test statistic of sum from i=1 to n of [(pi*lambda*(x_i)^2)/(n)] is calculated.
#' 
#' The test statistic serves as an alternative to that proposed by Clark and Evans (1954).  Instead of comparing event-to-event distances, the function makes use of sample point-to-event distances.
#' Knowing that the test statistic follows a normal distribution of mean 1 and variance 1/n, it can be used in the calculation of a z-statistic to evaluate the hypothesis of complete spatial randomness.
#' 
#' In order for the selection of sample points to not bias the results, it is suggested that the function be run multiple times to obtain an average test statistic .
#' 
#' The test statistic should be used with caution as the function does not account for edge effects.
#' As points along the border will have larger nearest neighbor distances, the normal approximation of the test statistic will underestimate the mean distance.  
#' When looking at event-to-point distances, it is expected that distances will be larger in a clustered process than in a process that exhibits complete spatial randomness.  
#' While the event-to-point distance will be small if the random event happens to be located in a cluster, there is a high probability that the sample point will be located in a sparsely populated region and therefore have a large nearest neighbor distance.  
#' By underestimating the mean distance, there is consequently more evidence to reject the null hypothesis of complete spatial randomness.
#' 
#' While the test statistic can provide a good first assessment of the null hypothesis of complete spatial randomness, it should not be relied upon as a definitive measure.  
#' More accurate conclusions can likely be drawn by comparing the observed process to simulations of a random process generated over the specific region of observation.  
#' 
#' @source
#' Cressie, Noel. "Chapter 8: Spatial Point Patterns." Statistics for Spatial Data. Revised ed. New York: John Wiley and Sons, 1993. N. pag. Print.
#' 
#' Gatrell, Anthony C. "Chapter 3: Introductory Methods for Point Patterns." Interactive Spatial Data Analysis. By Trevor C. Bailey. N.p.: Routledge, 1995. N. pag. Print.

#' @examples
#' # To load data corresponding to the location of earthquakes in California:
#' data(quake)
#' 
#' # To load data corresponding to the boundary:
#' data(boundary)
#' 
#' # To compute the one hundred values of the test statistic:
#' p <- mat.or.vec(100,1)
#' for (i in 1:100) {
#'    p[i] <- pielou(quake[,3], quake[,2], 30, boundary[,1], boundary[,2], 7.177) }
#' # To compute the average test statistic:
#' pavg <- mean(p)
#' # To calculate a z-statistic to evalute the null hypothesis of complete spatial randomness:
#' z <- (pavg-1)/sqrt(1/30) 
#' 
#' @export

pielou <- function(xcoord, ycoord, points, boundaryx, boundaryy, lambda)
{
  min <- mat.or.vec(points, 1)
  for (t in 1:points)
  {
    count <- 0
    results <- mat.or.vec(1, length(xcoord))
    while(count < 1)
    {
      # To propose a random point:
      xpoint <- runif(1, min=min(xcoord), max=max(xcoord))
      ypoint <- runif(1, min=min(ycoord), max=max(ycoord))
      # To determine if the proposed point falls within the region of interest:
      p <- point.in.polygon(xpoint, ypoint, boundaryx, boundaryy)
      if(p!=0) # If p!0 the point is in the region of interest
      {
        xpoint <- xpoint
        ypoint <- ypoint
        count <- count + 1
      }
    }
    # To find the minimum point-to-event distance:
    for (i in 1:length(xcoord))
    {
      results[i] <- sqrt((xpoint - xcoord[i])^2 + (ypoint - ycoord[i])^2)
    }
    min[t] <- min(results)
  }
  # To calculate the test statistic:
  pstat <- ((pi*lambda/points)*sum(min^2))
  return(pstat)
}
