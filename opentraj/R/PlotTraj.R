PlotTraj <-
function( traj, ... ) { 
  # This function is designed to plot hySplit Forward and Backward trajectories calculated 
  # by the function ProcTraj
  # 
  # Args:
  #   traj: SpatialLines or SpatialLinesDataFrame or SpatialPoints data frame
  #   ...: Arguments passed to or from methods
  #
  # Results:
  #    A plot with a map in the background background
  
  #get the old par configuration
  oldpar <- par(no.readonly=TRUE)
  
  # it reduces the margin's size
  par(mar = c(0,0,0,0) + 2.0)
  
  # gets the bounding box of the object
  # this information will help to focus on the trajectory
  bb <- bbox(traj)
  
  PlotBgMap(traj, xlim=bb[1,], ylim=bb[2,], axes=TRUE)
  
  # print the grid
  grid(col="white")
  
  # plot the trajectory lines
  plot(traj, add=TRUE, ...)
  
  # get the coordinates from the trajectory
  cc <- coordinates(traj)

  box()
  
  # restore the par configuration
  par(oldpar)
}
