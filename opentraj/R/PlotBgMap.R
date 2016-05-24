PlotBgMap <-
function( traj, ... ) {

  # extract the projection from the sp object
  hySplitProj <- CRS(proj4string(traj))
  
  # apply the projection to the background map
  canada <- spTransform(opentraj::worldmap, hySplitProj)
  
  # plot the map
  plot(opentraj::worldmap, border="white", col="lightgrey", ... )
}
