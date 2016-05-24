ComputeDistance <-
function(Lat1, Lat2, Lon1, Lon2)
{   
  rads <- 180 / pi
  gd <- 6371 * acos(sin(Lat1 / rads) * sin(Lat2 / rads) + cos(Lat1 / rads) * cos(Lat2 / rads) * cos((Lon1 - Lon2) / rads))
  return(gd)
}
