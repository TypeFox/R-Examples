ComputeAzimuth <-
function(lat1, lat2, lon1, lon2)
{   
  rads <- 180 / pi
  
  a_distance <- acos(sin(lat1 / rads) * sin(lat2 / rads) + cos(lat1 / rads) * cos(lat2 / rads) * cos((lon1 - lon2) / rads))
  
  acos((sin(lat2 / rads) - sin(lat1 / rads) * cos(a_distance)) / (cos(lat1 / rads) * sin(a_distance)))
}
