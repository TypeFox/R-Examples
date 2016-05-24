GetEdges <-
function(center, radius, nedges, orient = 45){
  # give little offset when lat = lon = 0
  if(center[1] == 0 | center[2] == 0) center = center + 1e-6

  # Earth constants
  earth.radius = 6378140.0
  radian = 180 / pi
  f = 1/298.257
  ecc = 0.08181922
  distance = radius ##still a bit messy to get an easy-to-use radius param.
  
  # Convert Lat / Lon to radians
  lat = center[2] / radian
  lon = center[1] / radian

  # Compute edges of desired fig, centered on Lat/Lon, and conforming to radius + orientation constrains
  res = NULL
  for(bearing in seq(orient, orient + 360, by = 360 / nedges)){
    b = bearing / radian
    R = earth.radius * (1 - ecc * ecc) / (1 - ecc * ecc * sin(lat)^2 )^1.5
#     $R = $EARTH_RADIUS_EQUATOR * (1 - $e * $e) / pow( (1 - $e*$e * pow(sin($lat),2)), 1.5);

    psi = distance / R
    phi = pi/2 - lat

    arccos = cos(psi) * cos(phi) + sin(psi) * sin(phi) * cos(b) 
    latA = (pi/2 - acos(arccos)) * radian
    arcsin = sin(b) * sin(psi) / sin(phi)
    lonA = (lon - asin(arcsin)) * radian
    res = rbind(res, c(lonA, latA))
    }
  res
  }

