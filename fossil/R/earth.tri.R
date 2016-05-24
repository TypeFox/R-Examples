`earth.tri` <-
function(long1, lat1, long2, lat2, long3, lat3) {
  if (((long1-long2)^2-(lat1-lat2))==32400) warning('site 1 and 2 are exactly opposite one another, and so the triangle is undefined')
  if (((long1-long3)^2-(lat1-lat3))==32400) warning('site 1 and 3 are exactly opposite one another, and so the triangle is undefined')
  if (((long2-long3)^2-(lat2-lat3))==32400) warning('site 2 and 3 are exactly opposite one another, and so the triangle is undefined')
  R <- 40041.47/(2*pi)
  cx <- deg.dist(long1,lat1,long2,lat2)/R
  bx <- deg.dist(long1,lat1,long3,lat3)/R
  ax <- deg.dist(long2,lat2,long3,lat3)/R
  A <- acos((cos(ax)-cos(bx)*cos(cx))/(sin(bx)*sin(cx)))
  B <- acos((cos(bx)-cos(cx)*cos(ax))/(sin(cx)*sin(ax)))
  C <- acos((cos(cx)-cos(ax)*cos(bx))/(sin(ax)*sin(bx)))
  SA <- R^2*((A+B+C)-pi)
  return(SA)
}

