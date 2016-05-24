dist2line=function(object.ppp, line.ends)
################################################################################
# Calculates perpendicular distance of a point process contained within a strip
# to center line of the strip they are contained in.  It also computes the position
# on the line. This is the inverse of the offset.points function.
#
# Arguments:
#
#   object.ppp - point process for observations in a strip
#   line.ends  - ends of line x0,y0,x1,y1
#
# Value:  list with elements
#    distVals  - vector of perpendicular distances
#    projection- dataframe of projected locations on the line
#
# Devin Johnson
# March 4 2008
################################################################################
{
  x.bar = object.ppp$x
  y.bar = object.ppp$y
  x0 = as.double(line.ends[1])
  y0 = as.double(line.ends[2])
  x1 = as.double(line.ends[3])
  y1 = as.double(line.ends[4])
  d = sum((x0-x1)^2 + (y0-y1)^2)
  u.bar = ((x.bar-x0)*(x1-x0) + (y.bar-y0)*(y1-y0))/d
  p.x = x0+u.bar*(x1-x0)
  p.y = y0+u.bar*(y1-y0)
  distVals = sqrt((x.bar-p.x)^2 + (y.bar-p.y)^2)
  return(list(distance=distVals, projection=data.frame(x=p.x,y=p.y)))
}

