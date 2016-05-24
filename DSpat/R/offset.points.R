offset.points=function(line,pts)
################################################################################
# Convert x,y point locations on the line and a distance (negative is left of
# line for the direction of travel) for a series of points in the strip.
#
# Arguments:
# line - vector of  (4 columns); line traverses from (x0,y0) to (x1,y1)
# pts  - dataframe of x,y,distance for each observed  point; x,y is the location
#        on the line that is perpendicular to the object; a negative distance
#        implies it is on the left side of the line as defined by the direction of
#        travel.
#
# Value: pts dataframe with x,y locations of the objects offset from the line at
#        the appropriate distance and side.
#
# Compute length of line, slope and set up x and y values for vertices
#
# 12/20/2007
# Jeff Laake
################################################################################
{
  d=sqrt((line[1]-line[2])^2+(line[3]-line[4])^2)
  slope=(line[4]-line[3])/(line[2]-line[1])
  if(!is.infinite(slope))
  {
     theta=pi/2-asin((line[4]-line[3])/d)
     deltay=abs(pts$distance)*sin(theta)
#   Vertices of strip must be listed so they are counter-clockwise
#   This depends on whether the line is going from left to right or right to left
       if(line[1]<line[2])
          y=pts$y-sign(pts$distance)*deltay
       else
          y=pts$y+sign(pts$distance)*deltay
       x=-slope*(y-pts$y)+pts$x
   }
#  Deal with special case of vertical line
   else
   {
       y=pts$y
       if(y[1]<y[2])
          x=pts$x+pts$distance
       else
          x=pts$x-pts$distance
   }
   pts$x=x
   pts$y=y
   return(pts)
}

