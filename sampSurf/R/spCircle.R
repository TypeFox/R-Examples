spCircle = function(radius,
                    spUnits = CRS(projargs=as.character(NA)),
                    centerPoint = c(x=0, y=0),   #centerPoint
                    nptsPerimeter = 100,
                    spID = paste('circle',.StemEnv$randomID(),sep=':'),
                    ...
                   )
{   
#---------------------------------------------------------------------------
#
#   This routine will create an "sp" circular SpatialPolygons object based on
#   the arguments...
#
#   Arguments...
#     radius = the radius of the circle
#     spUnits = the CRS units
#     centerPoint = the circle's center
#     nptsPerimeter = the number of points forming the perimeter of the polygon
#     spID = the spatial ID for the object
#
#   Returns...
#     a list with...
#       spCircle = the SpatialPolygons circular object
#       location =  the SpatialPoints center point
#
#   Please note that you may want to rename the components of, e.g., the
#   spCircle object to make more sense. An example where it is used for
#   dbh is in the standingTree constructor, where spDBH == spCircle...
#
#       names(spDBH@polygons) = 'pgsDBH'
#       names(spDBH@polygons$pgsDBH@Polygons) = 'pgDBH'
#
#
#Author...									Date: 24-Oct-2011
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   make sure the center is a named vector of length 2...
#
    if(any(is.na(match(c('x','y'),names(centerPoint)))))
      stop('Please use names x and y for circle centerPoint vector')
    if(length(centerPoint) != 2)
      stop('Please supply one set of (x,y) coordinates for the plot center location.')

  
#
#   some other checks...
#
    if(radius <= 0)
      stop('radius must be positive!')
    if(nptsPerimeter < 20) {
      warning('Using less than 20 points for the circle perimeter is not recommended--set to 20')
      nptsPerimeter = 20
    }

    area = pi*radius*radius

    ##location = centerPoint

#
#   left half of the circle, then right...
#
    circ =  seq(0, 2*pi, len=nptsPerimeter)

#
#   make the circle outline...
#
    circle = matrix(c(centerPoint['x'] + radius*cos(circ),
                    centerPoint['y'] + radius*sin(circ), rep(1,nptsPerimeter) ), nrow=nptsPerimeter)  
    
#    
#   any little difference between start & end pts with identical() can mess up the
#   the sp package Polygon routine, so set the end point exactly to start, then transform...
#
    circle = rbind(circle, circle[1,])
    dimnames(circle) = list(NULL,c('x','y','hc'))

#
#   and make a SpatialPolygons object...
#
    pgCircle = Polygon(circle[,-3])                               #sans hc
    pgsCircle = Polygons(list(circPlot = pgCircle), ID = spID)
    spCircle = SpatialPolygons(list(pgsCircle = pgsCircle))       #takes a list of Polygons objects
  
#
#   no id for center point, but it can be added to be the same as spID when
#   we make a container class for the center points elsewhere...
#
    loc = matrix(centerPoint, nrow=1)
    colnames(loc) = names(centerPoint)
    location = SpatialPoints(loc, proj4string = spUnits)
    

    return( list(spCircle=spCircle, location=location) )
}   #spCircle
