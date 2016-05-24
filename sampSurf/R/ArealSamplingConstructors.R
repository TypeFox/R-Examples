#---------------------------------------------------------------------------
#
#   This file holds the S4 definition for the constructor methods of the
#   ArealSampling class & subclasses...
#
#   The methods include...
#     1. a constructor for 'circularPlot'
#     2. a constructor for 'pointRelascope' (Jan 2011)
#     3. a constructor for 'perpendicularDistance' (Jan 2011)
#     4. a constructor for 'distanceLimit' (Mar 2011)
#     5. a constructor for 'angleGauge' (6-Dec-2011)
#
#   Note that the sp package should be loaded for the complete functionality. 
#
#Author...									Date: 20-Aug-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#   generic definition...
#
if(!isGeneric("circularPlot")) 
  setGeneric('circularPlot',  
             function(radius, ...) standardGeneric('circularPlot'),
             signature = c('radius')
            )

if(!isGeneric("pointRelascope")) 
  setGeneric('pointRelascope',  
             function(angleDegrees, ...) standardGeneric('pointRelascope'),
             signature = c('angleDegrees')
            )

if(!isGeneric("perpendicularDistance")) 
  setGeneric('perpendicularDistance',  
             function(kpds, ...) standardGeneric('perpendicularDistance'),
             signature = c('kpds')
            )

if(!isGeneric("distanceLimited")) 
  setGeneric('distanceLimited',  
             function(distanceLimit, ...) standardGeneric('distanceLimited'),
             signature = c('distanceLimit')
            )

if(!isGeneric("angleGauge")) 
  setGeneric('angleGauge',  
             function(baf, ...) standardGeneric('angleGauge'),
             signature = c('baf')
            )

if(!isGeneric("lineSegment")) 
  setGeneric('lineSegment',  
             function(length, orientation, ...) standardGeneric('lineSegment'),
             signature = c('length', 'orientation')
            )

          
#================================================================================
#  1. method for functions and class circularPlot...
#
setMethod('circularPlot',
          signature(radius = 'numeric'),
function(radius,
         units = 'metric',
         spUnits = CRS(projargs=as.character(NA)),
         centerPoint = c(x=0, y=0),   #centerPoint
         description = 'fixed area circular plot',
         nptsPerimeter = 100,
         #spID = unlist(strsplit(tempfile('cp:',''),'\\/'))[2],
         #spID = paste('cp',format(runif(1,0,10000),digits=8),sep=':'),
         spID = paste('cp',.StemEnv$randomID(),sep=':'),
         ...
        )
{
#------------------------------------------------------------------------------
#
#   make sure the center is a named vector of length 2...
#
    if(any(is.na(match(c('x','y'),names(centerPoint)))))
      stop('Please use names x and y for centerPoint vector')
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

    location = centerPoint

#
#   left half of the circle, then right...
#
    circ =  seq(0, 2*pi, len=nptsPerimeter)

#
#   make the plot outline...
#
    circPlot = matrix(c(centerPoint['x'] + radius*cos(circ),
                        centerPoint['y'] + radius*sin(circ), rep(1,nptsPerimeter) ), nrow=nptsPerimeter)  
    
#    
#   any little difference between start & end pts with identical() can mess up the
#   the sp package Polygon routine, so set the end point exactly to start, then transform...
#
    circPlot = rbind(circPlot, circPlot[1,])
    
#    circPlot = circPlot %*% trMat
    dimnames(circPlot) = list(NULL,c('x','y','hc'))
    

#
#   and make a SpatialPolygons object if sp is available...
#
    pgCircPlot = Polygon(circPlot[,-3])                             #sans hc
    pgsCircPlot = Polygons(list(circPlot=pgCircPlot), ID=spID)
    spCircPlot = SpatialPolygons(list(pgsCircPlot=pgsCircPlot))      #takes a list of Polygons objects
  
#
#   no id for center point, but it can be added to be the same as spID when
#   we make a container class for the center points elsewhere...
#
    loc = matrix(centerPoint, nrow=1)
    colnames(loc) = names(centerPoint)
    location = SpatialPoints(loc, proj4string = spUnits)

    cp = new('circularPlot', radius=radius, area=area, perimeter=spCircPlot,
             description=description, units=units,
             location = location, spID=spID, spUnits=spUnits )

    return(cp)
}   #circularPlot constructor
)   #setMethod
    




          
#================================================================================
#  2. constructor method for class pointRelascope...
#
setMethod('pointRelascope',
          signature(angleDegrees = 'numeric'),
function(angleDegrees,
         units = 'metric',
         description = 'point relascope method',
         ...
        )
{
#------------------------------------------------------------------------------
#
#   get the angle in radians and area factor...
#
    angleRadians = .StemEnv$deg2Rad(angleDegrees)
    phi = (pi - angleRadians + sin(angleRadians)*cos(angleRadians))/(2*sin(angleRadians)*sin(angleRadians))

#
#   squared-length and reach:width factors...
#
    if(units=='metric')
      slFactor = .StemEnv$smpHectare/phi
    else
      slFactor = .StemEnv$sfpAcre/phi

    rwFactor = 1/tan(angleRadians/2)

    prs = new('pointRelascope', angleDegrees=angleDegrees, angleRadians=angleRadians,
              phi=phi, slFactor=slFactor, rwFactor=rwFactor,
              description=description, units=units
             )

    return(prs)
}   #pointRelascope constructor
)   #setMethod
       




          
#================================================================================
#  3. constructor method for class perpendicularDistance...
#
setMethod('perpendicularDistance',
          signature(kpds = 'numeric'),
function(kpds,
         units = 'metric',
         description = 'perpendicular distance method',
         ...
        )
{
#------------------------------------------------------------------------------
#
#   get the volume, surface area, or coverage factor...
#
    if(units=='metric')
      factor = .StemEnv$smpHectare/(2*kpds)
    else
      factor = .StemEnv$sfpAcre/(2*kpds)

    pds = new('perpendicularDistance', kpds=kpds, factor=factor, units=units,
              description=description
             )

    return(pds)
}   #perpendicularDistance constructor
)   #setMethod
    




          
#================================================================================
#  4. constructor method for class distanceLimited...
#
setMethod('distanceLimited',
          signature(distanceLimit = 'numeric'),
function(distanceLimit,
         units = 'metric',
         description = 'distance limited method',
         ...
        )
{
#------------------------------------------------------------------------------
#
    dl = new('distanceLimited',
             distanceLimit=distanceLimit,
             units=units,
             description=description
            )

    return(dl)
}   #distanceLimited constructor
)   #setMethod
    




          
#================================================================================
#  5. constructor method for class angleGauge...
#
setMethod('angleGauge',
          signature(baf = 'numeric'),
function(baf,
         units = 'metric',
         description = 'angle gauge method',
         ...
        )
{
#------------------------------------------------------------------------------
#
    if(units == .StemEnv$msrUnits$metric) {
      unitArea = .StemEnv$smpHectare
      conv = .StemEnv$m2cm
    }
    else {
      unitArea = .StemEnv$sfpAcre
      conv = .StemEnv$ft2in
    }
  
#
#   get the angle in radians and plot radius factor...
#
    angleRadians = 2*asin(sqrt(baf/unitArea))
    angleDegrees = .StemEnv$rad2Deg(angleRadians)
    diopters = 100 * tan(angleRadians)
    k = 2*sin(angleRadians/2)
    
#
#   points & lines...
#
    alpha = sqrt(unitArea/baf)
    PRF = alpha/2                                    #ft/ft or m/m
    prf = PRF/conv                                   #ft/in or m/cm

#
#   lines -- per ft or m segments...
#
    DF = sqrt(unitArea) * sqrt(baf)
    df = DF*12
    

    ag = new('angleGauge', angleDegrees=angleDegrees, angleRadians=angleRadians, diopters=diopters, k=k,
              baf=baf, prf=prf, PRF=PRF, alpha=alpha,
              df=df, DF=DF,
              description=description, units=units
             )

    return(ag)
}   #angleGauge constructor
)   #setMethod





#================================================================================
#  6. method for functions and class lineSegment...
#
setMethod('lineSegment',
          signature(length = 'numeric', orientation = 'numeric'),
function(length,
         orientation,              #in degrees
         units = 'metric',
         spUnits = CRS(projargs=as.character(NA)),
         centerPoint = c(x=0, y=0),   #centerPoint
         description = 'line segment',
         spID = paste('ls',.StemEnv$randomID(),sep=':'),
         ...
        )
{
#------------------------------------------------------------------------------
#
#   make sure the center is a named vector of length 2...
#
    if(any(is.na(match(c('x','y'),names(centerPoint)))))
      stop('Please use names x and y for centerPoint vector')
    if(length(centerPoint) != 2)
      stop('Please supply one set of (x,y) coordinates for the plot center location.')
 
  
#
#   some other checks...
#
    if(!is.numeric(length) || length <= 0)
      stop('length must be positive!')
    if(!is.numeric(orientation) || orientation < 0)
      stop('line orientation must be in [0, 360] degrees!')
    orientation = .StemEnv$deg2Rad(orientation)  #this will handle >360 reduction

#
#   now make the line segment in mattrix form first--heading due east...
#
    halfLen = length/2
    lineSeg = cbind(c(-halfLen, halfLen), c(0, 0), c(1,1))
    
#
#   rotate to north, then to correct position...
#
    rotAng = pi/2 - orientation

    trMat = transfMatrix(rotAng, centerPoint)
    lineSeg = lineSeg %*% trMat
    dimnames(lineSeg) = list(NULL,c('x','y','hc'))

#
#   and make a SpatialLines object...
#
    LineSeg = Line(lineSeg[,-3])                             #sans hc
    LinesSeg = Lines(list(LineSeg=LineSeg), ID=spID)
    spLinesSeg = SpatialLines(list(LinesSeg=LinesSeg),      #takes a list of Polygons objects
                                proj4string = spUnits                       
                               )
    
#
#   no id for center point, but it can be added to be the same as spID when
#   we make a container class for the center points elsewhere...
#
    loc = matrix(centerPoint, nrow=1)
    colnames(loc) = names(centerPoint)
    location = SpatialPoints(loc, proj4string = spUnits)

    ls = new('lineSegment', length=length, orientation=orientation, segment=spLinesSeg,
             description=description, units=units,
             location = location, spID=spID, spUnits=spUnits )

    return(ls)
}   #lineSegment constructor
)   #setMethod

