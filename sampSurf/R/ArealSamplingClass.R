#---------------------------------------------------------------------------
#
#   This file holds the S4 class definitions for the Areal Sampling method
#   related classes.
#
#   Classes...
#     1. ArealSampling: virtual class for all methods
#     2. circularPlot: class for fixed-radius circular plot sampling
#     3. pointRelascope: class for point relascope sampling (PRS) (Jan 2011)
#     4. perpendicularDistance: class for PDS (Jan 2011)
#     5. distanceLimit: class for variable plot MC (Mar 2011)
#     6. angleGauge: class for Bitterlich sampling methods
#     7. lineSegment: class for numerous line-based sampling methods (3-Oct-2012)
#
#Author...									Date: 19-Aug-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#




#=================================================================================================
#
#  define the virtual ArealSampling class...
#
setClass('ArealSampling',
         
#
#  slots for the class and its subclasses...
#
    representation(description = 'character',      #more descriptive name
                   units = 'character'             #English or metric units
                  ),
    prototype = list(description = '',                    #some defaults for validity checking
                     units = 'metric'
                    ),
    contains = 'VIRTUAL',
    validity = function(object) {
                 if(!(object@units %in% c('English','metric')))
                   return('units of measure must be "English" or "metric"')

                 return(TRUE)
               } #validity check
) #class ArealSampling
#=================================================================================================











#=================================================================================================
#
#  the circular plot class is just a direct descendant of 'ArealSampling'...
#
#    location = the center of the circle
#
setClass('circularPlot',
    representation(radius = 'numeric',              #plot radius
                   area = 'numeric',                #plot area
                   perimeter = 'SpatialPolygons',   #sp polygon of the plot perimeter
                   location = 'SpatialPoints',      #plot center location
                   spID = 'character',              #short id name for polygon label
                   spUnits = 'CRS'                  #sp units, character will change
                  ),
    prototype = list(radius = -1,                    #
                     area = 0,
                     location = SpatialPoints(matrix(c(0,0), nrow=1, dimnames=list('1',c('x','y'))) ),
                     spID = paste('cp', format(runif(1, 0,10000),digits=8), sep=':'),
                     spUnits = CRS(projargs=as.character(NA)) 
                    ),
    contains = 'ArealSampling',                     #a subclass of the virtual 'ArealSampling' class
    validity = function(object) {
                 if(object@radius <= 0)
                   return('plot radius must be positive!')
                 if(object@area <= 0)
                   return('plot area must be positive!')
                 
                 locNames = match(colnames(object@location), c('x','y'))
                 if(any(is.na(locNames)))
                   return('location names must be "x" and "y"!')
      
                 if(!is.na(object@spUnits@projargs) && object@spUnits@projargs == '+proj=longlat')
                   return(paste('spUnits must be commensurate with units,',
                                'please convert to non-geographic coordinate system!')
                         )
      
                 return(TRUE)
               } #validity check
) #class circularPlot
         

#
# initialize is called after the prototype values are set, so we can use them to
# set flags for default initialization here, before validity checking...
#
setMethod('initialize', 'circularPlot',
  function(.Object, ...) {

    if(.Object@radius < 0)
      .Object@radius = runif(1, 1, 20)  #meters
    .Object@area = pi * .Object@radius * .Object@radius

               
    callNextMethod(.Object, ...)
 } #function
) #setMethod initialize circularPlot
#=================================================================================================








#=================================================================================================
#
#  the point relascope class is just a direct descendant of 'ArealSampling'...
#
#
setClass('pointRelascope',
    representation(angleDegrees = 'numeric',        #relascope angle in degrees
                   angleRadians = 'numeric',        #relascope angle in radians
                   phi = 'numeric',
                   slFactor = 'numeric',            #squared-length factor
                   rwFactor = 'numeric'             #reach:width factor (width=1 always)
                  ),
    prototype = list(angleDegrees = 90,
                     angleRadians = pi/2,
                     phi = 0.78539816,
                     slFactor = 55462.315,
                     rwFactor = 1
                    ),
    contains = 'ArealSampling',                     #a subclass of the virtual 'ArealSampling' class
    validity = function(object) {
                 if(object@angleDegrees <= 0 || object@angleDegrees > 90)
                   return('Relascope angle must be between 0 and 90 degrees!')
                 
                 return(TRUE)
               } #validity check
) #class pointRelascope








#=================================================================================================
#
#  the perpendicular distance class is just a direct descendant of 'ArealSampling'...
#
#
setClass('perpendicularDistance',
    representation(factor = 'numeric',              #volume, surface, or coverage factor in appropriate units
                   kpds = 'numeric'                 #pds factor in appropriate units
                  ),
    prototype = list(kpds = 43.56,                  #e.g., for volume: ft^{-1}
                     factor = 500,                  #e.g., for volume: ft^3/acre
                     units = 'English',
                     description = 'English PDS'
                    ),
    contains = 'ArealSampling',                     #a subclass of the virtual 'ArealSampling' class
    validity = function(object) {
                 if(object@kpds <= 0 || object@factor <= 0)
                   return('PDS factors must be greater than zero!')
                 
                 return(TRUE)
               } #validity check
) #class perpendicularDistance







#=================================================================================================
#
#  the distance limited class is just a direct descendant of 'ArealSampling'...
#
#
setClass('distanceLimited',
    representation(distanceLimit = 'numeric'        #the limiting distance
                  ),
    prototype = list(distanceLimit = 1, 
                     units = 'English',
                     description = 'English DL'
                    ),
    contains = 'ArealSampling',                     #a subclass of the virtual 'ArealSampling' class
    validity = function(object) {
                 if(is.na(object@distanceLimit) || object@distanceLimit <= 0)
                   return('distance limit  must be greater than zero!')
                 
                 return(TRUE)
               } #validity check
) #class distanceLimited





#=================================================================================================
#
#  the angle guage class is just a direct descendant of 'ArealSampling'...
#
#  this is for small angle sampling methods like horizontal point sampling for standing
#  trees--not point or transect relascope methods for logs
#
#
setClass('angleGauge',
    representation(angleDegrees = 'numeric',        #gauge angle in degrees
                   angleRadians = 'numeric',        #gauge angle in radians
                   diopters = 'numeric',            #prism angle diopters (Delta)
                   k = 'numeric',                   #gauge constant (dimensionless)
                   prf = 'numeric',                 #plot radius factor (ft/in or m/cm)
                   PRF = 'numeric',                 #prf (ft/ft or m/m) as in tree dbh units
                   alpha = 'numeric',               #proportionality constanty (ft/ft or m/m)
                   #points...
                   baf = 'numeric',                 #basal area factor (ft^2/ac or m^2/ha)
                   #lines: based on 1ft or 1m segments...
                   df = 'numeric',                  #diameter factor in inches or cm
                   DF = 'numeric'                   #diameter factor in ft or m
                  ),
    prototype = list(angleDegrees = 2.29,
                     angleRadians = 0.040,
                     baf = 4,                       #metric
                     prf = 0.25,
                     alpha = 50
                    ),
    contains = 'ArealSampling',                     #a subclass of the virtual 'ArealSampling' class
    validity = function(object) {
                 maxDegrees = .StemEnv$angleGaugeMaxDegrees
                 if(object@angleDegrees <= 0 || object@angleDegrees > maxDegrees)
                   return(paste('Gauge angle must be between 0 and',maxDegrees,'degrees!'))
                 
                 return(TRUE)
               } #validity check
) #class angleGauge







#=================================================================================================
#
#  the lineSegment class is just a direct descendant of 'ArealSampling'...
#
#  potentially, the lineSegment class might be used to build more complicated structures like
#  "L" or "Y" shaped transects sometimes used in LIS
#
setClass('lineSegment',
    representation(orientation = 'numeric',         #line orientation clockwise from North==0 **RADIANS**
                   length = 'numeric',              #line segment length
                   segment = 'SpatialLines',        #sp lines of the segment
                   location = 'SpatialPoints',      #line segment center location
                   spID = 'character',              #short id name for polygon label
                   spUnits = 'CRS'                  #sp units, character will change
                  ),
    prototype = list(orientation = -1,              #we don't want this to be a valid object here,
                     length = 0,                    #let initialize assign...
                     location = SpatialPoints(matrix(c(0,0), nrow=1, dimnames=list('1',c('x','y'))) ),
                     spID = paste('ls', format(runif(1, 0,10000),digits=8), sep=':'),
                     spUnits = CRS(projargs=as.character(NA)) 
                    ),
    contains = 'ArealSampling',                     #a subclass of the virtual 'ArealSampling' class
    validity = function(object) {
                 if(object@length <= 0)
                   return('line segment length must be positive!')
                 if(object@orientation < 0 || object@orientation > 2*pi)
                   return('line segment orientation must be within [0, 2*pi] radians ([0,360] degrees)!')
                 
                 locNames = match(colnames(object@location), c('x','y'))
                 if(any(is.na(locNames)))
                   return('location names must be "x" and "y"!')
      
                 if(!is.na(object@spUnits@projargs) && object@spUnits@projargs == '+proj=longlat')
                   return(paste('spUnits must be commensurate with units,',
                                'please convert to non-geographic coordinate system!')
                         )
      
                 return(TRUE)
               } #validity check
) #class lineSegment
 
         
#
# initialize is called after the prototype values are set, so we can use them to
# set flags for default initialization here, before validity checking...
#
setMethod('initialize', 'lineSegment',
  function(.Object, ...) {

    if(.Object@length <= 0)
      .Object@length = runif(1, 1, 10)  #meters
    if(.Object@orientation < 0 || .Object@orientation > 2*pi)
      .Object@orientation = runif(1)*2*pi

               
    callNextMethod(.Object, ...)
 } #function
) #setMethod initialize lineSegment
#=================================================================================================
