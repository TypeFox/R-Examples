#---------------------------------------------------------------------------
#
#   This file holds the S4 class definitions for the stem-related classes.
#
#   Classes...
#   1. Stem: virtual -- parent of all stem objects
#   2. downLog -- for "normal" down logs
#   3. standingTree -- for "normal" standing trees
#
#   See page 361-362 in Chambers (2008) for information on validity checking--
#   specifically, by the time the methods below are called, we can assume
#   that all validity checks have been made on the object's slots and
#   superclass constraints.
#
#Author...									Date: 6-Aug-2010
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
#  1.  define the virtual Stem class...
#
setClass('Stem',
         
#
#  slots for the class and its subclasses...
#
    representation(species = 'character',          #species of piece                   
                   units = 'character',            #English or metric units
                   location = 'SpatialPoints',     #object "central" location
                   spUnits = 'CRS',                #sp units, character will change
#                  other...                   
                   description = 'character',      #descriptive comment
                   userExtra = 'ANY'               #anything else the user wants to include--no checks
                  ),
    prototype = list(species = '',               #some defaults for validity checking
                     units = 'metric',
                     location = SpatialPoints(matrix(c(0,0), nrow=1, dimnames=list('1',c('x','y'))) ),
                     spUnits = CRS(projargs=as.character(NA)),
                     description = '',
                     userExtra = NULL
                    ),
    contains = 'VIRTUAL',
    validity = function(object) {
                 if(!(object@units %in% .StemEnv$msrUnits))
                   return('units of measure must be "English" or "metric"')
                 
                 return(TRUE)
               } #validity check
) #class Stem






#=================================================================================================
#
#  2. the downed log class is just a direct descendant of 'Stem'...
#
#  inherited slots interpretation...
#    location = the center of the log longitudinally and radially
#
setClass('downLog',
    representation(buttDiam = 'numeric',          #diameter at the large end in same units as length
                   topDiam = 'numeric',           #diameter at the small end in same units as length
                   logLen = 'numeric',            #log length
                   logAngle = 'numeric',          #log lie angle from top pointing due East **RADIANS**
                   solidType = 'numericNULL',     #used in generic taper function
                   logVol = 'numeric',            #this would be total volume in cubic units
                   surfaceArea = 'numeric',       #total log surface area
                   coverageArea = 'numeric',      #projected coverage area
                   biomass = 'numeric',           #total log biomass--green or dry as desired
                   carbon = 'numeric',            #total log carbon
                   conversions = 'numeric',       #biomass & carbon conversions
                   taper = 'data.frame',          #diameter and length values
                   profile = 'data.frame',        #log profile of both sides, lying North, not rotated
                   rotLog = 'matrix',             #same as profile, but rotated and translated as desired
                   spLog = 'SpatialPolygons',     #rotated by logAngle and translated as desired
                   slNeedleAxis = 'SpatialLines'  #the central needle axis translated & rotated correctly
                  ),
    prototype = list(buttDiam = -1,               #these are dummy values to trigger initialize
                     topDiam = -1,
                     logLen = -1,
                     logAngle = 0,               #this is a valid angle
                     solidType = -1,
                     logVol = -1,
                     surfaceArea = -1,
                     coverageArea = -1,
                     biomass = -1,
                     carbon = -1
                    ),
    contains='Stem',                         #a subclass of the virtual 'Stem' class
    validity = function(object) {
                 if(object@topDiam > object@buttDiam)
                   return('butt diameter must be at least as large as top diameter!')
                 if(object@buttDiam <=0 || object@topDiam <0 || object@logLen <=0)
                   return('logLen, butt and top diameters must have some size!')
#                 if(object@logAngle < 0 || object@logAngle > 360)
#                   return('log angle must be 0 <= logAngle <- 360 degrees')
                 if(object@logAngle <0 || object@logAngle > 2*pi)
                   return('log angle must be 0 <= logAngle <- 2*pi radians')
                 if(!is.null(object@solidType) && (object@solidType < .StemEnv$solidTypes[1] ||
                                                   object@solidType > .StemEnv$solidTypes[2]))
                   return( paste('solidType must be in: [', .StemEnv$solidTypes[1], ',',
                                 .StemEnv$solidTypes[2], ']' ) )
                 
                 taperNames = match(colnames(object@taper), c('diameter','length'))
                 if(any(is.na(taperNames)))
                   return('taper names bad--woops!')
                 if(nrow(object@taper) < 2 || ncol(object@taper) < 2)
                   return('taper data frame must have at least two rows and columns!')
                 nr = nrow(object@taper)
                 #identical is too strong here...
                 if(!isTRUE(all.equal(object@buttDiam, object@taper[1,'diameter'])) ||
                    !isTRUE(all.equal(object@topDiam, object@taper[nr,'diameter'])) )
                   return('buttDiam and topDiam must be equal to their respective measurements in the taper data!')
                 if(!isTRUE(all.equal(object@logLen, object@taper[nr,'length']-object@taper[1,'length'])))
                   return('logLen must be equal to the total log length in taper data!')
                 if(object@logVol < 0)
                   return('negative logVol not allowed!')
                 if(object@surfaceArea < 0)
                   return('negative surfaceArea not allowed!')
                 if(object@coverageArea < 0)
                   return('negative coverageArea not allowed!')
                 if(!is.na(object@biomass) && object@biomass<0)
                   return('negative biomass not allowed!')
                 if(!is.na(object@carbon) && object@carbon<0)
                   return('negative carbon not allowed!')                 
                 conversionsNames = match(colnames(object@conversions), c('volumeToWeight','weightToCarbon'))
                 if(any(is.na(conversionsNames)))
                   return('profile names bad--woops!')
                 
                 profileNames = match(colnames(object@profile), c('radius','length'))
                 if(any(is.na(profileNames)))
                   return('profile names bad--woops!')
                 rotNames = match(colnames(object@rotLog), c('x','y', 'hc'))
                 if(any(is.na(rotNames)))
                   return('rotLog names bad--woops!')

                 if(!is.na(object@spUnits@projargs) && object@spUnits@projargs == '+proj=longlat')
                   return(paste('spUnits must be commensurate with units,',
                                'please convert to non-geographic coordinate system!')
                         )

                 
                 return(TRUE)
               } #validity check
) #class downLog


#-------------------------------------------------------------------------------
# this is just for very simple initialization with new, for more complex
# logs with taper measurements, use downLog constructor methods to initialize...
#
# initialize is called after the prototype values are set, so we can use those to
# set flags for default initialization here, before validity checking...
#
# note that we can not call downLog from below because it calls new(), which
# would get us into infinite recursion
#-------------------------------------------------------------------------------
#
setMethod('initialize', 'downLog',
  function(.Object, ...) {
    if(.Object@units == .StemEnv$msrUnits$English) {  #initializer puts these in meters or feet...
      if(.Object@buttDiam < 0)
        .Object@buttDiam = runif(1, 0.4, 2.6)  #feet
      if(.Object@logLen < 0)
        .Object@logLen = runif(1, 1, 30) 
    }
    else{
      if(.Object@buttDiam < 0)
        .Object@buttDiam = runif(1, 0.1, 0.8)  #meters
      if(.Object@logLen < 0)
        .Object@logLen = runif(1, 1, 10)
    }
    if(.Object@topDiam < 0)
      .Object@topDiam = runif(1, 0, 0.9)*.Object@buttDiam
    if(.Object@solidType < 1 || .Object@solidType > 10)
      .Object@solidType = round( runif(1,1,10) )

    .Object@logAngle = with( .StemEnv, runif(1,  logAngles[1], logAngles[2]) )#

    taper = as.data.frame(matrix(c(.Object@buttDiam, .Object@topDiam, 0, .Object@logLen), nrow=2))
    colnames(taper) = c('diameter','length')
    .Object@taper = taper

    if(.Object@logVol < 0)
      .Object@logVol = .StemEnv$wbVolume(.Object@buttDiam, .Object@topDiam, .Object@logLen, .Object@solidType)
    if(.Object@surfaceArea < 0)
      .Object@surfaceArea = .StemEnv$wbSurfaceArea(.Object@buttDiam, .Object@topDiam,
                                                   .Object@logLen, .Object@solidType,
                                                   0, .Object@logLen)
    if(.Object@coverageArea < 0)
      .Object@coverageArea = .StemEnv$wbCoverageArea(.Object@buttDiam, .Object@topDiam,
                                                     .Object@logLen, .Object@solidType,
                                                     0, .Object@logLen)
    if(.Object@biomass < 0)
      .Object@biomass = NA_real_

    if(.Object@carbon < 0)
      .Object@carbon = NA_real_

               
    callNextMethod(.Object, ...)
 } #function
) #setMethod initialize










#=================================================================================================
#
#  3. the standing tree class is just a direct descendant of 'Stem'...
#
#  inherited slots interpretation...
#    location = the center of the tree==pith at the base
#
setClass('standingTree',
    representation(buttDiam = 'numeric',          #diameter at the large end in same units as height
                   topDiam = 'numeric',           #diameter at the small end in same units as height
                   height = 'numeric',            #tree height
                   dbh = 'numeric',               #dbh in same units as height
                   ba = 'numeric',                #basal area
                   solidType = 'numericNULL',     #used in generic taper function
                   treeVol = 'numeric',           #this would be total volume in cubic units
                   surfaceArea = 'numeric',       #total tree surface area
                   #coverageArea = 'numeric',      #projected coverage area
                   biomass = 'numeric',           #total stem biomass--green or dry as desired
                   carbon = 'numeric',            #total stem carbon
                   conversions = 'numeric',       #biomass & carbon conversions
                   taper = 'data.frame',          #diameter and height values
                   profile = 'data.frame',        #tree profile of both sides standing
                   transTree = 'matrix',          #same as profile, but translated as desired
                   spTree = 'SpatialPolygons',    #outline of the full standing tree from profile
                   spDBH =  'SpatialPolygons'     #outline of dbh at location
                   #pith?? slNeedleAxis = 'SpatialLines'  #the central needle/pith axis 
                  ),
    prototype = list(buttDiam = -1,               #these are dummy values to trigger initialize
                     topDiam = -1,
                     height = -1,
                     dbh = NA_real_,
                     ba = NA_real_,
                     solidType = -1,
                     treeVol = -1,
                     surfaceArea = -1,
                     biomass = -1,
                     carbon = -1
                    ),
    contains='Stem',                         #a subclass of the virtual 'Stem' class
    validity = function(object) {
                 if(object@height < .StemEnv$dbhHgt[object@units] && !is.na(object@dbh))
                   return('tree has a dbh but total tree height is less than breast height!')
                 if(object@topDiam > object@buttDiam)
                   return('butt diameter must be at least as large as top diameter!')
                 if(object@buttDiam <=0 || object@topDiam <0 || object@height <=0)
                   return('height, butt and top diameters must have some size!')
                 if(!is.null(object@solidType) && (object@solidType < .StemEnv$solidTypes[1] ||
                                                   object@solidType > .StemEnv$solidTypes[2]))
                   return( paste('solidType must be in: [', .StemEnv$solidTypes[1], ',',
                                 .StemEnv$solidTypes[2], ']' ) )
                 
                 taperNames = match(colnames(object@taper), c('diameter','height'))
                 if(any(is.na(taperNames)))
                   return('taper names bad--woops!')
                 if(nrow(object@taper) < 2 || ncol(object@taper) < 2)
                   return('taper data frame must have at least two rows and columns!')
                 nr = nrow(object@taper)
                 #identical is too strong here, but even names in one and not the other (i.e., in buttDiam
                 #but not taper) will cause this to fail, so be careful!...
                 if(!isTRUE(all.equal(object@buttDiam, object@taper[1,'diameter'])) ||
                    !isTRUE(all.equal(object@topDiam, object@taper[nr,'diameter'])) )
                   return('buttDiam and topDiam must be equal to their respective measurements in the taper data!')
                 if(!isTRUE(all.equal(object@height, object@taper[nr,'height']-object@taper[1,'height'])))
                   return('height must be equal to the total tree height in taper data!')
                 if(object@treeVol < 0)
                   return('negative treeVol not allowed!')
                 if(object@surfaceArea < 0)
                   return('negative surfaceArea not allowed!')
                 if(!is.na(object@biomass) && object@biomass<0)
                   return('negative biomass not allowed!')
                 if(!is.na(object@carbon) && object@carbon<0)
                   return('negative carbon not allowed!')                 
                 conversionsNames = match(colnames(object@conversions), c('volumeToWeight','weightToCarbon'))
                 if(any(is.na(conversionsNames)))
                   return('profile names bad--woops!')
 
                 if( object@height < .StemEnv$dbhHgt[object@units] && !is.na(object@dbh) )
                   return(paste('tree is too short for a dbh, yet one has been recorded=',object@dbh))
                 #the following is possible with reverse taper and short trees--might want to be made a warning...
                 if(object@dbh > object@buttDiam || object@dbh < object@topDiam)
                   return('dbh is larger than butt diameter, or smaller than top diameter!')  
                    
                 
                 profileNames = match(colnames(object@profile), c('radius','height'))
                 if(any(is.na(profileNames)))
                   return('profile names bad--woops!')
                 transNames = match(colnames(object@transTree), c('x','y', 'hc'))
                 if(any(is.na(transNames)))
                   return('transTree names bad--woops!')

                 if(!is.na(object@spUnits@projargs) && object@spUnits@projargs == '+proj=longlat')
                   return(paste('spUnits must be commensurate with units,',
                                'please convert to non-geographic coordinate system!')
                         )
                 
                 
                 return(TRUE)
               } #validity check
) #class standingTree


#-------------------------------------------------------------------------------
# this is just for very simple initialization with new, for more complex
# trees with taper measurements, use standingTree constructor methods to initialize...
#
# initialize is called after the prototype values are set, so we can use those to
# set flags for default initialization here, before validity checking...
#
# note that we can not call standingTree from below because it calls new(), which
# would get us into infinite recursion
#-------------------------------------------------------------------------------
#
setMethod('initialize', 'standingTree',
  function(.Object, ...) {
    if(.Object@units == .StemEnv$msrUnits$English) {  #initializer puts these in meters or feet...
      if(.Object@buttDiam < 0)
        .Object@buttDiam = runif(1, 0.4, 2.6)  #feet
      if(.Object@height < 0)
        .Object@height = runif(1, 1, 30) 
    }
    else{
      if(.Object@buttDiam < 0)
        .Object@buttDiam = runif(1, 0.1, 0.8)  #meters
      if(.Object@height < 0)
        .Object@height = runif(1, 1, 20)
    }
    if(.Object@topDiam < 0)
      .Object@topDiam = runif(1, 0, 0.9)*.Object@buttDiam
    if(.Object@solidType < 1 || .Object@solidType > 10)
      .Object@solidType = round( runif(1,1,10) )


    taper = as.data.frame(matrix(c(.Object@buttDiam, .Object@topDiam, 0, .Object@height), nrow=2))
    colnames(taper) = c('diameter','height')
    .Object@taper = taper

    if(.Object@treeVol < 0)
      .Object@treeVol = .StemEnv$wbVolume(.Object@buttDiam, .Object@topDiam, .Object@height, .Object@solidType)
    if(.Object@surfaceArea < 0)
      .Object@surfaceArea = .StemEnv$wbSurfaceArea(.Object@buttDiam, .Object@topDiam,
                                                   .Object@height, .Object@solidType,
                                                   0, .Object@height)
    if(is.na(.Object@dbh))
      .Object@dbh = .StemEnv$wbTaper(.Object@buttDiam, .Object@topDiam,
                                                   .Object@height, 1, .Object@solidType,
                                                   .Object@height, FALSE)$diameter
    if(is.na(.Object@ba))
      .Object@ba = .StemEnv$baFactor[.Object@units]*.Object@dbh^2
    if(.Object@biomass < 0)
      .Object@biomass = NA_real_

    if(.Object@carbon < 0)
      .Object@carbon = NA_real_

               
    callNextMethod(.Object, ...)
 } #function
) #setMethod initialize
