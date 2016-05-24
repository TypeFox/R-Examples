#---------------------------------------------------------------------------
#
#   This file holds the S4 definition for the constructors of downLog
#   classes via the downLog generic and methods.
#
#   The methods key off the argument 'data', which can be one of...
#     1. a data frame: for observed taper
#     2. measured variables like butt diameter, etc.; this will construct
#        a valid taper data frame and then call method 1 for the rest
#
#     3. a function: so one could pass one's own taper and volume functions
#     4. numeric: the butt and tip diameters & length?????????
#   **3. & 4. are not implemented, and may not need to be since a user can
#     make a data frame themselves from taper equations
#
#   The log's canonical position is centered at (0,0) with tip pointing east.
#   Translation and rotation are both relative to this position. Note in
#   particular that the log base/butt is not at (0,0), its center in (x,y)
#   is at (0,0).
#
#   Note that the sp package should be loaded for the complete functionality.
#
#   One little interesting idiosyncrasy is that when the tempfile() function
#   is used to generate log IDs for spatial polygons, they will not be keyed
#   to the .GlobalEnv random number stream, and so any comparison with identical
#   will yield FALSE for two given objects that were generated with the same
#   initial stream; therefore, I have elected to go back to simple random number
#   IDs below...
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
#   generic definition...
#
#if (!isGeneric("downLog")) 
  setGeneric('downLog',  
             function(object, ...) standardGeneric('downLog'),
             signature = c('object')
            )







          
#================================================================================
#  1. method for taper data frame object construction of class downLog...
#
setMethod('downLog',
          signature(object = 'data.frame'),
function(object,
         solidType = NULL,             #defaults to null for passed taper
         logAngle = 0,                 #canonical
         logVol = NULL,
         surfaceArea = NULL,
         coverageArea = NULL,
         biomass = NA,
         vol2wgt = NA,
         carbon = NA,
         wgt2carbon = NA,
         centerOffset = c(x=0, y=0),   #log center offset
         species = '',
         logID = paste('log',.StemEnv$randomID(),sep=':'),
         description = NULL,
         userExtra = NULL,
         units = 'metric',
         spUnits = CRS(projargs=as.character(NA)),
         runQuiet = TRUE,
         ...
        )
{
#------------------------------------------------------------------------------
#
#   some initial checks...
#
    if(any(is.na(match(c('x','y'),names(centerOffset)))))
      stop('Please use names x and y for centerOffset vector')
    if(length(centerOffset) < 2 || length(centerOffset) > 3)
      stop('Please supply one set of x,y[,z] coordinates for the object location.')
    if(is.na(match(units, .StemEnv$msrUnits)))
      stop('Illegal measurement units!')
    
#
#   check for a valid data frame since it will be required below, even though it
#   is redundantly checked in class validity checking...
#
    taper = object
    taperNames = match(colnames(taper), c('diameter','length'))
    if(any(is.na(taperNames)))
      stop('bad taper names--should be: "diameter" and "length"')
    nSegs = nrow(taper) - 1
    if(nSegs < 1)
      stop('there must be at least 2 rows in the taper data frame!')
    if(!isTRUE(all.equal(taper[1,'length'], 0.0)))
      stop('The bottom taper measurement of the log should have length=0')
    logLen = taper[nSegs+1,'length'] - taper[1,'length']
    if(logLen <= 0)
      stop('log length must be greater than zero!')

#
#   description field...
#
    if(is.null(description) || is.na(description))
      description = ''
    else
      description = as.character(description)

#
#   assign diameters, etc...
#
    buttDiam = taper[1, 'diameter']
    topDiam = taper[nSegs+1, 'diameter']
    solidType = solidType   
    
  
#
#   first use the object's validity check to check diameters, length, etc.
#   with dummy taper for now, etc...
#
    log = new('downLog', buttDiam=buttDiam, topDiam=topDiam, logLen=logLen,
              logAngle=logAngle, solidType=solidType, taper=taper, species=species
             )

#
#   okay, valid object to here, if the user passed log taper (solidType=NULL), then
#   use Smalian's to calculate volume if it is missing; otherwise, if solidType is
#   available, assume the taper came from the built-in equation...
#   
    if(is.null(logVol) || is.na(logVol)) {
      if(is.null(solidType))                         #user-defined taper
        logVol = .StemEnv$SmalianVolume(taper)$logVol
      else                                           #from default taper equation
        logVol = .StemEnv$wbVolume(buttDiam, topDiam, logLen, solidType)
    }

#
#   same type of test for surface area...
#
    if(is.null(surfaceArea) || is.na(surfaceArea)) {
      if(is.null(solidType))                        #user-defined taper
        surfaceArea = .StemEnv$splineSurfaceArea(taper, lenBot=0, lenTop=logLen)   #spline function
      else                                          #default taper equation
        surfaceArea = .StemEnv$wbSurfaceArea(buttDiam, topDiam, logLen, solidType, lenBot=0, lenTop=logLen)
    }

#
#   and again for coverage area...
#
    if(is.null(coverageArea) || is.na(coverageArea)) {
      if(is.null(solidType))                        #user-defined taper
        coverageArea = .StemEnv$splineCoverageArea(taper, lenBot=0, lenTop=logLen)   #spline function
      else                                          #default taper equation
        coverageArea = .StemEnv$wbCoverageArea(buttDiam, topDiam, logLen, solidType, lenBot=0, lenTop=logLen)
    }
    
#
#   biomass and carbon are a little different, they can be NA...
#
    if(is.na(biomass) && !is.na(vol2wgt))
      biomass = logVol * vol2wgt
    else if(!is.na(biomass) && is.na(vol2wgt))
      vol2wgt = biomass/logVol
    else if(is.na(biomass) && is.na(vol2wgt)) {   #cast to real NAs
      biomass = NA_real_
      vol2wgt = NA_real_
    }
    else if (!is.na(biomass) && !is.na(vol2wgt))
      stop('Please specify either biomass or vol2wgt, but not both!')
    if(!is.na(biomass) && is.na(carbon) && !is.na(wgt2carbon))
      carbon = biomass * wgt2carbon
    else if(!is.na(biomass) && !is.na(carbon) && is.na(wgt2carbon))
      wgt2carbon = carbon/biomass
    else if(is.na(carbon) && is.na(wgt2carbon)) {   #cast to real NAs
      carbon = NA_real_
      wgt2carbon = NA_real_
    }
    else if (!is.na(carbon) && !is.na(wgt2carbon))
      stop('Please specify either carbon or wgt2carbon, but not both!')
      

    
#
#   now, get the full profile for the outline as if the tree were standing upright
#   along the positive y-axis with the diameters center on x=0;
#   keep the duplicate point at the tip for truncation...
#
    rad = taper$diameter/2
    profile = data.frame( rad = c(-rad, rev(rad)) )          #center at diameter=x=0 
    profile$length = with(taper, c(length, rev(length)) )    #duplicate lengths too
    profile = rbind(profile, profile[1,])                    #close the polygon
    colnames(profile) = c('radius','length')
    np = nrow(profile)

#    
#   put the log on the ground, aligned with tip towards positive x-axis...
#
    halfLen = logLen/2
    rotLog = profile
    rotLog$hc = rep(1, np)                         #homogeneous coordinates
    rotLog = as.matrix(rotLog)
    trMat = transfMatrix(offset = c(0,-halfLen))   #to center at (0,0)
    rotLog = rotLog %*% trMat                      #translate
    trMat = transfMatrix(angle = -pi/2)            #lay the log down due East centered at (0,0)
    rotLog = rotLog %*% trMat

#
#   now rotate and translate to the desired position...
#
    trMat = transfMatrix(logAngle, offset = centerOffset) 
    rotLog = rotLog %*% trMat
    dimnames(rotLog) = list(NULL, c('x','y','hc'))
    needleAxis = matrix(c(-halfLen, halfLen, 0, 0, 1, 1), nrow=2) #row-major order
    needleAxis = needleAxis %*% trMat
    dimnames(needleAxis) = list(NULL, c('x','y','hc'))    

#
#   and make a SpatialPolygons object if sp is available...
#
    pgLog = Polygon(rotLog[,-3])        #sans hc
    pgsLog = Polygons(list(log=pgLog), ID=logID)
    spLog = SpatialPolygons(list(pgsLog=pgsLog),
                            proj4string = spUnits                       
                           )     #takes a list of Polygons objects

    
#
#   make the needle axis and the center point for the log...
#
    naLine = Line(needleAxis[,-3])  #no hc
    naLines = Lines(list(naLine=naLine), ID=paste(logID,'Needle',sep=':'))
    slNeedleAxis = SpatialLines(list(naLines=naLines),
                                proj4string = spUnits                                                    
                               )

    loc = matrix(centerOffset, nrow=1)
    colnames(loc) = names(centerOffset)
    location = SpatialPoints(loc,
                             proj4string = spUnits,        
                             bbox = bbox(spLog)              #same bbox as the log for consistency
                            )

#
#   now create it for real, it should be fine by this point...
#
    log = new('downLog',
              buttDiam = buttDiam,
              topDiam = topDiam,
              logLen = logLen,
              logAngle = logAngle,
              solidType = solidType,
              logVol = logVol,
              surfaceArea = surfaceArea,
              coverageArea = coverageArea,
              biomass = biomass,
              carbon = carbon,
              conversions = c(volumeToWeight=vol2wgt, weightToCarbon=wgt2carbon),
              taper = taper,
              profile = profile,
              rotLog = rotLog,
              spLog = spLog,
              slNeedleAxis = slNeedleAxis,
              units = units,
              spUnits = spUnits,
              location = location,
              description = description,
              userExtra = userExtra,
              species = species
             )

#
#   check for consistency between slot attribute values and corresponding taper values...
#
    checkStemDimensions(log)   #to within 1% by default

    if(!runQuiet)
      cat('\n')
    
    return(log)
}   #downLog method for data.frames
)   #setMethod






          
#================================================================================
#  2. method for using simple measurements in constructing class downLog...
#
#  This method just sets up the taper data frame from the measured variables and
#  passes things on to the data frame constructor from there...
#
setMethod('downLog',
          signature(object = 'missing'),
function(#object,
         buttDiam = 5,                 #cm
         topDiam = 0,                  #cm
         logLen = 5,                   #meters
         nSegs = 20,
         solidType = 3,                #defaults to 3  
         logAngle = 0,                 #canonical position
         logVol = NULL,
         surfaceArea = NULL,
         coverageArea = NULL,
         biomass = NA,
         vol2wgt = NA,
         carbon = NA,
         wgt2carbon = NA,
         centerOffset = c(x=0, y=0),   #log center offset
         species = '',
         logID = paste('log',.StemEnv$randomID(),sep=':'),
         description = NULL,
         userExtra = NULL,
         units = 'metric',
         spUnits = CRS(projargs=as.character(NA)),
         runQuiet = TRUE,
         ...
        )
{
#------------------------------------------------------------------------------
#
#   some initial checks...
#
    if(nSegs < 1)
      stop('A log must have at least one segment!')

#
#   convert diameters to meters or feet as appropriate...
#
    if(units == .StemEnv$msrUnits$English) {
      buttDiam = .StemEnv$in2ft * buttDiam
      topDiam = .StemEnv$in2ft * topDiam
    }
    else {
      buttDiam = .StemEnv$cm2m * buttDiam
      topDiam =  .StemEnv$cm2m * topDiam
    }

    
#
#   first use the object's validity check to check diameters, length, etc.
#   with dummy taper for now, etc...
#
    taper = data.frame(diameter=c(buttDiam, topDiam), length=c(0, logLen))
    log = new('downLog', buttDiam=buttDiam, topDiam=topDiam, logLen=logLen,
              logAngle=logAngle, solidType=solidType, taper=taper, species=species
             )

#
#   okay, valid object to here, just get the taper in terms of diameter and height
#   if it is missing...
#   
    taper = .StemEnv$wbTaper(buttDiam, topDiam, logLen, nSegs, solidType)
    if(is.null(logVol) || is.na(logVol))
      logVol = .StemEnv$wbVolume(buttDiam, topDiam, logLen, solidType)


#
#   all other checks will be made in the data frame constructor or during
#   validity checking on the proposed object...
#
    theLog = downLog(taper, logAngle = logAngle, logVol = logVol, solidType=solidType,
                     surfaceArea = surfaceArea, coverageArea = coverageArea,
                     biomass = biomass, vol2wgt = vol2wgt,
                     carbon = carbon, wgt2carbon = wgt2carbon,
                     centerOffset = centerOffset, species = species,
                     logID = logID, description = description,
                     userExtra = userExtra, units = units,
                     spUnits = spUnits, runQuiet = runQuiet,
                     ...
                    )

    return(theLog)
}   #downLog method for measured variables
)   #setMethod

 



#showMethods('downLog')
