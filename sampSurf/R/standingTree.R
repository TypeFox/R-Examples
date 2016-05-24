#---------------------------------------------------------------------------
#
#   This file holds the S4 definition for the constructors of standingTree
#   classes via the standingTree generic and methods.
#
#   The methods key off the argument 'data', which can be one of...
#     1. a data frame: for observed taper
#     2. measured variables like dbh, etc.; this will construct
#        a valid taper data frame and then call method 1 for the rest
#
#   These constructors are patterned very closely after the downLog
#   constructors functions, which were written first.
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
#   generic definition...
#
#if (!isGeneric("standingTree")) 
  setGeneric('standingTree',  
             function(object, ...) standardGeneric('standingTree'),
             signature = c('object')
            )







          
#================================================================================
#  1. method for taper data frame object construction of class standingTree...
#
setMethod('standingTree',
          signature(object = 'data.frame'),
function(object,
         solidType = NULL,             #defaults to null for passed taper
         treeVol = NULL,
         surfaceArea = NULL,
         biomass = NA,
         vol2wgt = NA,
         carbon = NA,
         wgt2carbon = NA,
         centerOffset = c(x=0, y=0),   #tree base-pith center offset
         species = '',
         treeID = paste('tree',.StemEnv$randomID(),sep=':'),
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
    taperNames = match(colnames(taper), c('diameter','height'))
    if(any(is.na(taperNames)))
      stop('bad taper names--should be: "diameter" and "height"')
    nSegs = nrow(taper) - 1
    if(nSegs < 1)
      stop('there must be at least 2 rows in the taper data frame!')
    if(!isTRUE(all.equal(taper[1,'height'], 0.0)))
      stop('The bottom taper measurement of the tree should have length=0')
    height = taper[nSegs+1,'height'] - taper[1,'height']
    if(height <= 0)
      stop('tree height must be greater than zero!')

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
    if(height < .StemEnv$dbhHgt[units])
      dbh = NA_real_
    else
      if(is.null(solidType))                                   #user-defined taper
        dbh = spline(taper$height, taper$diameter, xout=.StemEnv$dbhHgt[units])$y
      else
        dbh = .StemEnv$wbTaper(buttDiam, topDiam, height, nSegs, solidType, .StemEnv$dbhHgt[units],
                               isLog=FALSE)$diameter
  
#
#   first use the object's validity check to check diameters, height, etc.
#
    tree = new('standingTree', buttDiam=buttDiam, topDiam=topDiam, dbh=dbh, height=height,
               solidType=solidType, taper=taper, species=species
              )

#
#   okay, valid object to here, if the user passed tree taper (solidType=NULL), then
#   use Smalian's to calculate volume if it is missing; otherwise, if solidType is
#   available, assume the taper came from the built-in equation...
#   
    if(is.null(treeVol) || is.na(treeVol)) {
      if(is.null(solidType))                                   #user-defined taper
        treeVol = .StemEnv$SmalianVolume(taper, FALSE)$logVol  #logVol is treeVol here
      else                                                     #from default taper equation
        treeVol = .StemEnv$wbVolume(buttDiam, topDiam, height, solidType)
    }

#
#   same type of test for surface area...
#
    if(is.null(surfaceArea) || is.na(surfaceArea)) {
      if(is.null(solidType))                        #user-defined taper
        surfaceArea = .StemEnv$splineSurfaceArea(taper, lenBot=0, lenTop=height, isLog=FALSE)   #spline function
      else                                          #default taper equation
        surfaceArea = .StemEnv$wbSurfaceArea(buttDiam, topDiam, height, solidType, lenBot=0, lenTop=height)
    }
    
#
#   biomass and carbon are a little different, they can be NA...
#
    if(is.na(biomass) && !is.na(vol2wgt))
      biomass = treeVol * vol2wgt
    else if(!is.na(biomass) && is.na(vol2wgt))
      vol2wgt = biomass/treeVol
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
#   now, get the full profile for the outline for the tree standing upright
#   along the positive y-axis with the diameters center on x=0;
#   keep the duplicate point at the tip for truncation...
#
    rad = taper$diameter/2
    profile = data.frame( rad = c(-rad, rev(rad)) )          #center at diameter=x=0; assume straight 
    profile$height = with(taper, c(height, rev(height)) )    #duplicate heights too
    profile = rbind(profile, profile[1,])                    #close the polygon
    colnames(profile) = c('radius','height')
    np = nrow(profile)


#
#   now translate to the desired position...
#
    transTree = profile
    transTree$hc = rep(1, np)                         #homogeneous coordinates
    transTree = as.matrix(transTree)
    trMat = transfMatrix(offset = centerOffset) 
    transTree = transTree %*% trMat
    dimnames(transTree) = list(NULL, c('x','y','hc'))

#
#   and make a SpatialPolygons object if sp is available...
#
    pgTree = Polygon(transTree[,-3])        #sans hc
    pgsTree = Polygons(list(tree = pgTree), ID = treeID)
    spTree = SpatialPolygons(list(pgsTree = pgsTree),
                             proj4string = spUnits                       
                             )     #takes a list of Polygons objects

#
#   now for dbh--assume circular cross-section...
#
    sp.dbh = spCircle(radius=dbh/2, spUnits=spUnits, centerPoint=centerOffset,
                      spID=paste('dbh',treeID,sep='.'), ...)
    spDBH = sp.dbh$spCircle
    names(spDBH@polygons) = 'pgsDBH'
    names(spDBH@polygons$pgsDBH@Polygons) = 'pgDBH'

#
#   now create it for real, it should be fine by this point...
#
    tree = new('standingTree',
               buttDiam = buttDiam,
               topDiam = topDiam,
               height = height,
               dbh = dbh,
               ba = .StemEnv$baFactor[units]*dbh*dbh,
               solidType = solidType,
               treeVol = treeVol,
               surfaceArea = surfaceArea,
               biomass = biomass,
               carbon = carbon,
               conversions = c(volumeToWeight=vol2wgt, weightToCarbon=wgt2carbon),
               taper = taper,
               profile = profile,
               transTree = transTree,
               spTree = spTree,
               spDBH = spDBH,
               units = units,
               spUnits = spUnits,
               location = sp.dbh$location,
               description = description,
               userExtra = userExtra,
               species = species
              )

#
#   check for consistency between slot attribute values and corresponding taper values...
#
    checkStemDimensions(tree)   #to within 1% by default

    
    return(tree)
}   #standingTree method for data.frames
)   #setMethod






          
#================================================================================
#  2. method for using simple measurements in constructing class standingTree...
#
#  This method just sets up the taper data frame from the measured variables and
#  passes things on to the data frame constructor from there...
#
setMethod('standingTree',
          signature(object = 'missing'),
function(#object,
         dbh = 20,                     #cm
         topDiam = 0,                  #cm
         height = 15,                  #meters
         nSegs = 20,
         solidType = 3,                #must have some taper model for butt diam
         treeVol = NULL,
         surfaceArea = NULL,
         biomass = NA,
         vol2wgt = NA,
         carbon = NA,
         wgt2carbon = NA,
         centerOffset = c(x=0, y=0),   #tree base-pith center offset
         species = '',
         treeID = paste('tree',.StemEnv$randomID(),sep=':'),
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
      stop('A tree must have at least one segment!')
    if(is.null(solidType) || is.na(solidType))
      stop('For simple tree measurements, a valid solidType must be specified!')
    if(dbh <=0 || dbh < topDiam)
      stop('dbh must be a real positive value at least as large as topDiam!')
    if(height < .StemEnv$dbhHgt[units] && !is.na(dbh))
      stop('tree has a dbh but total tree height is less than breast height!')

#
#   convert diameters to meters or feet as appropriate...
#
    if(units == .StemEnv$msrUnits$English) {
      topDiam = .StemEnv$in2ft * topDiam
      dbh = .StemEnv$in2ft * dbh
    }
    else {
      topDiam =  .StemEnv$cm2m * topDiam
      dbh =  .StemEnv$cm2m * dbh
    }

#
#   back butt diameter out from the default taper equation...
#
    H.d = ((height - .StemEnv$dbhHgt[units])/height)^(2/solidType)
    buttDiam = (topDiam*(H.d-1) + dbh)*H.d^{-1}      #will be in correct height units here!
    buttDiam = as.numeric(buttDiam)  #strip id from .StemEnv$dbhHgt conversion above for new() below

#
#   first use the object's validity check to check diameters, length, etc.
#   with dummy taper for now, etc...
#
    taper = data.frame(diameter=c(buttDiam, topDiam), height=c(0, height))
    log = new('standingTree', buttDiam=buttDiam, topDiam=topDiam, height=height, 
              dbh=dbh, solidType=solidType, taper=taper, species=species
             )

#
#   okay, valid object to here, just get the taper in terms of diameter and height
#   if it is missing...
#   
    taper = .StemEnv$wbTaper(buttDiam, topDiam, height, nSegs, solidType, isLog=FALSE)
    if(is.null(treeVol) || is.na(treeVol))
      treeVol = .StemEnv$wbVolume(buttDiam, topDiam, height, solidType)

#
#   all other checks will be made in the data frame constructor or during
#   validity checking on the proposed object...
#
    theTree = standingTree(taper, solidType=solidType, treeVol=treeVol,
                           surfaceArea = surfaceArea, biomass = biomass, 
                           vol2wgt = vol2wgt, carbon = carbon,  wgt2carbon = wgt2carbon,
                           centerOffset = centerOffset, species = species,
                           treeID = treeID, description = description, userExtra = userExtra,
                           units = units, spUnits = spUnits, runQuiet = runQuiet,
                           ...
                          )
    
    return(theTree)
}   #standingTree method for measured variables
)   #setMethod
