#---------------------------------------------------------------------------
#
#   Methods for generic summary() for class...
#     (1) Stem virtual class
#     (2) downLog subclass
#     (3) standingTree subclass
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



#================================================================================
#  1. method for class Stem...
#
setMethod('summary',
          signature(object = 'Stem'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   just a simple summary of common items from virtual class...
#------------------------------------------------------------------------------
    cat('\nObject of class:', class(object))
    .StemEnv$underLine(60)
    if(!is.na(object@description))
      cat(object@description, fill=60)
    .StemEnv$underLine(60, prologue='')

    cat('\nStem...')
    cat('\n  Species: ', object@species)
    cat('\n  units of measurement: ', object@units)
    cat('\n  spatial units: ', object@spUnits@projargs)
    cat('\n  location...')
    cat('\n    x coord: ', coordinates(object@location)[,'x'])
    cat('\n    y coord: ', coordinates(object@location)[,'y'])

    cat('\n')
    
    return(invisible())
}   #summary for 'Stem'
) #setMethod




#================================================================================
#  2. method for subclass "downLog"...
#
setMethod('summary',
          signature(object = 'downLog'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to 'Stem' method for 'downLog'...
#------------------------------------------------------------------------------
    callNextMethod()
    cat('    (Above coordinates are for log center)')
    cat('\n  Spatial ID:', object@spLog@polygons$pgsLog@ID)

    cat('\n\ndownLog...')
    if(object@units == .StemEnv$msrUnits$metric) {
      cat('\n  Butt diameter = ', object@buttDiam, ' meters (', object@buttDiam*.StemEnv$m2cm, ' cm)',sep='')
      cat('\n  Top diameter = ', object@topDiam, ' meters (', object@topDiam*.StemEnv$m2cm, ' cm)', sep='')
      cat('\n  Log length =', object@logLen, 'meters')
      cat('\n  Log volume =', object@logVol, 'cubic meters')
      cat('\n  Log surface area =', object@surfaceArea, 'square meters')
      cat('\n  Log coverage area =', object@coverageArea, 'square meters')
    }
    else {
      cat('\n  Butt diameter = ', object@buttDiam, ' feet (', object@buttDiam*.StemEnv$ft2in, ' in)',sep='')
      cat('\n  Top diameter = ', object@topDiam, ' feet (', object@topDiam*.StemEnv$ft2in, ' in)', sep='')
      cat('\n  Log length =', object@logLen, 'feet')
      cat('\n  Log volume =', object@logVol, 'cubic feet')
      cat('\n  Log surface area =', object@surfaceArea, 'square feet')
      cat('\n  Log coverage area =', object@coverageArea, 'square feet')
    }
    if(!is.na(object@biomass))
      cat('\n  Log biomass =', object@biomass)
    if(!is.na(object@carbon))
      cat('\n  Log carbon =', object@carbon)
    if(!is.na(object@conversions['volumeToWeight']))
      cat('\n  Volume to weight conversion =', object@conversions['volumeToWeight'])
    if(!is.na(object@conversions['weightToCarbon']))
      cat('\n  Weight to carbon conversion =', object@conversions['weightToCarbon'])

    
    cat('\n  Log angle of lie =', object@logAngle, 'radians')
    cat(' (', .StemEnv$rad2Deg(object@logAngle), ' degrees)', sep='')
    #cat('\n  Log angle of lie = ', .StemEnv$rad2Deg(object@logAngle), 'degrees')
    cat('\n  Taper parameter = ', ifelse(is.null(object@solidType), 'NULL', object@solidType) )
    
    cat('\n\nTaper (in part)...\n')
    print(head(object@taper))

    if(!is.null(object@userExtra))
      cat('\n  "Note: userExtra" slot is non-NULL')
    
#
#   important check to see if any valid SpatialPolygon exists for the object...
#
    if(length(object@spLog@polygons) == 0)  #check for object made with new()
      cat('\n\n***No spLog "SpatialPolygons" -- please use downLog constructor!\n')

    cat('\n')
        
    return(invisible())
}   #summary for 'downLog'
) #setMethod







#================================================================================
#  3. method for subclass "standingTree"...
#
setMethod('summary',
          signature(object = 'standingTree'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to 'Stem' method for 'standingTree'...
#------------------------------------------------------------------------------
    callNextMethod()
    cat('    (Above coordinates are for dbh center)')
    cat('\n  Spatial ID:', object@spTree@polygons$pgsTree@ID)

    cat('\n\nstandingTree...')
    if(object@units == .StemEnv$msrUnits$metric) {
      cat('\n  Butt diameter = ', object@buttDiam, ' meters (', object@buttDiam*.StemEnv$m2cm, ' cm)',sep='')
      cat('\n  Top diameter = ', object@topDiam, ' meters (', object@topDiam*.StemEnv$m2cm, ' cm)', sep='')
      cat('\n  DBH = ', object@dbh, ' meters (', object@dbh*.StemEnv$m2cm, ' cm)', sep='')
      cat('\n  Basal area = ', object@ba, 'square meters')
      cat('\n  Height =', object@height, 'meters')
      cat('\n  Tree volume =', object@treeVol, 'cubic meters')
      cat('\n  Tree surface area =', object@surfaceArea, 'square meters')
    }
    else {
      cat('\n  Butt diameter = ', object@buttDiam, ' feet (', object@buttDiam*.StemEnv$ft2in, ' in)',sep='')
      cat('\n  Top diameter = ', object@topDiam, ' feet (', object@topDiam*.StemEnv$ft2in, ' in)', sep='')
      cat('\n  DBH = ', object@dbh, ' feet (', object@dbh*.StemEnv$ft2in, ' in)', sep='')
      cat('\n  Basal area = ', object@ba, 'square feet')
      cat('\n  Height =', object@height, 'feet')
      cat('\n  Tree volume =', object@treeVol, 'cubic feet')
      cat('\n  Tree surface area =', object@surfaceArea, 'square feet')
    }
    if(!is.na(object@biomass))
      cat('\n  Tree biomass =', object@biomass)
    if(!is.na(object@carbon))
      cat('\n  Tree carbon =', object@carbon)
    if(!is.na(object@conversions['volumeToWeight']))
      cat('\n  Volume to weight conversion =', object@conversions['volumeToWeight'])
    if(!is.na(object@conversions['weightToCarbon']))
      cat('\n  Weight to carbon conversion =', object@conversions['weightToCarbon'])
    cat('\n  Taper parameter = ', ifelse(is.null(object@solidType), 'NULL', object@solidType) )
    
    cat('\n\nTaper (in part)...\n')
    print(head(object@taper))

    if(!is.null(object@userExtra))
      cat('\n  "Note: userExtra" slot is non-NULL')
    
#
#   important check to see if any valid SpatialPolygon exists for the object...
#
    if(length(object@spTree@polygons) == 0)  #check for object made with new()
      cat('\n\n***No spTree "SpatialPolygons" -- please use standingTree constructor!\n')
    if(length(object@spDBH@polygons) == 0)  #check for object made with new()
      cat('\n\n***No spDBH "SpatialPolygons" -- please use standingTree constructor!\n')

    cat('\n')
        
    return(invisible())
}   #summary for 'standingTree'
) #setMethod
    



#showMethods('summary')
