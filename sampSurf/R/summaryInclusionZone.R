#---------------------------------------------------------------------------
#
#   Methods for generic summary() for InclusionZone class...
#     (1) InclusionZone base class
#
#         ...downLogIZ subclasses...
#     (2) downLog component class--on a per unit area basis
#     (3) standUpIZ class
#     (4) chainSawIZ class
#     (5) sausageIZ class
#     (6) pointRelscopeIZ class
#     (7) perpendicularDistanceIZ class
#     (8) distanceLimitedPDSIZ class
#     (9) distanceLimitedIZ class
#
#         ...standingTreeIZ subclasses...
#      1. standingTreeIZ
#      2. circularPlotIZ
#      3. horizontalPointIZ
#
#Author...									Date: 24-Aug-2010
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
#  method for virtual class InclusionZone...
#
setMethod('summary',
          signature(object = 'InclusionZone'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   just a simple summary of common items from virtual class...
#------------------------------------------------------------------------------
    cat('\nObject of class:', class(object))
    if(.hasSlot(object, 'antithetic'))                    #for Monte Carlo Sampling derivatives
      if(object@antithetic)                               #distinguish antithetic variants
        cat(' (antithetic)')
    .StemEnv$underLine(60)
    if(!is.na(object@description))
      cat(object@description, fill=60)
    .StemEnv$underLine(60, prologue='')

    cat('\nInclusionZone...')
    cat('\n  Units of measurement: ', object@units)    
    cat('\n  Per unit area blowup factor:', object@puaBlowup)
    cat(ifelse(object@units == .StemEnv$msrUnits$metric, ' per hectare', ' per acre'))
    cat('\n\n  Object bounding box...\n');print(object@bbox)

    if(!is.null(object@userExtra))
      cat('\n  "Note: userExtra" slot is non-NULL')

    #cat('\n')
    
    return(invisible())
}   #summary for 'InclusionZone'
) #setMethod





#================================================================================
#  method for class "downLogIZ"...
#
setMethod('summary',
          signature(object = 'downLogIZ'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to 'InclusionZone' method for 'downLogIZ'...
#------------------------------------------------------------------------------
    callNextMethod()

    cat('\ndownLog component estimates...')
    cat('\n  Spatial ID:', object@downLog@spLog@polygons$pgsLog@ID)
    if(object@units == .StemEnv$msrUnits$metric) {
      cat('\n  Number of logs:', object@puaEstimates$Density, 'per hectare')
      cat('\n  Volume:', object@puaEstimates$volume, 'cubic meters per hectare')
      cat('\n  Surface area:', object@puaEstimates$surfaceArea, 'square meters per hectare')
      cat('\n  Coverage area:', object@puaEstimates$coverageArea, 'square meters per hectare')
      cat('\n  Length:', object@puaEstimates$Length, 'meters per hectare')
      cat('\n  Biomass (woody):', object@puaEstimates$biomass, 'per hectare')
      cat('\n  Carbon content:', object@puaEstimates$carbon, 'per hectare')
    }
    else {
      cat('\n  Number of logs:', object@puaEstimates$Density, 'per acre')
      cat('\n  Volume:', object@puaEstimates$volume, 'cubic feet per acre')
      cat('\n  Surface area:', object@puaEstimates$surfaceArea, 'square feet per acre')
      cat('\n  Coverage area:', object@puaEstimates$coverageArea, 'square feet per acre')
      cat('\n  Length:', object@puaEstimates$Length, 'feet per acre')
      cat('\n  Biomass (woody):', object@puaEstimates$biomass, 'per acre')
      cat('\n  Carbon content:', object@puaEstimates$carbon, 'per acre')
    }

    cat('\n')
    
    return(invisible())
}   #summary for 'downLogIZ'
) #setMethod


    



#================================================================================
#  method for class "standUpIZ"...
#
setMethod('summary',
          signature(object = 'standUpIZ'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to 'InclusionZone' & 'downLogIZ' methods for 'standUpIZ'...
#------------------------------------------------------------------------------
    callNextMethod()
    
    cat('\nstandUpIZ...')
    cat('\n  use \"summary\" on the circularPlot slot for details')
    #summary(object@circularPlot)

    cat('\n')
    
    return(invisible())
}   #summary for 'standUpIZ'
) #setMethod

    



#================================================================================
#  method for class "chainSawIZ"...
#
setMethod('summary',
          signature(object = 'chainSawIZ'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to 'InclusionZone' & 'downLogIZ' methods for 'chainSawIZ'...
#------------------------------------------------------------------------------
    callNextMethod()
    cat('  The above estimates are based on the expanded sliver portion.\n\n')
    cat('  The following are unexpanded...')
    if(object@units == .StemEnv$msrUnits$metric) {
      cu = 'cubic meters'
      su = 'square meters'
      lu = 'meters'
    }
    else {
      cu = 'cubic feet'
      su = 'square feet'
      lu = 'feet'
    }
    
    cat('\n    Sliver area is',object@bolt$area['propArea'],'of bounding bolt')
    cat('\n    Total log volume:', object@downLog@logVol, cu)
    cat('\n      Bounding bolt volume:',object@bolt$boltVol, cu)
    cat('\n      Sliver volume:', object@bolt$sectVol, cu)
    cat('\n    Total log surface area:', object@downLog@surfaceArea, su)
    cat('\n      Bounding bolt surface area:',object@bolt$boltSA, su)
    cat('\n      Sliver surface area:', object@bolt$sectSA, su)
    cat('\n    Total log coverage area:', object@downLog@coverageArea, su)
    cat('\n      Bounding bolt coverage area:',object@bolt$boltCA, su)
    cat('\n      Sliver coverage area:', object@bolt$sectCA, su)
    cat('\n    Total log length:', object@downLog@logLen, lu)
    cat('\n      Bounding bolt length:',object@bolt$boltLen, lu)
    cat('\n    Total log biomass:', object@downLog@biomass)
    cat('\n      Bounding bolt biomass:',object@bolt$boltBms)
    cat('\n      Sliver biomass:', object@bolt$sectBms)
    cat('\n    Total log carbon:', object@downLog@carbon)
    cat('\n      Bounding bolt carbon:',object@bolt$boltCarbon)
    cat('\n      Sliver carbon:', object@bolt$sectCarbon)
    
    cat('\n\nchainSawIZ...')
    cat('\n  use \"summary\" on the circularPlot slot for sample plot details')
    #summary(object@circularPlot)

    cat('\n')
    
    return(invisible())
}   #summary for 'chainSawIZ'
) #setMethod






#================================================================================
#  method for class "sausageIZ"...
#
setMethod('summary',
          signature(object = 'sausageIZ'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to 'InclusionZone' & 'downLogIZ' methods for 'sausageIZ'...
#------------------------------------------------------------------------------
    callNextMethod()
    
    cat('\nsausageIZ...')
    cat('\n  Spatial ID:', object@perimeter@polygons$pgsSausage@ID)
    if(object@units == .StemEnv$msrUnits$metric) {
      cat('\n  radius = ', object@radius, ' meters',sep='')
      cat('\n  area = ', object@area, ' square meters', sep='')
      cat(' (', format(object@area/.StemEnv$smpHectare, digits=4), ' hectares)')
    }
    else {
      cat('\n  radius = ', object@radius, ' feet',sep='')
      cat('\n  area = ', object@area, ' square feet', sep='')
      cat(' (', format(object@area/.StemEnv$sfpAcre, digits=4), ' acres)', sep='')
    }
    cat('\n  Number of perimeter points:', dim(object@sausage)[1], '(closed polygon)')


    cat('\n')
    
    return(invisible())
}   #summary for 'sausageIZ'
) #setMethod






#================================================================================
#  method for class "pointRelascopeIZ"...
#
setMethod('summary',
          signature(object = 'pointRelascopeIZ'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to 'InclusionZone' & 'downLogIZ' methods for 'pointRelascopeIZ'...
#------------------------------------------------------------------------------
    callNextMethod()
    
    cat('\npointRelascopeIZ...')
    cat('\n  Spatial ID:', object@perimeter@polygons$pgsPRS@ID)
    if(object@units == .StemEnv$msrUnits$metric) {
      cat('\n  dual circle radius = ', object@radius, ' meters',sep='')
      cat('\n  area = ', object@area, ' square meters', sep='')
      cat(' (', format(object@area/.StemEnv$smpHectare, digits=4), ' hectares)')
    }
    else {
      cat('\n  dual circle radius = ', object@radius, ' feet',sep='')
      cat('\n  area = ', object@area, ' square feet', sep='')
      cat(' (', format(object@area/.StemEnv$sfpAcre, digits=4), ' acres)', sep='')
    }
    cat('\n  Number of perimeter points:', dim(object@dualCircle)[1], '(closed polygon)')


    cat('\n')
    
    return(invisible())
}   #summary for 'pointRelascopeIZ'
) #setMethod




#================================================================================
#  method for class "perpendicularDistanceIZ"...
#
setMethod('summary',
          signature(object = 'perpendicularDistanceIZ'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to 'InclusionZone' & 'downLogIZ' methods for 'perpendicularDistanceIZ'...
#------------------------------------------------------------------------------
    callNextMethod()
    
  #  if(is(object, 'omnibusPDSIZ') || is(object, 'omnibusDLPDSIZ'))
  #    cat('\nomnibusPDSIZ...')
  #  else                                              
  #    cat('\nperpendicularDistanceIZ...')
    cat('\n',class(object),'...',sep='')
    cat('\n  PDS type:', object@pdsType)
    cat('\n  Spatial ID:', object@perimeter@polygons$pgs@ID)
    if(object@units == .StemEnv$msrUnits$metric) {
      cat('\n  area = ', object@area, ' square meters', sep='')
      cat(' (', format(object@area/.StemEnv$smpHectare, digits=4), ' hectares)',sep='')
    }
    else {
      cat('\n  area = ', object@area, ' square feet', sep='')
      cat(' (', format(object@area/.StemEnv$sfpAcre, digits=4), ' acres)', sep='')
    }
    cat('\n  Number of perimeter points:', dim(object@izPerim)[1], '(closed polygon)')


    cat('\n')
    
    return(invisible())
}   #summary for 'perpendicularDistanceIZ'
) #setMethod





#================================================================================
#  method for class "distanceLimitedPDSIZ"...
#
setMethod('summary',
          signature(object = 'distanceLimitedPDSIZ'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to 'InclusionZone' & 'downLogIZ' methods for 'distanceLimitedPDSIZ'...
#------------------------------------------------------------------------------
    callNextMethod()
    cat('  (The above summary is for the entire DLPDS region)\n')

    #if(is(object, 'omnibusDLPDSIZ'))
    #  cat('\nomnibusDLPDSIZ...')
    #else if (is(object, 'hybridDLPDSIZ'))
    #  cat('\ndistanceLimitedMCPDSIZ...')
    #else   
    #  cat('\ndistanceLimitedPDSIZ...')
    cat('\n',class(object),' details...',sep='')
    if(object@units == .StemEnv$msrUnits$metric) {
      cat('\n  distance limit = ', object@dls@distanceLimit, ' meters', sep='')
      cat('\n  limiting diameter = ', object@dlsDiameter, ' meters', sep='')
      cat(' (', format(object@dlsDiameter*.StemEnv$m2cm, digits=4), ' cm)',sep='')
    }
    else {
      cat('\n  distance limit = ', object@dls@distanceLimit, ' feet', sep='')
      cat('\n  limiting diameter = ', object@dlsDiameter, ' feet', sep='')
      cat(' (', format(object@dlsDiameter*.StemEnv$ft2in, digits=4), ' in)',sep='')
    }
    cat('\n  distance limited component available =',ifelse(is.null(object@dlsPart),FALSE, TRUE)) 
    cat('\n  PDS component available =',ifelse(is.null(object@pdsPart),FALSE, TRUE)) 
    cat('\n  Summaries of individual DLPDS components can be viewed seperately')

    cat('\n')
    
    return(invisible())
}   #summary for 'distanceLimitedPDSIZ'
) #setMethod





#================================================================================
#  method for class "distanceLimitedIZ"...
#
setMethod('summary',
          signature(object = 'distanceLimitedIZ'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to 'InclusionZone' & 'downLogIZ' methods for 'distanceLimitedIZ'...
#------------------------------------------------------------------------------
    callNextMethod()
    cat('  (Note: NAs signify location-dependent attributes.)\n')
    
    cat('\n',class(object),'...',sep='')
    cat('\n  Spatial ID:', object@perimeter@polygons$pgs@ID)
    if(object@units == .StemEnv$msrUnits$metric) {
      cat('\n  distance limit = ', object@dls@distanceLimit, ' meters', sep='')
      cat('\n  area = ', object@area, ' square meters', sep='')
      cat(' (', format(object@area/.StemEnv$smpHectare, digits=4), ' hectares)',sep='')
    }
    else {
      cat('\n  distance limit = ', object@dls@distanceLimit, ' feet', sep='')
      cat('\n  area = ', object@area, ' square feet', sep='')
      cat(' (', format(object@area/.StemEnv$sfpAcre, digits=4), ' acres)', sep='')
    }
    cat('\n  Number of perimeter points:', dim(object@izPerim)[1], '(closed polygon)')

    cat('\n')
    
    return(invisible())
}   #summary for 'distanceLimitedIZ'
) #setMethod




#---------------------------------------------------------------------------
#
#   For "standingTreeIZ" related methods
#
#   1. "standingTreeIZ"
#   2. "circularPlotIZ"
#
#Author...									Date: 1-Dec-2011
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
#  method for class "standingTreeIZ"...
#
setMethod('summary',
          signature(object = 'standingTreeIZ'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to 'InclusionZone' method for 'standingTreeIZ'...
#------------------------------------------------------------------------------
    callNextMethod()

    cat('\nstandingTree component estimates...')
    cat('\n  Spatial ID:', object@standingTree@spTree@polygons$pgsTree@ID)
    if(object@units == .StemEnv$msrUnits$metric) {
      cat('\n  Number of trees:', object@puaEstimates$Density, 'per hectare')
      cat('\n  Basal area:', object@puaEstimates$basalArea, 'square meters per hectare')
      cat('\n  Volume:', object@puaEstimates$volume, 'cubic meters per hectare')
      cat('\n  Surface area:', object@puaEstimates$surfaceArea, 'square meters per hectare')
      cat('\n  Biomass (stem):', object@puaEstimates$biomass, 'per hectare')
      cat('\n  Carbon content:', object@puaEstimates$carbon, 'per hectare')
    }
    else {
      cat('\n  Number of trees:', object@puaEstimates$Density, 'per acre')
      cat('\n  Basal area:', object@puaEstimates$basalArea, 'square feet per acre')
      cat('\n  Volume:', object@puaEstimates$volume, 'cubic feet per acre')
      cat('\n  Surface area:', object@puaEstimates$surfaceArea, 'square feet per acre')
      cat('\n  Biomass (stem):', object@puaEstimates$biomass, 'per acre')
      cat('\n  Carbon content:', object@puaEstimates$carbon, 'per acre')
    }

    cat('\n')
    
    return(invisible())
}   #summary for 'standingTreeIZ'
) #setMethod


#================================================================================
#  method for class "circularPlotIZ"...
#
setMethod('summary',
          signature(object = 'circularPlotIZ'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to 'InclusionZone' & 'standingTreeIZ' methods for 'circularPlotIZ'...
#------------------------------------------------------------------------------
    callNextMethod()
    
    cat('\ncircularPlotIZ...')
    cat('\n  use \"summary\" on the circularPlot slot for details')
    #summary(object@circularPlot)

    cat('\n')
    
    return(invisible())
}   #summary for 'circularPlotIZ'
) #setMethod




#================================================================================
#  method for class "horizontalPointIZ"...
#
setMethod('summary',
          signature(object = 'horizontalPointIZ'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to 'circularPlotIZ' method for 'horizontalPointIZ'...
#------------------------------------------------------------------------------
    callNextMethod()
    
    cat('\nhorizontalPointIZ...')
    cat('\n  use \"summary\" on the angleGauge slot for details')
    #summary(object@circularPlot)

    cat('\n')
    
    return(invisible())
}   #summary for 'horizontalPointIZ'
) #setMethod
