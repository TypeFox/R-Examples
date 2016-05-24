#---------------------------------------------------------------------------
#
#   Methods for generic summary() for class...
#     (1) StemContainer virtual class
#     (2) downLogs (plural) subclass
#     (3) standingTrees (plural) subclass
#
#Author...									Date: 26-Oct-2011
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
#  1. method for class StemContainer...
#
setMethod('summary',
          signature(object = 'StemContainer'),
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
    cat('Container class object...')
    cat('\n  Units of measurement: ', object@units)
    cat('\n\n  Encapulating bounding box...\n')
    print(object@bbox)

    cat('\n')
    return(invisible())
}   #summary for 'StemContainer'
) #setMethod







#================================================================================
#  2. method for class "downLogs" (plural!)...
#
setMethod('summary',
          signature(object = 'downLogs'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   just a simple summary of items in the "downLogs" object...
#------------------------------------------------------------------------------
    callNextMethod()
    numLogs = length(object@logs)
    cat('  There are',numLogs,'logs in the population')

#
#   totals over all logs...
#
    totVol = object@stats['total','volume']
    totSA = object@stats['total','surfaceArea']
    totCA = object@stats['total','coverageArea']
    totBms = object@stats['total','biomass']
    totC = object@stats['total','carbon']
    if(object@logs[[1]]@units == .StemEnv$msrUnits$metric) {
      cat('\n  Population log volume = ', totVol, 'cubic meters')
      cat('\n  Population log surface area = ', totSA, 'square meters')
      cat('\n  Population log coverage area = ', totCA, 'square meters')
    }
    else {
      cat('\n  Population log volume = ', totVol, 'cubic feet')
      cat('\n  Population log surface area = ', totSA, 'square feet')
      cat('\n  Population log coverage area = ', totCA, 'square feet')
    }
    if(!is.na(totBms))
      cat('\n  Population log biomass = ', totBms)
    if(!is.na(totC))
      cat('\n  Population log carbon = ', totC)


#
#   averages per log...
#
    avgVol = object@stats['mean','volume']
    avgSA = object@stats['mean','surfaceArea']
    avgCA = object@stats['mean','coverageArea']
    avgBms = object@stats['mean','biomass']
    avgC = object@stats['mean','carbon']
    avgLen = object@stats['mean','length']
    if(object@logs[[1]]@units == .StemEnv$msrUnits$metric) {
      cat('\n  Average volume/log = ', avgVol, 'cubic meters')
      cat('\n  Average surface area/log = ', avgSA, 'square meters')
      cat('\n  Average coverage area/log = ', avgCA, 'square meters')
      cat('\n  Average length/log = ', avgLen, 'meters')
    }
    else {
      cat('\n  Average volume/log = ', avgVol, 'cubic feet')
      cat('\n  Average surface area/log = ', avgSA, 'square feet')
      cat('\n  Average coverage area/log = ', avgCA, 'square feet')
      cat('\n  Average length/log = ', avgLen, 'feet')
    }
    if(!is.na(avgBms))
      cat('\n  Average biomass/log =', avgBms)
    if(!is.na(avgC))
      cat('\n  Average carbon/log =', avgC)
    
    cat('\n(**All statistics exclude NAs)')
    

    cat('\n')
    return(invisible())
}   #summary for 'downLogs'
) #setMethod








#================================================================================
#  3. method for class "standingTrees" (plural!)...
#
setMethod('summary',
          signature(object = 'standingTrees'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   just a simple summary of items in the "standingTrees" object...
#------------------------------------------------------------------------------
#
    callNextMethod()
    numTrees = length(object@trees)
    cat('  There are',numTrees,'trees in the population')

#
#   totals over all trees...
#
    totVol = object@stats['total','volume']
    totSA = object@stats['total','surfaceArea']
    totBms = object@stats['total','biomass']
    totC = object@stats['total','carbon']
    if(object@trees[[1]]@units == .StemEnv$msrUnits$metric) {
      cat('\n  Population tree volume = ', totVol, 'cubic meters')
      cat('\n  Population tree surface area = ', totSA, 'square meters')
    }
    else {
      cat('\n  Population tree volume = ', totVol, 'cubic feet')
      cat('\n  Population tree surface area = ', totSA, 'square feet')
    }
    if(!is.na(totBms))
      cat('\n  Population tree biomass = ', totBms)
    if(!is.na(totC))
      cat('\n  Population tree carbon = ', totC)

#
#   averages per tree...
#
    avgVol = object@stats['mean','volume']
    avgSA = object@stats['mean','surfaceArea']
    avgBms = object@stats['mean','biomass']
    avgC = object@stats['mean','carbon']
    avgHgt = object@stats['mean','height']
    if(object@trees[[1]]@units == .StemEnv$msrUnits$metric) {
      cat('\n  Average volume/tree = ', avgVol, 'cubic meters')
      cat('\n  Average surface area/tree = ', avgSA, 'square meters')
      cat('\n  Average height/tree = ', avgHgt, 'meters')
    }
    else {
      cat('\n  Average volume/tree = ', avgVol, 'cubic feet')
      cat('\n  Average surface area/tree = ', avgSA, 'square feet')
      cat('\n  Average height/tree = ', avgHgt, 'feet')
    }
    if(!is.na(avgBms))
      cat('\n  Average biomass/tree =', avgBms)
    if(!is.na(avgC))
      cat('\n  Average carbon/tree =', avgC)
    
    cat('\n(**All statistics exclude NAs)')

    
    cat('\n')
    return(invisible())
}   #summary for 'standingTrees'
) #setMethod
  
