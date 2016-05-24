#---------------------------------------------------------------------------
#
#   Methods for generic summary() for sampSurf class...
#
#   Note: in some cases such as volume in the Monte Carlo methods for example,
#         cells can be within inclusion zones and have small (near zero)
#         values. In determining which cells are background cells (zero-valued)
#         vs. surface cells, we must therefore use the "digits" argument to
#         the raster freq() function because it in turn uses round() for the
#         actual comparisons and we don't want small real values rounded to
#         zero to be counted as background--see below.
#
#   Extended for standingTreeIZs-based sampling surfaces 5-Dec-2011, JHG.
#
#   Returns...
#     summary information invisibly
#
#Author...									Date: 5-Oct-2010
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
#  method for class Tract & subclasses...
#
setMethod('summary',
          signature(object = 'sampSurf'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   just a simple summary of common items...
#------------------------------------------------------------------------------
    cat('\nObject of class:', class(object))
    .StemEnv$underLine(60)
    if(!is.na(object@description))
      cat(object@description, fill=60)
    .StemEnv$underLine(60, prologue='')

    if(object@tract@units == .StemEnv$msrUnits$metric) {
      unitLen = 'meters'
      unitVol = 'cubic meters'
      unitSA = 'square meters'
    }
    else {
      unitLen = 'feet'
      unitVol = 'cubic feet'
      unitSA = 'square feet'
    }

#
#   let's see what we are dealing with...
#
    if(class(object@izContainer) == 'downLogIZs') {
      isLogs = TRUE
      stemName = 'log'
    }
    else {
      isLogs = FALSE
      stemName = 'tree'
    }
    
    cat('\nInclusion zone objects:', class(object@izContainer@iZones[[1]]) )
    if(.hasSlot(object@izContainer@iZones[[1]], 'pdsType'))                  #for PDS PPS derivatives
      cat(' (with PP to: ',object@izContainer@iZones[[1]]@pdsType,')',sep='')
    if(.hasSlot(object@izContainer@iZones[[1]], 'antithetic'))               #for Monte Carlo Sampling derivatives
      if(object@izContainer@iZones[[1]]@antithetic)                          #distinguish antithetic variants
        cat(' (antithetic)')
    cat('\nMeasurement units =', object@tract@units)

#
#   make a collection of Stems and get its statistics (population truth)...
#
    numStems = length(object@izContainer@iZones)
    cat('\nNumber of ',stemName,'s = ', numStems, sep='')
    if(isLogs)
      Stems = as(object@izContainer, 'downLogs')
    else
      Stems = as(object@izContainer, 'standingTrees')
    
    
    cat('\nTrue',stemName,'volume =', Stems@stats['total','volume'],unitVol)
    if(isLogs)
      cat('\nTrue',stemName,'length =', Stems@stats['total','length'],unitLen)
    else
      cat('\nTrue',stemName,'basal area =', Stems@stats['total','basalArea'],unitSA)
    cat('\nTrue',stemName,'surface area =', Stems@stats['total','surfaceArea'],unitSA)
    if(isLogs)
      cat('\nTrue',stemName,'coverage area =', Stems@stats['total','coverageArea'],unitSA)
    cat('\nTrue',stemName,'biomass =', Stems@stats['total','biomass'])
    cat('\nTrue',stemName,'carbon =', Stems@stats['total','carbon'])
    cat('\n\nEstimate attribute:', object@estimate)


#
#   and the surface stats for comparison; send everything back in the summary...
#
    cat('\nSurface statistics...')
    cat('\n  mean =', object@surfStats$mean)
    cat('\n  bias =', object@surfStats$bias )

    truth = switch(object@estimate,
                   volume =  Stems@stats['total','volume'],
                   Density = numStems,
                   Length =  Stems@stats['total','length'],
                   surfaceArea = Stems@stats['total','surfaceArea'],
                   coverageArea = Stems@stats['total','coverageArea'],
                   basalArea = Stems@stats['total','basalArea'],
                   biomass = Stems@stats['total','biomass'],
                   carbon = Stems@stats['total','carbon'],
                   NA
                  )
    summaryNames = c('estimate','truth','mean','pctBias','sum','var','stDev','pctCV','max',
                     'gcTot','gcRes','gcBack','gcIZ')
    vals = vector('list', length(summaryNames))            #return summary values
    names(vals) = summaryNames
    vals$estimate = object@estimate
    vals$truth = truth
    vals$mean = object@surfStats$mean

    vals$pctBias = object@surfStats$bias/truth*100
    cat('\n  bias percent =', vals$pctBias)
    vals$sum = object@surfStats$sum
    cat('\n  sum =', object@surfStats$sum)
    vals$var = object@surfStats$var
    cat('\n  var =', object@surfStats$var)
    vals$stDev = object@surfStats$stDev
    cat('\n  st. dev. =', object@surfStats$stDev)
    vals$pctCV = 100*object@surfStats$stDev/object@surfStats$mean
    cat('\n  cv % =', vals$pctCV)
    vals$max = object@surfStats$max
    cat('\n  surface max =', object@surfStats$max)
    #cat('\n  st. error =', object@surfStats$se)
    vals$gcTot = object@surfStats$nc
    cat('\n  total # grid cells =', object@surfStats$nc)
    vals$gcRes = xres(object@tract)
    cat('\n  grid cell resolution (x & y) =', xres(object@tract), unitLen)
    ncellZero = freq(object@tract, 0, digits=15)               #zero cells, note freq rounds, use digits
    vals$gcBack = ncellZero
    cat('\n  # of background cells (zero) =', ncellZero)
    vals$gcIZ = object@surfStats$nc - ncellZero
    cat('\n  # of inclusion zone cells =', vals$gcIZ)
    cat('\n')


    
    #summary(object@tract)
    cat('\n')
    
    return(invisible(vals))
}   #summary for 'sampSurf'
) #setMethod

