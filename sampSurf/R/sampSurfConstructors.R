#---------------------------------------------------------------------------
#
#   This file holds the S4 definition for the constructor methods of the
#   sampSurf class...
#
#   The methods include signatures for...
#     Signature: "object", "tract"
#     1. "izContainer", "Tract": Takes a collection of inclusion
#        zones already generated for the tract argument.
#     2. "numeric", "Tract": This will allow generating a sampling surface
#        from scratch. The object is the number of logs/trees, you can specify any
#        argument for generating logs/trees that will be passed on to downLogs or
#        standingTrees in the "..." argument list, similar for other methods.
#
#   Revamped to handle standing tree methods 5-Dec-2011, JHG.
#
#Author...									Date: 1-Oct-2010
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
setGeneric('sampSurf',  
           function(object, tract, ...) standardGeneric('sampSurf'),
             signature = c('object', 'tract')
            )




          
#================================================================================
#
# 1. Takes a collection of Stem inclusion zones and a "Tract" object...
#
setMethod('sampSurf',
          signature(object = 'izContainer', tract='Tract'), 
function(object, 
         tract,
         estimate = unlist(c(.StemEnv$puaEstimates, .StemEnv$ppEstimates)),
         description = 'sampling surface object',
         runQuiet = FALSE,
         ...
        )
{
#---------------------------------------------------------------------------
#
#   This is the main routine for calculating stats on the surface etc, any
#   other constructors should eventually call this one.
#
#   Arguments...
#
#---------------------------------------------------------------------------
#
#   it doesn't make sense to use the point-based chainSawIZ, use fullChainSawIZ instead...
#
    if(is(object@iZones[[1]], 'chainSawIZ'))                 #this request does not make sense
      stop('You must use \"fullChainSawIZ\" for the chainSaw method!')

#
#   check to see if we are using a boundary slopover correction method...
#
    bndCorrect = is(tract, 'mirageTract')   #|| is(tract, 'walkthroughTract')

#
#   throw a warning if any of the inclusion zones land outside the tract boundary...
#
    if(!bndCorrect) {  #*jhg*
      bb.iz = bbox(object)
      bb.tr = bbox(tract)
      if(any(bb.iz[,'min']<bb.tr[,'min']) || any(bb.iz[,'max']>bb.tr[,'max']))
        warning('Some object inclusion zones lie outside the tract--this will impart a bias!!')
    }

#
#   let's see what we are dealing with here...
#
    estimate = match.arg(estimate)
    if(is(object, 'downLogIZs')) {
      isLogs = TRUE
      stemName = 'log'
      if(!estimate %in% c(.StemEnv$validEstimates$downLogs, .StemEnv$ppEstimates))
        stop(paste(estimate,'is not a valid attribute for downLogs'))
    }
    else {
      isLogs = FALSE
      stemName = 'tree'
      if(!estimate %in% c(.StemEnv$validEstimates$standingTrees, .StemEnv$ppEstimates))
        stop(paste(estimate,'is not a valid attribute for standingTrees'))
    }

    
#
#   heap each inclusion zone in the collection...
#
    nStems = length(object@iZones)
    if(!runQuiet) {
      cat('\nNumber of ',stemName,'s in collection = ', nStems, sep='')
      cat('\nHeaping ',stemName,': ',sep='')
    }
    for(i in seq_len(nStems)) {
      if(!runQuiet)
        cat(i,',',sep='')
        if(bndCorrect) #*jhg*
          izg = izGridMirage(object@iZones[[i]], tract, ...)        #InclusionZoneGrid mirage
        else  
          izg = izGrid(object@iZones[[i]], tract, ...)              #InclusionZoneGrid
      tract = heapIZ(izg, tract, estimate = estimate, ...)        #heap it up
    }


#
#   get the true population attribute value for the collection...
#   note that it is possible for biomass and carbon to have logs with
#   no estimates (NA), so we must account for that since these quantities are
#   optional...
#
    if(isLogs)
      Stems = as(object, 'downLogs')
    else
      Stems = as(object, 'standingTrees')

    truth = switch(estimate,
                   volume =  Stems@stats['total','volume'],
                   Density = nStems,
                   Length =  Stems@stats['total','length'],
                   surfaceArea = Stems@stats['total','surfaceArea'],
                   coverageArea = Stems@stats['total','coverageArea'],
                   basalArea = Stems@stats['total','basalArea'],
                   biomass = Stems@stats['total','biomass'],
                   carbon = Stems@stats['total','carbon'],
                   NA
                  )

#
#   we must adjust the area of the tract in case it is different from one hectare
#   or one acre--this then adjusts the /acre or /hectare values to total values...
#
    unitArea = ifelse(object@units==.StemEnv$msrUnits$English, .StemEnv$sfpAcre, .StemEnv$smpHectare) 
    areaAdjust = area(tract)/unitArea                             #m^2/m^2 or ft^2/ft^2
    tract = setValues(tract, getValues(tract) * areaAdjust)
    
#
#   some surface stats using raster...
#
    surfStats = list( mean = cellStats(tract, mean), 
                      sum = cellStats(tract,sum),    
                      var = cellStats(tract, var), 
                      nc = ncell(tract),
                      max = maxValue(tract)       
                    )
    surfStats$stDev = sqrt(surfStats$var)
    surfStats$se = surfStats$stDev/sqrt(surfStats$nc)
    surfStats$bias = surfStats$mean - truth
    surfStats$popTotal = truth
             

#
#   create the object...
#
    ss = new('sampSurf',
             description = description,
             izContainer = object,
             tract = tract,
             estimate = estimate,
             surfStats = surfStats
            )

    if(!runQuiet)
      cat('\n')

    return(ss)
}   #sampSurf for "izContainer"
)   #setMethod






          
#================================================================================
#
# 2. Takes the number of stems and a "Tract" object, other arguments for, e.g.,
#    downLogs, can be passed via "..."
#
#
setMethod('sampSurf',
          signature(object = 'numeric', tract='Tract'),
function(object, 
         tract,
         iZone,
         estimate = unlist(c(.StemEnv$puaEstimates, .StemEnv$ppEstimates)),
         description = 'sampling surface object',
         runQuiet = FALSE,
         ...
        )
{
#---------------------------------------------------------------------------
#
#   Arguments...
#     iZone = a character string (or name) of a legal InclusionZone object,
#             the routine will check to make sure it is applicable to the
#             correct "Stem" subclass
#
#   a few checks...
#
    nStems = round(object)
    if(nStems<1)
      stop('You must specify a positive number of stems in "object"!')

#
#   make sure the inclusion zone constructor is a valid available type...
#
    if(!is.character(iZone))                        #extends takes a character name
      iZone = deparse(substitute(iZone))
    if(extends(iZone, 'downLogIZ')) {
      isLogs = TRUE                                 #proper English
      papa = getClass('downLogIZ')
    }
    else if(extends(iZone, 'standingTreeIZ')) {
      isLogs = FALSE
      papa = getClass('standingTreeIZ')
    }
    else                                  #catch non-InclusionZone subclass values...
      stop('Invalid inclusion zone constructor name supplied: iZone = ',iZone)
    #above test is not quite enough, iZone must actually be a subclass, not the parent itself...     
    validNames = names(papa@subclasses)
    if(is.na(match(iZone, validNames)))
      stop('Invalid inclusion zone constructor name supplied: iZone = ',iZone)
    if(iZone=='chainSawIZ')               #catch this error too before going to default constructor
      stop('You must use \"fullChainSawIZ\" zones for the chainSaw method!')
 

#
#   get the logs/trees and the inclusion zones...
#
    if(isLogs) {
      dlogs = downLogs(nStems, tract, ...)
      izs = downLogIZs(lapply(dlogs@logs, iZone, ...))
    }
    else {
      strees = standingTrees(nStems, tract, ...)
      izs = standingTreeIZs(lapply(strees@trees, iZone, ...))
    }

#
#   just apply the default constructor now...
#
    ss = sampSurf(izs, tract, estimate=estimate, 
                  description=description, runQuiet=runQuiet, ...)

    return(ss)
}   #sampSurf for "numeric"
)   #setMethod
    
