boltDimensions = function(dlog, runQuiet = FALSE)
{
#---------------------------------------------------------------------------
#
#   This function will calculate segment/bolt dimensions like volume
#   and surface area for an object of class "downLog". It will do this
#   for each bolt in the taper slot's data frame object. A summary is printed
#   if desired.
#
#   Arguments...
#     dlog = a sampSurf "downLog" object
#     runQuiet = TRUE: no feedback; FALSE: print a summary
#
#   Returns...
#     a data frame with self-evident column names invisibly.
#
#Author...									Date: 15-Feb-2011
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
    require(sampSurf)
    
    if(!is(dlog, 'downLog'))
      stop('must pass a downLog argument!')

#
#   get taper info, etc...
#
    taper = dlog@taper
    nTaper = nrow(taper)
    buttDiam = dlog@buttDiam
    topDiam = dlog@topDiam
    logLen = dlog@logLen
    solidType = dlog@solidType

#
#   bolt diameters and lengths...
#
    nSegs = nTaper - 1
    boltDim = matrix(NA, nrow=nSegs, ncol=5)
    boltDim[, 1] = taper[1:nSegs, 'diameter']   #bottom diameter of bolt
    boltDim[, 2] = taper[2:nTaper, 'diameter']  #top diameter of bolt
    boltDim[, 3] = taper[1:nSegs, 'length']     #length to bottom of bolt
    boltDim[, 4] = taper[2:nTaper, 'length']    #length to top of bolt
    boltDim[, 5] = diff(taper$length)           #bolt length

#
#   bolt volumes, biomass and carbon...
#
    vol = matrix(NA, nrow=nSegs, ncol=1)
    biomass = vol
    carbon = vol
    if(!is.null(solidType))
      for(i in seq_len(nSegs)) 
        vol[i, 1] = .StemEnv$wbVolume(buttDiam, topDiam, logLen, solidType, boltDim[i,4]) -
                    .StemEnv$wbVolume(buttDiam, topDiam, logLen, solidType, boltDim[i,3])
    else 
      vol[, 1] = .StemEnv$SmalianVolume(taper)$boltVol
    if(!is.na(dlog@biomass))
      biomass[, 1] = vol[,1]*dlog@conversions['volumeToWeight']
    if(!is.na(dlog@carbon))
      carbon[, 1] = biomass[,1]*dlog@conversions['weightToCarbon']

#
#   bolt surface areas...
#
    sa = matrix(NA, nrow=nSegs, ncol=1)
    for(i in seq_len(nSegs)) {
      if(!is.null(solidType))
        sa[i, 1] = .StemEnv$wbSurfaceArea(buttDiam, topDiam, logLen, solidType, boltDim[i,3],
                                        boltDim[i,4] )
      else
        sa[i, 1] = .StemEnv$splineSurfaceArea(taper, boltDim[i,3], boltDim[i,4] )
    }

#
#   bolt coverage areas...
#
    ca = matrix(NA, nrow=nSegs, ncol=1)
    for(i in seq_len(nSegs)) {
      if(!is.null(solidType))
        ca[i, 1] = .StemEnv$wbCoverageArea(buttDiam, topDiam, logLen, solidType, boltDim[i,3],
                                        boltDim[i,4] )
      else
        ca[i, 1] = .StemEnv$splineCoverageArea(taper, boltDim[i,3], boltDim[i,4] )
    }
   
      
#
#   return in a data frame...
#
    df = data.frame(boltDim, vol, sa, ca, biomass, carbon)
    colnames(df) = c('botDiam', 'topDiam', 'botLen', 'topLen', 'boltLen',
                      .StemEnv$puaEstimates[c('volume', 'surfaceArea', 'coverageArea',
                                              'biomass', 'carbon')])


#
#   some results if interested...
#
    if(!runQuiet) {
      sums = colSums(df)
      if(dlog@units == 'metric') {
        lu = 'meters'
        cu = 'cubic meters'
        su = 'square meters'
      }
      else {
        lu = 'feet'
        cu = 'cubic feet'
        su = 'square feet'
      }
      if(is.null(solidType)) {
        fromArea = '(from spline fit)'
        fromVol = "(from Smalian's)"
      }
      else {
        fromArea = '(from taper equation)'
        fromVol = fromArea
      }
      cat('\nSummary of bolts in taper data frame...')
      .StemEnv$underLine(40,postfix='')
      cat('\n  Units =', dlog@units)
      cat('\n  Number of segments =', nSegs)
      cat('\n  Solid type =', ifelse(is.null(solidType), 'NULL', solidType))
      cat('\n  Total Length =', sums['boltLen'], lu)
      cat('\n  Total volume =', sums['volume'], cu, fromVol)
      cat('\n  Total biomass =', sums['biomass'])
      cat('\n  Total carbon =', sums['carbon'])
      cat('\n  Total surface area =', sums['surfaceArea'], su, fromArea)
      cat('\n  Total coverge area =', sums['coverageArea'], su, fromArea)
      cat('\n')
    }

    return(invisible(df))
}   #boltDimensions

