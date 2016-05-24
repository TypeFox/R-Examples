#---------------------------------------------------------------------------
#
#   These methods will plot a histogram of objects of different classes as
#   e.g.,
#
#   1. "izContainer" -- inclusion zone areas
#   2. "sampSurf" -- the surface estimate at each grid cell
#   3. "downLogs" -- several different attibutes available
#   4. "standingTrees" -- same
#   5. "InclusionZoneGrid" -- like "sampSurf" (7-Feb-2013)
#
#Author...									Date: 10-May-2011
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
#  1. method for class izContainer...
#
setMethod('hist',
          signature(x = 'izContainer'),
function(x,
         xlab = 'Inclusion Zone Area',
         main = NA,
         col = 'gray90',
         ...
        )
{
#------------------------------------------------------------------------------
#   this just plots the histo of inclusion zone areas...
#------------------------------------------------------------------------------
#

    izAreas = sapply(x@iZones, area)
    hg = hist(izAreas, xlab=xlab, main=main, col=col, ...)

    return(invisible(hg))

}    #hist for 'izContainer'
) #setMethod
     


#================================================================================
#  2. method for class sampSurf...
#
setMethod('hist',
          signature(x = 'sampSurf'),
function(x,
         zeroTrunc = TRUE,              #exclude zeros?
         xlab = x@estimate,
         main = NA,
         col = 'gray90',
         ...
        )
{
#------------------------------------------------------------------------------
#   this just plots the sampling distribution histogram w/ or w/o zeros
#------------------------------------------------------------------------------
#
    values = getValues(x@tract)
    if(zeroTrunc) 
      values = values[values>0]
    hg = hist(values, xlab=xlab, main=main, col=col, ...)
    if(zeroTrunc)                                         #zero cells, note freq rounds(), use digits
      cat('\nHistogram is zero-truncated:',freq(x@tract, 0, digits=15),'zeros excluded.\n')  

    return(invisible(hg))

}    #hist for 'sampSurf'
) #setMethod
     


#================================================================================
#  3. method for class "downLogs" container...
#
setMethod('hist',
          signature(x = 'downLogs'),
function(x,
         logAttr = c('logLen','buttDiam','topDiam','logAngle','solidType','logVol',
                     'surfaceArea','coverageArea','biomass','carbon'),
         xlab = logAttr,
         main = NA,
         col = 'gray90',
         ...
        )
{
#------------------------------------------------------------------------------
#   this just plots one of the desired metrics for the logs...
#------------------------------------------------------------------------------
#
    logAttr = match.arg(logAttr)

    vals = sapply( x@logs, function(z) slot(z, logAttr) )
    #check for no carbon or biomass conversions in the collection...
    if(all(is.na(vals)))
      stop('All values for ',logAttr,' are NA!')

#
#   convert diameters to usual units...
#
    if(logAttr == 'buttDiam' || logAttr == 'topDiam')
       if(x@units == .StemEnv$msrUnits$metric) 
         vals = vals*.StemEnv$m2cm
       else
         vals = vals*.StemEnv$ft2in

#
#   and angles too...
#
    if(logAttr == 'logAngle')
      vals = sapply(vals, .StemEnv$rad2Deg)    #not vectorized at present V0.52

    
    hg = hist(vals, main=main, xlab=xlab, col=col, ...)

    return(invisible(hg))

}    #hist for 'downLogs'
) #setMethod
     
     


#================================================================================
#  4. method for class "standingTrees" container...
#
setMethod('hist',
          signature(x = 'standingTrees'),
function(x,
         treeAttr = c('height','buttDiam','dbh','topDiam','solidType','treeVol',
                      'surfaceArea','biomass','carbon'),
         xlab = treeAttr,
         main = NA,
         col = 'gray90',
         ...
        )
{
#------------------------------------------------------------------------------
#   this just plots one of the desired metrics for the trees...
#------------------------------------------------------------------------------
#
    treeAttr = match.arg(treeAttr)

    vals = sapply( x@trees, function(z) slot(z, treeAttr) )
    #check for no carbon or biomass conversions in the collection...
    if(all(is.na(vals)))
      stop('All values for ',treeAttr,' are NA!')

#
#   convert diameters to usual units...
#
    if(treeAttr == 'buttDiam' || treeAttr == 'topDiam' || treeAttr == 'dbh')
       if(x@units == .StemEnv$msrUnits$metric) 
         vals = vals*.StemEnv$m2cm
       else
         vals = vals*.StemEnv$ft2in
    
    hg = hist(vals, main=main, xlab=xlab, col=col, ...)

    return(invisible(hg))

}    #hist for 'standingTrees'
) #setMethod
     




#================================================================================
#  5. method for class "InclusionZoneGrid"...
#
setMethod('hist',
          signature(x = 'InclusionZoneGrid'),
function(x,
         zeroTrunc = TRUE,              #exclude zeros? (background cells)
         estimate = 'volume',           #see matching below
         main = NA,
         xlab = NA,
         col = 'gray90',
         ...
        )
{
#------------------------------------------------------------------------------
#   this just plots the sampling distribution histogram w/ or w/o zeros
#------------------------------------------------------------------------------
#
#   note that for the legal names below, we are not using, e.g.,
#   estimate = names(c(.StemEnv$puaEstimates, .StemEnv$ppEstimates)),
#   but are allowing whatever the creator of the particular sampling object put in
#   the data frame; thus, we are not using match.arg() here either...
#
    legalNames =  colnames(x@data)          
    edx = pmatch(estimate,legalNames)
    if(is.na(edx))
      stop('estimate must be one of:',paste(legalNames,collapse=','))
    else
      estimate = legalNames[edx]
    if(is.na(xlab))
      xlab = estimate
    
    values = x@data[,estimate]
    if(zeroTrunc)                     #background cells
      values = values[values>0]
    if(all(is.na(values)))
      stop('All values for ',estimate,' are NA')
    hg = hist(values, xlab=xlab, main=main, col=col, ...)
    if(zeroTrunc)                                         #zero cells, note freq rounds(), use digits
      cat('\nHistogram is zero-truncated:',length(x@data[,estimate])-length(values),'zeros excluded.\n')  

    return(invisible(hg))

}    #hist for 'InclusionZoneGrid'
) #setMethod
