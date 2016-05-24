#---------------------------------------------------------------------------
#
#   Methods for generic summary() for MonteCarloSampling class...
#     (1) 'MonteCarloSampling' virtual base class
#     (2) 'crudeMonteCarlo' and subclasses
#     (3) 'antitheticSampling'
#
#
#Author...									Date: 1-May-2013
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
#  method for virtual class MonteCarloSampling...
#
setMethod('summary',
          signature(object = 'MonteCarloSampling'),
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
    
    return(invisible())
}   #summary for 'MonteCarloSampling'
) #setMethod







#================================================================================
#  summary method for class 'crudeMonteCarlo' and subclasses...
#
setMethod('summary',
          signature(object = 'crudeMonteCarlo'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   just a simple summary of common items--added from virtual class...
#------------------------------------------------------------------------------
    callNextMethod()

    if(is(object@stem, 'downLog')) {
      objName = 'log'
      hgtName = 'length'
    }
    else {
      objName = 'tree'
      hgtName = 'height'
    }
    
    cat('Original Stem object class:', class(object@stem) )
    cat('\nProxy taper function:', object@proxy)
    taper = object@stem@taper
    height = taper[,hgtName]
    cat('\nFull',objName, hgtName, '=', max(taper[,hgtName]))
    cat('\nSegment ', hgtName, ' bounds = ', object@segBnds[1], ' to ', object@segBnds[2],sep='')
    cat('\nTrue volume =', object@trueVol)
    cat('\nVolume estimate =', object@volEst)
    cat('\nRelative error % =', object@relErrPct)
    cat('\nVariance estimate =', object@volVar)
    cat('\n', 1-object@alphaLevel,'% confidence Interval = ', format(object@ci.lo, digits=3),
        ' to ', format(object@ci.up, digits=4), sep='')
    n.s = length(object@hgt.s)
    cat('\nNumber of samples n =', n.s)

#
#   if you get this message, then modify the proxy you are using so that it does not
#   taper to zero at the tip...
#
    if(any(is.na(object@vol.s))) {                          #NaNs have been converted to NAs
      na.is = length(object@vol.s[is.na(object@vol.s)])
      cat('\n***>Warning,', na.is,
          'proxy zero diameters at tip produced NaNs for volume & have been removed: adjusted n =',
          n.s-na.is)
      if(object@relErrPct > 100)
        cat('\n            Remaining small diameters may have contributed to inflated volume estimates!')
      cat('\nConsider modifying your proxy so that it does not taper to zero!')
    }
    cat('\n')
        
    return(invisible())
}   #summary for 'crudeMonteCarlo'
) #setMethod








#================================================================================
#  method for class antitheticSampling...
#
setMethod('summary',
          signature(object = 'antitheticSampling'),
function(object,
         succinct = TRUE,
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

    if(!succinct) {
      cat(summary(object@mcsObj))
      cat(summary(object@mcsAnti))
    }
    else{
      mcsObj = object@mcsObj
      if(is(mcsObj@stem, 'downLog')) {
        objName = 'log'
        hgtName = 'length'
      }
      else {
        objName = 'tree'
        hgtName = 'height'
      }    
      cat('Monte Carlo Class objects:', class(mcsObj))
      cat('\nOriginal Stem object class:', class(mcsObj@stem) )
      cat('\nProxy taper function:', mcsObj@proxy)
      taper = mcsObj@stem@taper
      height = taper[,hgtName]
      cat('\nFull',objName, hgtName, '=', max(taper[,hgtName]))
      cat('\nSegment ',hgtName,' bounds = ', mcsObj@segBnds[1], ' to ', mcsObj@segBnds[2],sep='')
      cat('\nTrue volume =', mcsObj@trueVol)
    }
    
    cat('\nAntithetic estimates...')
    cat('\n  Volume estimate =', object@volEst)
    cat('\n  Relative error % =', object@relErrPct)
    cat('\n  Variance estimate =', object@volVar)
    cat('\n  ', 1-object@alphaLevel,'% confidence Interval = ', format(object@ci.lo, digits=4),
        ' to ', format(object@ci.up, digits=4), sep='')
    
    if(succinct) 
      cat('\n  Number of samples n =', mcsObj@n.s)
    
    cat('\n')
    
    return(invisible())
}   #summary for 'antitheticSampling'
) #setMethod

