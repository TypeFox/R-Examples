#---------------------------------------------------------------------------
#
#   Summary methods for "monte" and associated classes...
#
#   1. "montePop"
#   2. "monteSample"
#   3. "monte"
#
#Author...									Date: 16-Feb-2012
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
#  1. method for class montePop...
#
setMethod('summary',
          signature(object = 'montePop'),
function(object,
         succinct = FALSE,
         ...
        )
{
#------------------------------------------------------------------------------
#   just a simple summary...
#------------------------------------------------------------------------------
    if(!succinct) {
      cat('\nObject of class:', class(object))
      if(nchar(object@description) > 0) {
        .StemEnv$underLine(60)
        if(!is.na(object@description))
          cat(object@description, fill=60)
        .StemEnv$underLine(60, prologue='')
      }
    }
 
    cat('\nPopulation...')
    cat('\n  Mean =', object@mean)
    cat('\n  Variance =', object@var)
    cat('\n  Standard Deviation =', object@stDev)
    cat('\n  Total =', object@total)
    cat('\n  Size (N) =', object@N)
    cat('\n  Zero-truncated =', object@zeroTruncated)

#   these are optionally defined in the object...    
    if(!all(is.na(object@n))) {
      cat('\n  Sample sizes (n) = ')
      cat(object@n, sep=', ')
    }
    if(!all(is.na(object@fpc))) {
      cat('\n  Finite population corrections = ')
      cat(format(object@fpc, digits=4), sep=', ')
    }
    if(!all(is.na(object@varMean))) {
      cat('\n  Variance of the mean = ')
      cat(object@varMean, sep=', ')
    }
    if(!all(is.na(object@stErr))) {
      cat('\n  Standard error of the mean = ')
      cat(object@stErr, sep=', ')
    }
    

    cat('\n')
    return(invisible())
}   #summary for 'montePop'
) #setMethod


#================================================================================
#  2. method for class montePopn...
#
#setMethod('summary',
#          signature(object = 'montePopn'),
#function(object,
         #succinct = FALSE,
#         ...
#        )
#{
#------------------------------------------------------------------------------
#   just a simple summary...
#------------------------------------------------------------------------------
#    callNextMethod(object, ...)

#    cat('  Sample sizes (n) = '); cat(object@n, sep=', ')
#    cat('\n  Finite population corrections = '); cat(object@fpc, sep=', ')
#    cat('\n  Variance of the mean = '); cat(object@varMean, sep=', ')
#    cat('\n  Standard error of the mean = '); cat(object@stErr, sep=', ')
    
#    cat('\n')
#    return(invisible())    
#}   #summary for 'montePopn'
#) #setMethod









  


#================================================================================
#  3. method for class "monteSample" and subclasses...
#
setMethod('summary',
          signature(object = 'monteSample'),
function(object,
         succinct = FALSE,
         ...
        )
{
#------------------------------------------------------------------------------
#   just a simple summary...
#------------------------------------------------------------------------------
    if(!succinct)
      cat('\nObject of class:', class(object))
    cat('\nNumber of Monte Carlo samples =', object@mcSamples)
    cat('\nSample sizes: n = '); cat(object@n, sep=', ')

    cat('\nSample summary statistics (mean values)...\n')
    print(object@stats)
    
    confLevel = (1.0 - object@alpha)*100
    cat('\nPercentage of confidence intervals (',confLevel,'%) that caught the population mean...\n',sep='')
    print(object@caughtPct)

    return(invisible())
}   #summary for 'monteSample'
) #setMethod
    







#================================================================================
#  4. method for class "monteBSSample" and subclasses...
#
setMethod('summary',
          signature(object = 'monteBSSample'),
function(object,
         succinct = FALSE,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to the default summary...
#------------------------------------------------------------------------------
    cat('\nNumber of bootstrap samples =',object@R)
    callNextMethod(object, succinct=succinct, ...)
    if(any(object@degenerate > 0)) {
      cat('\nNumber of degenerate bootstrap samples realized...\n')
      print(object@degenerate)
    }

    return(invisible())
}   #summary for 'monteBSSample'
) #setMethod
    


  

#================================================================================
#  5. method for class monte...
#
setMethod('summary',
          signature(object = 'monte'),
function(object,
         succinct = FALSE,
         ...
        )
{
#------------------------------------------------------------------------------
#   just a simple summary...
#------------------------------------------------------------------------------
#
    if(!succinct) {
      cat('\nObject of class:', class(object))
      if(nchar(object@description) > 0) {
        .StemEnv$underLine(60)
        if(!is.na(object@description))
          cat(object@description, fill=60)
        .StemEnv$underLine(60, prologue='')
      }
    }
    
    cat('Estimate attribute =', object@estimate)
    cat('\n')
    summary(object@pop, succinct = succinct, ...)

    cat('\nNormal theory results...')
    if(!is.null(object@NTsamples)) 
      summary(object@NTsamples, succinct = succinct, ...)
    else
      cat('\n  No normal theory information available.\n')


    
    cat('\nBootstrap results...')
    if(!is.null(object@BSsamples)) {
      summary(object@BSsamples, succinct = succinct, ...)
    }
    else
      cat('\n  No bootstrap information available.\n')
    

    cat('\n')
    return(invisible())
}   #summary for 'monte'
) #setMethod
    
