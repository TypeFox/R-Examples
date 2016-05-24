#---------------------------------------------------------------------------
#
#    This file holds the definitions of the container classes for objects of
#    class "MonteCarloSampling" or "antitheticSampling".
#
#    Monte Carlo Sampling container...
#
#    1.  "mcsContainer": The base class definition
#    2a. "mcsContainer" constructor generic definition
#    2b. "mcsContainer" constructor method itself
#    3.  "mcsContainer" show method
#    4.  "mcsContainer" summary method
#    5.  "mcsContainer" hist method
#    6.  "mcsContainer" plot method
#
#    Antithetic sampling container...
#
#    7.  "antitheticContainer" class definition
#    8a. "antitheticContainer" constructor generic definition
#    8b. "antitheticContainer" constructor method
#
#Author...									Date: 16-May-2013
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#



#=================================================================================================
#
#  1. the mcsContainer class structure for collections of objects...
#
setClass('mcsContainer',
    representation(mcsObjs = 'list',                   #list of some subclass of "MonteCarloSampling"
                   stats = 'matrix',                   #MC summary stats for each tree or log
                   description = 'character'           #hmmm?
                  ),
         
    prototype = list(mcsObjs = list(),                  #empty, zero-length list
                     description = ''
                    ),
         
    validity = function(object) {

                 numObjs = length(object@mcsObjs)
                 if(numObjs < 1)
                   return('no "InclusionZone" objects found in mcsObjs slot!')

                 for(i in seq_len(numObjs))
                   validObject(object@mcsObjs[[i]])

#                consistent class check...
                 class = class(object@mcsObjs[[1]])
                 for(i in seq_len(numObjs))
                   if(class(object@mcsObjs[[i]]) != class)  #could us is() for softer comparison w/ inheritance?
                     return('Please do not mix Monte Carlo classes in the collection!')
                                 
                 return(TRUE)
               } #validity check
) #class mcsContainer 




#================================================================================
#
#  2a. generic definitions...
#
if(!isGeneric("mcsContainer")) 
  setGeneric('mcsContainer',  
             function(object, ...) standardGeneric('mcsContainer'),
             signature = c('object')
            )

          
#================================================================================
#  2b. base method for a previously made collection of mcs objects in a list...
#
setMethod('mcsContainer',
          signature(object = 'list'),
function(object,
         description = 'Monte Carlo Sampling container object',
         ...
        )
{
#------------------------------------------------------------------------------
#
#
    numObjs = length(object)
    if(numObjs < 1)
      stop('Error: there must be at least one object in the list for a collection!')


#
#   stem-based summary stats...
#
    volEst = sapply(object, function(x) x@volEst)
    volVar = sapply(object, function(x) x@volVar)
    ci.lo = sapply(object, function(x) x@ci.lo)
    ci.up = sapply(object, function(x) x@ci.up)
    trueVol = sapply(object, function(x) x@trueVol)
    relErrPct = sapply(object, function(x) x@relErrPct)
    
    stats = rbind(trueVol, volEst, relErrPct, volVar, ci.lo, ci.up)

    mcobj = new('mcsContainer',
                mcsObjs = object,
                stats = stats,
                description = description
               )

    return(mcobj)
    
}   #mcsContainer method for list
)   #setMethod
    



#================================================================================
#  3. showmethod for base class mcsContainer...
#
#  this will handle any subclass
#
setMethod('show',
          signature(object = 'mcsContainer'),
function(object)
{
    return(summary(object))
}   #show for 'mcsContainer'
) #setMethod




#================================================================================
#  4. summary method for base class mcsContainer & antitheticContainer...
#
#  Note that I cheated here and tested for class 'antitheticContainer' to
#  invoke some differences rather than inherit from this for that subclass...
#
setMethod('summary',
          signature(object = 'mcsContainer'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   just a simple summary of stats from collection...
#------------------------------------------------------------------------------
    .StemEnv$underLine(60)
    cat(object@description, fill=60)
    .StemEnv$underLine(60, prologue='')
    
    numObjs = length(object@mcsObjs)
    cat('There are ', numObjs, ' ', class(object@mcsObjs[[1]]), ' objects in the collection',sep='')
    if(is(object, 'antitheticContainer')) 
      cat('\n--Antithetic sampling from original', class(object@mcsObjs[[1]]@mcsObj),'objects.')

    cat('\n\nSummary stats over all objects...\n')
    statSum = apply(object@stats, 1, summary)
    print(statSum)

    cat('\nProxy tabulation...')
    if(is(object, 'antitheticContainer')) 
      print(table(sapply(object@mcsObjs, function(x)x@mcsObj@proxy)))
    else 
      print(table(sapply(object@mcsObjs,function(x)x@proxy)))    

    cat('\n')
    return(invisible(statSum))
}   #summary for 'mcsContainer'
) 




#================================================================================
#  5. hist method for class "mcsContainer"...
#
setMethod('hist',
          signature(x = 'mcsContainer'),
function(x,
         stat = c('relErrPct', 'volVar'),
         xlab = stat,
         main = NA,
         col = 'gray90',
         ...
        )
{
#------------------------------------------------------------------------------
#   simple for now...
#------------------------------------------------------------------------------
    stat = match.arg(stat)

    if(stat == 'relErrPct') 
      vals = sapply(x@mcsObjs,function(z) z@relErrPct)
    else
      vals = sapply(x@mcsObjs,function(z) z@volVar)
    
    hg = hist(vals, main=main, xlab=xlab, col=col, ...)

    return(invisible(hg))

}    #hist for 'mcsContainer'
) #setMethod



#================================================================================
#   6. plot method for class "mcsContainer"...
#
setMethod('plot',
          signature(x = 'mcsContainer', y='missing'),
function(x,
         xlab = 'Estimated volume',
         ylab = 'True volume',
         showDiagonal = TRUE,
         ...                                      #for plot
        )
{
#------------------------------------------------------------------------------
#   just a plot of true versus estimated volumes...
#------------------------------------------------------------------------------
#
    xvals = sapply(x@mcsObjs, function(z) z@volEst)
    yvals = sapply(x@mcsObjs, function(z) z@trueVol)

    plot(xvals, yvals, xlab=xlab, ylab=ylab, ...)
    if(showDiagonal)
      abline(0,1, lty='dashed', col='gray50')
        
    return(invisible())

}   #plot for 'mcsContainer'
) #setMethod
    
    







#################################################################################
#
#                 Antithetic Sampling container section...
#
#################################################################################




#=================================================================================================
#
#  7. the antitheticContainer class structure for collections of objects...
#
setClass('antitheticContainer',

    contains = 'mcsContainer',      #a subclass 
         
    validity = function(object) {

#                consistent class check...
                 class = class(object@mcsObjs[[1]])
                 if(class != 'antitheticSampling')
                   return('antitheticContainer can only have objects of class \"antitheticSampling\" !')
                                 
                 return(TRUE)
               } #validity check
) #class antitheticContainer 




#================================================================================
#
#  8a. generic definitions...
#
if(!isGeneric("antitheticContainer")) 
  setGeneric('antitheticContainer',  
             function(object, ...) standardGeneric('antitheticContainer'),
             signature = c('object')
            )

          
#================================================================================
#  8b. method for a previously made collection of antithetic mcs objects in a list...
#
setMethod('antitheticContainer',
          signature(object = 'list'),
function(object,
         description = 'Antithetic Sampling container object',
         ...
        )
{
#------------------------------------------------------------------------------
#
#   all this does is call the mase mcsContainer method, then cast it to the
#   antitheticContainer class...
#
    mcobj = mcsContainer(object, description, ...)

    anti = as(mcobj, 'antitheticContainer')

    return(anti)
    
}   #antitheticContainer method for mcsContainer
)   #setMethod
