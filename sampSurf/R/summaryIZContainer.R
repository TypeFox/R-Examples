#---------------------------------------------------------------------------
#
#   summary methods for izContainer...
#
#   1. izContainer class
#   2. downLogIZs class
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
#  method for virtual class izContainer...
#
setMethod('summary',
          signature(object = 'izContainer'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   just a simple summary of common items from virtual class...
#------------------------------------------------------------------------------
    cat('\nContainer object of class:', class(object))
    .StemEnv$underLine(60)
    cat(object@description, fill=60)
    .StemEnv$underLine(60, prologue='')
    numIZs = length(object@iZones)
    cat('There are',numIZs,'inclusion zones in the population')
    cat('\nInclusion zones are of class:', class(object@iZones[[1]]))
    if(.hasSlot(object@iZones[[1]], 'antithetic'))                    #for Monte Carlo Sampling derivatives
      if(object@iZones[[1]]@antithetic)                               #distinguish antithetic variants
        cat(' (antithetic)')
    cat('\nUnits of measurement: ', object@units)
    
    izAreas = sapply(object@iZones, area)
    cat('\nSummary of inclusion zone areas')
    if(object@units=='English')
      cat(' in square feet...\n')
    else
      cat(' in square meters...\n')
    asum = summary(izAreas)
    nas = names(asum)
    asum = c(asum, var(izAreas), sd(izAreas))
    names(asum) = c(nas, 'var', 'SDev')
    print(asum)
 
    cat('\nEncapulating bounding box...\n')
    print(object@bbox)

    cat('\n')
    return(invisible())
}   #summary for 'izContainer'
) 


#not required right now, perhaps later...
#================================================================================
#  method for data frames and class "downLogIZs" (plural!)...
#
#setMethod('summary',
#          signature(object = 'downLogIZs'),
#function(object,
#         ...
#        )
#{
#------------------------------------------------------------------------------
#   just a simple summary of items in the "downLogIZs" object...
#------------------------------------------------------------------------------
#    callNextMethod()
   
#    cat('\n')
#    return(invisible())
#}   #summary for 'downLogIZs'
#) 
