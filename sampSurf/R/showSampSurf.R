#---------------------------------------------------------------------------
#
#   Methods for generic show() for class...
#     (1) sampSurf and subclasses
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
#  method for base class sampSurf...
#
setMethod('show',
          signature(object = 'sampSurf'),
function(object)
{
    cat('\nObject of class:', class(object))
    .StemEnv$underLine(60)
    if(!is.na(object@description))
      cat(object@description, fill=60)
    .StemEnv$underLine(60, prologue='')

    cat('\nInclusion zone objects:', class(object@izContainer@iZones[[1]]) )
    cat('\nEstimate:', object@estimate)
    cat('\nNumber of logs =', length(object@izContainer@iZones) )

    cat('\n')
    print(object@tract)
    return(invisible())
}   #show for 'sampSurf'
) #setMethod
