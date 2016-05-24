#---------------------------------------------------------------------------
#
#   Methods for generic summary() for Tract class...
#     (1) Tract class
#
#Author...									Date: 31-Aug-2010
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
          signature(object = 'Tract'),
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
      
    print(callNextMethod(object, ...) )
    cat('\n')
    
    return(invisible())
}   #summary for 'Tract'
) #setMethod



