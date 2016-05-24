#---------------------------------------------------------------------------
#
#   Methods for generic show() for class...
#     (1) Tract and subclasses
#
#Author...									Date: 23-Sept-2010
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
#  method for base class Tract...
#
setMethod('show',
          signature(object = 'Tract'),
function(object)
{
    .StemEnv$underLine(60)
    if(!is.na(object@description))
      cat(object@description, fill=60)
    .StemEnv$underLine(60, prologue='')
    
    cat('Measurement units =', object@units)
    #area = nrow(object)*ncol(object)*xres(object)^2
    area = object@area
    if(object@units == .StemEnv$msrUnits$metric) {
      cat('\nArea in square meters = ', area,' (',area/.StemEnv$smpHectare,' hectares)',sep='')
    }
    else {
      cat('\nArea in square feet = ', area,' (',area/.StemEnv$sfpAcre,' acres)',sep='')
    }
    cat('\n\n')
    
    callNextMethod()
    return(invisible())
}   #show for 'Tract'
) #setMethod


#================================================================================
#  method for base class bufferedTract & subclasses...
#
setMethod('show',
          signature(object = 'bufferedTract'),
function(object)
{
    callNextMethod()

    bw = object@bufferRect['x',1] - bbox(object)['x',1]
    cat('Buffer width = ', bw)
    cat('\n\n')

    return(invisible())
}   #show for 'Tract'
) #setMethod
