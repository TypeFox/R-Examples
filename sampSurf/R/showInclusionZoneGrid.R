#---------------------------------------------------------------------------
#
#   Methods for generic show() for class...
#     (1) InclusionZoneGrid and subclasses
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
#  method for base class InclusionZoneGrid & subclasses...
#
setMethod('show',
          signature(object = 'InclusionZoneGrid'),
function(object)
{
    return(summary(object))
}   #show for 'InclusionZoneGrid'
) #setMethod
