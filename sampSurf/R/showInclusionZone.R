#---------------------------------------------------------------------------
#
#   Methods for generic show() for class...
#     (1) InclusionZone and subclasses
#     (2) izContainer and subclass
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
#  method for base class InclusionZone & subclasses...
#
#  this will handle any subclass
#
setMethod('show',
          signature(object = 'InclusionZone'),
function(object)
{
    return(summary(object))
}   #show for 'InclusionZone'
) #setMethod





#================================================================================
#  method for base class izContainer...
#
#  this will handle any subclass
#
setMethod('show',
          signature(object = 'izContainer'),
function(object)
{
    return(summary(object))
}   #show for 'izContainer'
) #setMethod
