#---------------------------------------------------------------------------
#
#   Methods for generic show() for class...
#     (1) Stem and subclasses (e.g., downLog, standingTree)
#     (2) StemContainer and subclasses (e.g., downLogs, standingTrees)
#
#   For now, these all just call the corresponding summary method.
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
#  method for base class Stem & subclasses...
#
setMethod('show',
          signature(object = 'Stem'),
function(object)
{
    return(summary(object))
}   #show for 'Stem'
) #setMethod


#================================================================================
#  method for data frames and class "StemContainer" & subclasses...
#
setMethod('show',
          signature(object = 'StemContainer'),
function(object)
{
    return(summary(object))
}   #show for 'StemContainer'
) #setMethod
