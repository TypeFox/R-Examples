#---------------------------------------------------------------------------
#
#   Methods for generic show() for class...
#     (1) "montePop" and subclasses
#     (2) "monte" and subclasses
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
#  method for base class montePop & subclasses...
#
setMethod('show',
          signature(object = 'montePop'),
function(object)
{
    return(summary(object, succinct=TRUE))
}   #show for 'montePop'
) #setMethod



#================================================================================
#  method for base class montePop & subclasses...
#
setMethod('show',
          signature(object = 'monteSample'),
function(object)
{
    return(summary(object, succinct=TRUE))
}   #show for 'monteSample'
) #setMethod




#================================================================================
#  method for base class monte & subclasses...
#
setMethod('show',
          signature(object = 'monte'),
function(object)
{
    return(summary(object, succinct=TRUE))
}   #show for 'monte'
) #setMethod
