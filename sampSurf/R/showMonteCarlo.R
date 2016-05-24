#################################################################################
#
#  show and summary routines...
#
#  Just two methods will handle all their subclasses...
#
#  1. "MonteCarloSampling"
#  2. "anititheticSampling" 
#
#
#Author...									Date: 9-May-2013
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
#  show method for class 'MonteCarloSampling'...
#
setMethod('show',
          signature(object = 'MonteCarloSampling'),
function(object)
{
    return(summary(object))
}   #show for 'MonteCarloSampling'
) #setMethod



#================================================================================
#  show method for class 'antitheticSampling'...
#
setMethod('show',
          signature(object = 'antitheticSampling'),
function(object)
{
    return(summary(object))
}   #show for 'antitheticSampling'
) #setMethod


