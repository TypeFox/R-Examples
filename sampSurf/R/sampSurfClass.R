#---------------------------------------------------------------------------
#
#   This file holds the S4 class definitions for the "sampSurf" class.
#
#   Note that we could not simply extend class Tract here, as we then would
#   have no mechanism to deal with subclasses of Tract (like bufferedTract)
#   in this class.
#
#Author...									Date: 1-Oct-2010
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
#  define the sampSurf class...
#
setClass('sampSurf',
#
#  slots for the class and its subclasses...
#
    representation(description = 'character',
                   izContainer = 'izContainer',      #collection of iz objects
                   tract = 'Tract',                  #the underlying "Tract" or subclass
                   estimate = 'character',           #the estimate for the tract layer
                   surfStats = 'list'                #sampling surface summary statistics
                  ),
    prototype = list(description = 'sampling surface object',  #some defaults for validity checking
                     estimate = '',
                     surfStats = list()
                    ),
    validity = function(object) {

                 if(!(object@estimate %in% c(.StemEnv$puaEstimates, .StemEnv$ppEstimates) ))
                   return('Invalid per unit area estimate in class sampSurf!')

#                check comparable units within tract and inclusion zones...
                 numIZs = length(object@izContainer@iZones)
                 for(i in seq_len(numIZs))
                   if(object@tract@units != object@izContainer@iZones[[i]]@units)
                     return('At least one inclusion zone object does not match tract object units!')
                   
                 return(TRUE)
               } #validity check
) #class sampSurf 
