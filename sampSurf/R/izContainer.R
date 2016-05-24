#---------------------------------------------------------------------------
#
#    This file holds the definitions of the inclusion zone container classes.
#    It contains the following hierarchy, which can be added to...
#
#    1. izContainer: A virtual class that is the superclass of the following...       
#    2. downLogIZs: a container class for multiple objects of any subclass
#                    of 'downLogIZ' (original 20-Aug-2010)
#    3. standingTreeIZs: for objects of subclass of 'standingTreeIZ' (Nov-2011)
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






#=================================================================================================
#=================================================================================================
#
#  1. the virtual izContainer class structure for inclusion zone objects...
#
setClass('izContainer',
    representation(iZones = 'list',                    #list of some subclass of "InclusionZone"
                   units = 'character',                #English or metric units
                   bbox = 'matrix',                    #overall bounding box
                   description = 'character'           #hmmm?
                  ),
    prototype = list(iZones = list(),                  #empty, zero-length list
                     units = .StemEnv$msrUnits$metric,
                     bbox = matrix(rep(0,4), nrow=2, dimnames=list(c('x','y'), c('min','max'))),
                     description = ''
                    ),
    validity = function(object) {
                 if(!(object@units %in% .StemEnv$msrUnits))
                   return('units of measure must be "English" or "metric"')

                 numIZs = length(object@iZones)
                 if(numIZs < 1)
                   return('no "InclusionZone" objects found in iZones slot!')

                 for(i in seq_len(numIZs))
                   validObject(object@iZones[[i]])

                 for(i in seq_len(numIZs))
                   if(object@units != object@iZones[[i]]@units)
                     return('At least one inclusion zone has the wrong units!')

#                check on bbox matrix format...                 
                 if(!class(object@bbox) == 'matrix')
                   return('bbox slot must be a 2x2 matrix')
                 bboxNames = match(rownames(object@bbox), c('x','y'))
                 if(any(is.na(bboxNames)))
                   return('slot bbox rownames must be "x", "y"!')
                 bboxNames = match(colnames(object@bbox), c('min','max'))
                 if(any(is.na(bboxNames)))
                   return('slot bbox colnames must be "min", "max"!')
                                      
#                consistent units check...
                 units = object@iZones[[1]]@units
                 for(i in seq_len(numIZs))
                   if(object@iZones[[i]]@units != units)
                     return('You can not mix measurement units within a population of inclusion zones!')

#                consistent class check...
                 class = class(object@iZones[[1]])
                 for(i in seq_len(numIZs))
                   if(class(object@iZones[[i]]) != class)  #could us is() for softer comparison w/ inheritance?
                     return('You can not mix inclusion classes in the population!')
                                 
                 return(TRUE)
               } #validity check
) #class izContainer 




#=================================================================================================
#=================================================================================================
#
#  2. the downLogIZs class (plural) is a container class for a number of "downLogIZ" objects...
#
setClass('downLogIZs',
    contains = 'izContainer',
    validity = function(object) {
#                we just need to check one izone, since the izContainer validity checks for all
#                of the same type...
                 #if(!extends(class, 'InclusionZone'))  #last test passed for all, so just check first
                 if(!is(object@iZones[[1]], 'downLogIZ'))  
                   return('Classes of objects in iZones must be a subclass of "downLogIZ"!')                 
                 
                 return(TRUE)
               } #validity check
) #class downLogIZs 



#=================================================================================================
#=================================================================================================
#
#  2. the standingTreeIZs class (plural) is a container class for a number of "standingTreeIZ" objects...
#
setClass('standingTreeIZs',
    contains = 'izContainer',
    validity = function(object) {
#                we just need to check one izone, since the izContainer validity checks for all
#                of the same type...
                 if(!is(object@iZones[[1]], 'standingTreeIZ'))  
                   return('Classes of objects in iZones must be a subclass of "standingTreeIZ"!')                 
                 
                 return(TRUE)
               } #validity check
) #class standingTreeIZs 
