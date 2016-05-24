#---------------------------------------------------------------------------
#
#   This class setup is for container objects for the objects that are
#   subclasses of the "Stem" virtual class.
#
#   1. StemContainer -- virtual
#   2. downLogs -- for objects of class "downLog"
#   3. standingTrees -- for objects of class "StandingTree"
#
#Author...									Date: 25-Oct-2011
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
#  1. The StemContainer virtual class is a container class for any number of subclass objects...
#
setClass('StemContainer',
    representation(units = 'character',              #English or metric units
                   bbox = 'matrix',                  #the overall containing bbox matrix limits
                   stats = 'data.frame',             #summary of volume, etc. of Stems in collection
                   description = 'character'         #a short description of the object
                  ),
    prototype = list(units = .StemEnv$msrUnits$metric,
                     bbox = matrix(rep(0,4), nrow=2, dimnames=list(c('x','y'), c('min','max'))),
                     description = ''
                    ),
    contains = 'VIRTUAL',
    validity = function(object) {
                 if(!(object@units %in% .StemEnv$msrUnits))
                   return('units of measure must be "English" or "metric"')

#                check on bbox matrix format...                 
                 if(!class(object@bbox) == 'matrix')
                   return('bbox slot must be a 2x2 matrix')
                 bboxNames = match(rownames(object@bbox), c('x','y'))
                 if(any(is.na(bboxNames)))
                   return('slot bbox rownames must be "x", "y"!')
                 bboxNames = match(colnames(object@bbox), c('min','max'))
                 if(any(is.na(bboxNames)))
                   return('slot bbox colnames must be "min", "max"!')
                                                       
                 return(TRUE)
               } #validity check
) #class StemContainer








#=================================================================================================
#=================================================================================================
#
#  2. The downLogs class (plural) is a container class for any number of "downLog" objects...
#
setClass('downLogs',
    representation(logs = 'list'                    #the log objects as a list
                   ##units = 'character',              #English or metric units
                   ##bbox = 'matrix',
                   ##stats = 'data.frame'              #summary of volume, etc. of logs in collection
                   #numLogs = 'numeric'#,              #number of log objects in logs
                   #spLogs = 'SpatialPolygons'        #for simplicity in plotting
                  ),
    prototype = list(logs = list(),                  #empty, zero-length list
                     bbox = matrix(rep(0,4), nrow=2, dimnames=list(c('x','y'), c('min','max')))
                    ),
    contains='StemContainer',                         #a subclass of the virtual 'StemContainer' class
    validity = function(object) {
                 if(!(object@units %in% .StemEnv$msrUnits))
                   return('units of measure must be "English" or "metric"')
                 
                 numLogs = length(object@logs)
                 if(numLogs < 1)
                   return('no logs in collection!')

                 for(i in seq_len(numLogs))
                   validObject(object@logs[[i]])

                 for(i in seq_len(numLogs))
                   if(object@units != object@logs[[i]]@units)
                     return('At least one log has the wrong units!')

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
                 units = object@logs[[1]]@units
                 for(i in seq_len(numLogs))
                   if(object@logs[[i]]@units != units)
                     return('You can not mix measurement units within a population of logs!')
                 
                 return(TRUE)
               } #validity check
) #class downLogs 




#=================================================================================================
#=================================================================================================
#
#  3. The standingTrees class (plural) is a container class for any number of "standingTree"
#     objects...
#
setClass('standingTrees',
    representation(trees = 'list'                     #the standingTree objects as a list
                   #numTrees = 'numeric'#,              #number of standingTree objects in trees
                  ),
    prototype = list(trees = list(),                  #empty, zero-length list
                     description = 'container of standingTree objects'
                    ),
    contains='StemContainer',                         #a subclass of the virtual 'StemContainer' class
    validity = function(object) {
                 numTrees = length(object@trees)
                 if(numTrees < 1)
                   return('no "standingTree" objects found in "trees" slot!')

                 for(i in seq_len(numTrees))
                   validObject(object@trees[[i]])

                 for(i in seq_len(numTrees))
                   if(object@units != object@trees[[i]]@units)
                     return('At least one "standingTree" object has the wrong units!')
                 
#                consistent class check...
                 class = class(object@trees[[1]])
                 for(i in seq_len(numTrees))
                   if(class(object@trees[[i]]) != class)   #could us is() for softer comparison w/ inheritance?
                     return('You can not mix "Stem" classes in the population!')
                     
#                consistent units check...
                 units = object@trees[[1]]@units
                 for(i in seq_len(numTrees))
                   if(object@trees[[i]]@units != units)
                     return('You can not mix measurement units within a population of trees!')
                 
                 return(TRUE)
               } #validity check
) #class standingTrees 
