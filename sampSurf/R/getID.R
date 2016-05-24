#---------------------------------------------------------------------------
#
#   This generic will simply return the ID associated with the class of object
#   as determined in the signature "object"
#
#   Current methods include...
#     1. "standingTree" objects
#     2. "downLog" objects
#     3. "standingTrees" objects
#     4. "downLogs" objects
#
#   We could add others in the future, in which case if there are multiple spatial
#   objects, the IDS could be returned in a data frame. We could add an "all"
#   argument to the methods to select more than the base ID desired.
#
#   Note that standingTrees have both spTree and spDBH objects, right now, this
#   generic only returns the spTree ID. The spDBH ID is the same, but with "dbh."
#   prepended to it.
#
#
#Author...									Date: 11-Apr-2013
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#   generic definition...
#
if(!isGeneric("getID")) 
  setGeneric('getID',  
             function(object, ...) standardGeneric('getID'),
             signature = c('object')
            )




          
#================================================================================
#  1. for "standingTree" class objects...
#
setMethod('getID',
          signature(object = 'standingTree'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#
    spID = slot( slot(object@spTree,"polygons")[[1]], "ID")
    return(spID)
}   #standingTree
)   #setMethod


          
#================================================================================
#  2. for "downLog" class objects...
#
setMethod('getID',
          signature(object = 'downLog'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#
    spID = slot( slot(object@spLog,"polygons")[[1]], "ID")
    return(spID)
}   #downLog
)   #setMethod

    

          
#================================================================================
#  3. for "standingTrees" class objects...
#
setMethod('getID',
          signature(object = 'standingTrees'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#
    spIDs = sapply(object@trees, getID)
    return(spIDs)
}   #standingTrees
)   #setMethod
    

          
#================================================================================
#  4. for "downLogs" class objects...
#
setMethod('getID',
          signature(object = 'downLogs'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#
    spIDs = sapply(object@logs, getID)
    return(spIDs)
}   #downLogs
)   #setMethod
