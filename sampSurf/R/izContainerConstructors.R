#---------------------------------------------------------------------------
#
#   This file holds the S4 definition for the constructors of izContainer
#   subclasses. This is used for a collection or population of
#   any "InclusionZone" subclass, it has different methods for downLogIZ and
#   standingTreeIZ related inclusion zone collections...
#
#   Constructors include...
#     1. izContainer: base method functionality
#     2. downLogIZs: for collections of "downLogIZ" subclass objects
#     3. standingTreeIZs: for collections of standingTreeIZ subclass objects
#
#   Note that in both cases 2 & 3 above, there are different methods for
#   different classes of objects.
#
#   This was revamped from the original "downLogIZs" code, which no
#   longer exists. The changes were to move all of the common functionality
#   to the izContainer method, and then just create the objects in the
#   sublass methods. 1-Dec-2011, JHG.
#
#   Note that because "object" is always a list in the signature (albeit 
#   containing different class objects for the two subclasses), we can not 
#   make the subclass constructors methods of the base generic. So they 
#   all must be generics.
#
#Author...									Date: 24-Aug-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   generic definitions...
#
if(!isGeneric("izContainer")) 
  setGeneric('izContainer',  
             function(object, ...) standardGeneric('izContainer'),
             signature = c('object')
            )

if(!isGeneric("downLogIZs")) 
  setGeneric('downLogIZs',  
             function(object, ...) standardGeneric('downLogIZs'),
             signature = c('object')
            )

if(!isGeneric("standingTreeIZs")) 
  setGeneric('standingTreeIZs',  
             function(object, ...) standardGeneric('standingTreeIZs'),
             signature = c('object')
            )



          
#================================================================================
#  1. base method for a previously made collection of IZs in a list...
#
setMethod('izContainer',
          signature(object = 'list'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#
#   this handles the calculation of the overall bbox for each of the
#   subclasses...
#
    numIZs = length(object)
    if(numIZs < 1)
      stop('error in "object": must be at least one inclusion zone in the list')


#
#   make an array of the perimeter bbox matrices for each IZ object, then determine their
#   new overall extent...
#
    bboxArray = array(dim=c(2,2,numIZs))
    for(i in seq_len(numIZs)) {
      if(!is(object[[i]], 'InclusionZone'))     #catch objects that may not be in the correct form
        stop('All list elements must be a subclass of "InclusionZone"!')
      ###bboxArray[,,i] = bbox(perimeter(object[[i]])) #changed to below 21-Mar-2013, SU method could chop logs
      bboxArray[,,i] = bbox(object[[i]])               #that were larger than the IZ and near the plot border
    }
    dimnames(bboxArray) = dimnames(bbox(object[[1]])) #page dim doesn't matter
    bbox = bboxSum(bboxArray)                         #extend the bboxes to overall

    
    return(bbox)
}   #izContainer method for list
)   #setMethod






#********************************************
#  for downLogIZs objects...
#********************************************


          
#================================================================================
#  2. method previously made collection of downLog IZs in a list...
#
setMethod('downLogIZs',
          signature(object = 'list'),
function(object,
         description = '',
         ...
        )
{
#------------------------------------------------------------------------------
#
#   get the overall bbox for the items in the list from the class method...
#
    bbox = izContainer(object, ...)
  
#
#   create the object...
#    
    dl.iz = new('downLogIZs', iZones = object, bbox = bbox,
               units = object[[1]]@units, description = description)
    
    return(dl.iz)
}   #downLogIZs method for list
)   #setMethod




#================================================================================
#  2.1 method for  a previously made collection of downLog IZs from
#      a collection of "downLogs" objects...
#
setMethod('downLogIZs',
          signature(object = 'downLogs'),
function(object,
         iZone,
         description = '',
         ...
        )
{
#------------------------------------------------------------------------------
#
#   make sure the inclusion zone constructor is a valid available type;
#   note that I don't think the extended test for a valid subclass that is
#   in the sampSurf constructor is necessary here as we are only dealing
#   with downLogs...
#
    if(!is.character(iZone))                           #lapply will take a character name
      iZone = deparse(substitute(iZone))
    if(!extends(iZone, 'downLogIZ'))
      stop(paste('The inclusion zone specification',iZone,'is not for downLog objects!'))

#
#   now apply the inclusion zone build to each log and then create the container...
#
    dl.l = lapply(object@logs, iZone, ... )
    dl.iz = downLogIZs(dl.l, description=description, ...)
    
    return(dl.iz)
}   #downLogIZs method for "downLogs" objects
)   #setMethod

          







#********************************************
#  for standingTreeIZs objects...
#********************************************




          
#================================================================================
#  3. method for  a previously made collection of standingTree IZs in a list...
#
setMethod('standingTreeIZs',
          signature(object = 'list'),
function(object,
         description = '',
         ...
        )
{
#------------------------------------------------------------------------------
#
#   get the overall bbox for the items in the list from the class method...
#
    bbox = izContainer(object, ...)
  
#
#   create the object...
#    
    st.iz = new('standingTreeIZs', iZones = object, bbox = bbox,
               units = object[[1]]@units, description = description)
    
    return(st.iz)
}   #standingTreeIZs method for list
)   #setMethod



#================================================================================
#  3.1 method for  a previously made collection of standingTree IZs from
#      a collection of "standingTrees" objects...
#
setMethod('standingTreeIZs',
          signature(object = 'standingTrees'),
function(object,
         iZone,
         description = '',
         ...
        )
{
#------------------------------------------------------------------------------
#
#   make sure the inclusion zone constructor is a valid available type;
#   note that I don't think the extended test for a valid subclass that is
#   in the sampSurf constructor is necessary here as we are only dealing
#   with standingTrees...
#
    if(!is.character(iZone))                           #lapply will take a character name
      iZone = deparse(substitute(iZone))
    if(!extends(iZone, 'standingTreeIZ'))
      stop(paste('The inclusion zone specification',iZone,'is not for standingTree objects!'))

#
#   now apply the inclusion zone build to each tree and then create the container...
#
    st.l = lapply(object@trees, iZone, ... )
    st.iz = standingTreeIZs(st.l, description=description, ...)
    
    return(st.iz)
}   #standingTreeIZs method for "standingTrees" objects
)   #setMethod
