#---------------------------------------------------------------------------
#
#   This file holds the S4 class definitions for the class that defines the
#   minimal bounding grid at given resolution for an InclusionZone object.
#   These objects may be used later in generating a sampling surface by
#   "piling" or "heaping" them one on top of another within a "Tract"
#   object.
#
#
#Author...									Date: 17-Sept-2010
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
#  define the InclusionZoneGrid class...
#
setClass('InclusionZoneGrid',
#
#  slots for the class and its subclasses...
#
    representation(description = 'character',
                   iz = 'InclusionZone',             #iz object
                   grid = 'RasterLayer',             #for the grid
                   data = 'data.frame',              #pua estimates over the grid
                   bbox = 'matrix'                   #overall bounding box
                  ),
    prototype = list(description = 'gridded inclusion zone',  #some defaults for validity checking
                     bbox = matrix(rep(0,4), nrow=2, dimnames=list(c('x','y'), c('min','max'))),
                     data = data.frame(matrix(NA, nr=0, nc=length(.StemEnv$puaEstimates),
                            dimnames=list(character(0), names(.StemEnv$puaEstimates))) )
                    ),
    validity = function(object) {

                 #essentially the same checks as in bboxCheck()...
                 if(!nrow(object@bbox)==2 || !ncol(object@bbox)==2)
                   return('bbox slot must be a 2x2 matrix')
                 bboxNames = match(rownames(object@bbox), c('x','y'))
                 if(any(is.na(bboxNames)))
                   return('slot bbox rownames must be "x", "y"!')
                 bboxNames = match(colnames(object@bbox), c('min','max'))
                 if(any(is.na(bboxNames)))
                   return('slot bbox colnames must be "min", "max"!')
                 if(any( apply(object@bbox,1,function(x) if(x['min'] >= x['max']) TRUE else FALSE) ))
                   return('in slot bbox, "min" must be less than "max" for x and y!')

                 dfNames = match(colnames(object@data), c(names(.StemEnv$puaEstimates),
                                                          names(.StemEnv$ppEstimates)) )
                 if(any(is.na(dfNames)))
                   return('slot data colnames must contain all the per unit area estimate names')
                   
                 return(TRUE)
               } #validity check
) #class InclusionZoneGrid 
         



#=================================================================================================
#
#  define the InclusionZoneGrid class for the full chainsaw object where all possible
#  cuts are made within the sausage inclusion zone--a very specific class, but related
#  to the above; that is, for each grid cell within the inclusion zone, we apply the
#  chainSawIZ method and record the value of that cell...
#
setClass('csFullInclusionZoneGrid',
#
#  slots for the class; note that we need a list of "InclusionZoneGrid" objects, one for
#  each chainSaw estimate within the overall sausage inclusion zone...
#
    representation(chiz = 'list'                      #a list of InclusionZoneGrid objects
                  ),
    contains = 'InclusionZoneGrid',
    prototype = list(description = 'full chainsaw-sausage gridded inclusion zone',
                     chiz = list(),
                     bbox = matrix(rep(0,4), nrow=2, dimnames=list(c('x','y'), c('min','max'))),
                     data = data.frame(matrix(NA,
                                              nrow = 0,
                                              ncol = length(c(.StemEnv$puaEstimates,.StemEnv$ppEstimates)),
                                              dimnames = list(character(0),
                                                         names(c(.StemEnv$puaEstimates,.StemEnv$ppEstimates)))
                                             ) #matrix
                                      ) #df
                    ),
    sealed = TRUE,                           #no further changes or subclasses
    validity = function(object) {
                 #a check for "sausageIZ" would work below, but force it to be "fullChainSawIZ"...
                 if(!is(object@iz, 'fullChainSawIZ'))
                   return('The underlying inclusion zone must be of class "fullChainSawIZ".')

                 chizLen = length(object@chiz)
                 if(chizLen > 0) {
                   for(i in seq_len(chizLen)) {
                     if(isS4(object@chiz[[i]])) {               #can't check is.na on S4 objects!
                       if(!is(object@chiz[[i]], 'InclusionZoneGrid'))
                           return('All internal sausage grid cells must be InclusionZoneGrid objects!')
                       if(!is(object@chiz[[i]]@iz, 'chainSawIZ'))
                         return('Each internal sausage grid cell must be from a chainSawIZ object!')
                     }
                     else if(!is.na(object@chiz[[i]]))
                       return('External sausage grid cells must have value "NA".')
                   }
                 }
 
                   
                 return(TRUE)
               } #validity check
) #class csFullInclusionZoneGrid 
         
