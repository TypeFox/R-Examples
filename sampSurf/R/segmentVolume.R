#---------------------------------------------------------------------------
#
#   A simple little generic that will allow us to get volumes for segments
#   within a "Stem". If the object has a solidType, then the default taper
#   equation is used, otherwise, a spline function is fitted to the taper
#   slot data of the object (solidType=NULL).
#
#   Different methods are available for argument signatures of...
#
#   1. 'list' -- the guts, called from the below wrapper methods, please use 
#                those instead to minimize potential errors
#   2. 'downLog' objects
#   3. 'standingTree' objects
#
#   Arguments (common to all)...
#     object = a Stem object (see 1-3 above)
#     segBnds = the segment height/length lower and upper bounds, if they
#               match the entire stem, then the stem volume is returned, if NULL
#               or NA, then they are set to the full stem
#
#   Returns...
#     the volume of the segment
#
#   **For an alternate way to calculate segment/bolt volumes for solidType=NULL
#     stems before splineVolume was written, see chainsawSliver code for bolt
#     volume.
#
#Author...									Date: 25-Apr-2013
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   generic definition...
#
if(!isGeneric("segmentVolume")) 
  setGeneric('segmentVolume',  
             function(object, ...) standardGeneric('segmentVolume'),
             signature = c('object')
            )





#================================================================================
#  1. method for class 'list' -- please do not use this method, it is the "guts"
#     for the other two wrapper methods which call it and are more convenient...
#
setMethod('segmentVolume',
          signature(object = 'list'),
          
function(object,
         segBnds = NULL,   #will return total stem volume
         ...
        )
{
#---------------------------------------------------------------------------
#
#   unpack the list first...
#
    height = object$height
    buttDiam = object$buttDiam
    topDiam = object$topDiam
    solidType = object$solidType
    stemVol = object$stemVol
    stemObject = object$stemObject

#
#   check to see if segment bounds are available--set to whole stem if not...
#
    if( any(is.null(segBnds)) || any(is.na(segBnds)) ) {
      segBnds[1] = 0
      segBnds[2] = height
    }

    if(segBnds[1] < 0 || segBnds[2] > height || segBnds[2] <= segBnds[1])
      stop('Illegal height bounds on stem for segment volume determination!')

#
#   if the segment length is the whole stem, we are done, return the stem volume...
#
    segLength = segBnds[2] - segBnds[1]
    if(identical(segLength, height))
      return(stemVol)

#
#   if taper comes from the default equation, use the corresponding volume equation;
#   otherwise use a spline, which may be better than Smalian's, especially if
#   the taper is not close to linear within sections...
#
    if(is.null(solidType))     
      segVol = .StemEnv$splineVolume(stemObject@taper, segBnds[1], segBnds[2], is(stemObject,'downLog'),
                                     stemObject@units)
    else                                                      
      segVol = .StemEnv$wbVolume(buttDiam, topDiam, height, solidType, segBnds[2]) -
               .StemEnv$wbVolume(buttDiam, topDiam, height, solidType, segBnds[1])
    
    return(segVol)
}   #segmentVolume
)   #setMethod


      

       
#================================================================================
#  2. method for class 'downLog'...
#
setMethod('segmentVolume',
          signature(object = 'downLog'),
          
function(object,
         segBnds = c(low=0,  up=object@logLen),
         ...
        )
{
#---------------------------------------------------------------------------
#
#   just set the list up, and call the version for that signature...
#
    stem = list()
    stem$height = object@logLen
    stem$buttDiam = object@buttDiam
    stem$topDiam = object@topDiam
    stem$solidType = object@solidType
    stem$stemVol = object@logVol
    stem$stemObject = object 

    return( segmentVolume(stem, segBnds, ...) )
}   #segmentVolume for "downLog"
)   #setMethod

      

       
#================================================================================
#  3. method for class 'standingTree'...
#
setMethod('segmentVolume',
          signature(object = 'standingTree'),
          
function(object,
         segBnds = c(low=0,  up=object@height),
         ...
        )
{
#---------------------------------------------------------------------------
#
#   just set the list up, and call the version for that signature...
#
    stem = list()
    stem$height = object@height
    stem$buttDiam = object@buttDiam
    stem$topDiam = object@topDiam
    stem$solidType = object@solidType
    stem$stemVol = object@treeVol
    stem$stemObject = object 

    return( segmentVolume(stem, segBnds, ...) )
}   #segmentVolume for "standingTree"
)   #setMethod
          
