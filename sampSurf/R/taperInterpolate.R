#---------------------------------------------------------------------------
#
#   A simple little generic that will allow us to interpolate either
#   diameters or lengths/heights from a Stem object at certain points.
#   This function uses either the built-in taper equation if the taper
#   slot was generated from this, or cubic spline interpolation otherwise
#   (in other words, the user has supplied the taper data from an unknown--to
#   sampSurf--taper equation).
#
#   Different methods are available for argument signatures of...
#
#   1. 'list' -- the guts, called from the below wrapper methods, please use 
#                those instead to minimize potential errors
#   2. 'downLog' objects
#   3. 'standingTree' objects
#
#   The use of 'height' in the taper slot and for tree height in standingTree
#   objects, versus 'length' in downLog objects, is the only difference that
#   requires different wrapper routines in 2 & 3 above.
#
#   Arguments...
#     object = a Stem object (see 1-3 above)
#     whichSense = 'diameter' to interpolated diameters, pts must be lengths/heights
#                  'length'/'height' to interpolate lengths, pts must be diameters
#     pts = see above, to interpolate diameters, these are lengths/heights;
#           to interpolate lengths/heights, these are diameters
#
#   Returns...
#     a vector of interpolated points.
#
#Author...									Date: 30-Jan-2013 (update)
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
if(!isGeneric("taperInterpolate")) 
  setGeneric('taperInterpolate',  
             function(object, ...) standardGeneric('taperInterpolate'),
             signature = c('object')
            )





#================================================================================
#  1. method for class 'list' -- please do not use this method, it is the "guts"
#     for the other two wrapper methods which call it and are more convenient...
#
setMethod('taperInterpolate',
          signature(object = 'list'),
          
function(object,
         whichSense = c('diameter', 'length', 'height'),
         pts = NULL,
         ...
        )
{
#---------------------------------------------------------------------------
#
#   deconstruct the list passed, otherwise everything is the same for
#   standingTrees or downLogs...
#
#---------------------------------------------------------------------------
#
    whichSense = match.arg(whichSense)

    length = object$length
    diameter = object$diameter
    solidType = object$solidType
    buttDiam = object$buttDiam
    topDiam = object$topDiam
    totLength = object$totLength
    

    if(is.null(pts) || is.na(pts) || length(pts) < 1)
      stop('You must provide some points for interpolation!')

#
#   interpolate either diameter or length as desired...
#
    if(whichSense == 'diameter') {                    #interpolate diameters at given lengths
      rlen = range(length)
      rpts = range(pts)
      if(rpts[1]<rlen[1] || rpts[2]>rlen[2])
        stop('Some interpolated points are off the stem!')

      if(!is.null(solidType)) 
        diams = .StemEnv$wbTaper(buttDiam, topDiam, totLength, nSegs=1,
                                 solidType=solidType, hgt=pts)$diameter 
      else
        diams = spline(length, diameter, xout=pts)$y
      return(diams)
    }
    else {                                           #interpolate lengths at given diameters
      rdiam = range(diameter)
      rpts = range(pts)
      if(rpts[1]<rdiam[1] || rpts[2]>rdiam[2])
        stop('Some interpolated points are off the stem!')

      if(!is.null(solidType))                        #use inverted .StemEnv$wbTaper equation
        lengths = totLength*(1 - ((pts - topDiam)/(buttDiam - topDiam))^(solidType/2)) 
      else
        lengths = spline(diameter, length, xout=pts)$y 
      return(lengths)
    }
      

    return(NULL)
}   #taperInterpolate for "list"
)   #setMethod


      

       
#================================================================================
#  2. method for class 'downLog'...
#
setMethod('taperInterpolate',
          signature(object = 'downLog'),
          
function(object,
         whichSense = c('diameter', 'length'),
         pts = NULL,
         ...
        )
{
#---------------------------------------------------------------------------
#
#   just set the list up, and call the version for that signature...
#
    whichSense = match.arg(whichSense)

    stem = list()
    stem$length = object@taper$length
    stem$diameter = object@taper$diameter
    stem$buttDiam = object@buttDiam
    stem$topDiam = object@topDiam
    stem$totLength = object@logLen
    stem$solidType = object@solidType

    return( taperInterpolate(stem, whichSense, pts, ...) )
}   #taperInterpolate for "downLog"
)   #setMethod
    

      

       
#================================================================================
#  3. method for class 'standingTree'...
#
setMethod('taperInterpolate',
          signature(object = 'standingTree'),
          
function(object,
         whichSense = c('diameter', 'height'),
         pts = NULL,
         ...
        )
{
#---------------------------------------------------------------------------
#
#   just set the list up, and call the version for that signature...
#
    whichSense = match.arg(whichSense)

    stem = list()
    stem$length = object@taper$height
    stem$diameter = object@taper$diameter
    stem$buttDiam = object@buttDiam
    stem$topDiam = object@topDiam
    stem$totLength = object@height
    stem$solidType = object@solidType

    return( taperInterpolate(stem, whichSense, pts, ...) )
}   #taperInterpolate for "standingTree"
)   #setMethod
    
  
