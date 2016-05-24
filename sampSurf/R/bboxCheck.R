bboxCheck = function(bbox)
{
#---------------------------------------------------------------------------
#
#   A bbox is just a special matrix, unfortunately, the sp or raster
#   packages did not make it a class, so there is no check on validity
#   anywhere. This little function attempts to alleviate that inconsitency
#   without introducing it as a new class structure, which really would not
#   help since it is deep in the respective package code now. The raster
#   extent object is another possible workaround.
#
#   Arguments...
#     bbox = the bbox object
#
#   Returns...
#     TRUE if all tests were passed.
#
#Author...									Date: 1-Dec-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#
#---------------------------------------------------------------------------
#
#   these names are static and will not change...
#
    bbox.colnames = c('min','max')
    bbox.rownames = c('x','y')

#
#   some simple tests for a valid bbox object...
#
    if(!is(bbox, 'matrix'))
      stop(paste(deparse(substitute(bbox)),': bbox must be an object of class "matrix"!',sep=''))

    if(nrow(bbox)!=2 || ncol(bbox)!=2)
      stop(paste(deparse(substitute(bbox)),': bbox must have only two rows and two columns!',sep=''))

    rn = rownames(bbox)
    if(is.null(rn) || any(is.na( match(rn, bbox.rownames) )))
      stop(paste(deparse(substitute(bbox)),': bbox rownames must be "x", "y"!',sep=''))
    
    cn = colnames(bbox)
    if(is.null(cn) || any(is.na( match(cn, bbox.colnames) )))
      stop(paste(deparse(substitute(bbox)),': bbox colnames must be "min", "max"!',sep=''))

    if(any( apply(bbox,1,function(x) if(x['min'] >= x['max']) TRUE else FALSE) ))
      stop(paste(deparse(substitute(bbox)),': bbox "min" must be less than "max" for x and y!',sep=''))

    return(TRUE)
}   #bboxCheck

