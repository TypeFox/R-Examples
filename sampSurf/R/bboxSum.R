bboxSum = function(arr.bbox, ...)
{
#---------------------------------------------------------------------------
#
#   This little routine will calculate an overall bounding box matrix from
#   a 3D array of bbox matrices. Note that the array should be set up so
#   rows and columns are individual bbox matrices, and the third dimension
#   is the index to the page of individual matrices. See downLogIZs for an
#   example of its use.
#
#   Arguments...
#     arr.bbox = the array of bbox matrices
#
#   Returns...
#     a bbox matrix with overall min/max values
#
#Author...									Date: 31-Mar-2011
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   catch some errors...
#
    nn = dim(arr.bbox)
    if(length(nn)!=3)
      stop('Invalid array of bbox objects passed!')

    ce = tryCatch(apply(arr.bbox,3,bboxCheck), error=function(e)e)
    if( is(ce, 'simpleError') ) 
       stop('One or more bbox objects poorly formed... ', ce$message)

#
#   new overall bbox is simple...
#
    amin = apply(arr.bbox, 1, min)
    amax = apply(arr.bbox, 1, max)
    bbox = cbind(amin, amax)
    dimnames(bbox) = dimnames(arr.bbox)[-3]

    return(bbox)
}   #bboxSum

    
    
