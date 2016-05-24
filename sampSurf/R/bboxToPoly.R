bboxToPoly = function(object, ...)
{
#---------------------------------------------------------------------------
#
#   This function simply takes an object's bbox and turns it into a 
#   "SpatialPolygons" object for plotting.
#
#   Arguments...
#    object = an object with a bbox slot from package sampSurf
#
#   Returns...
#    the perimeter object invisibly
#
#Author...									Date: 17-Dec-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   a check to make sure the bbox is valid first...
#
   # slotNames = names(getSlots(class(object)))
   # if(!length(grep('bbox', slotNames)))
   #   stop('"object" does not have a bbox slot!')

    if(is.matrix(object)) 
      bb = object
    else
      bb = bbox(object)
    stopifnot(bboxCheck(bb))

#
#   check for a matrix, things come out transposed if so, and someone (like me) is bound
#   to pass this a mtrix sometime...
#
 #   if(is.matrix(object)) {
 #     tmp = bb
 #     bb = t(bb)
 #     dimnames(bb) = dimnames(tmp)
 #   }

#
#   now just do the conversion...
#
    sr = rbind(bb[,'min'],
               c(bb['x','max'], bb['y','min']),
               c(bb['x','max'], bb['y','max']),
               c(bb['x','min'], bb['y','max']),
               bb[,'min'])  #closed polygon
    colnames(sr) = c('x','y')
    rownames(sr) = 1:5
    pgSR = Polygon(sr)
    perim.ID = paste('perimeter',.StemEnv$randomID(),sep=':')
    pgsSR = Polygons(list(sr=pgSR), ID=perim.ID)
    perimeter = SpatialPolygons(list(pgsSR=pgsSR)) 
  
    return(invisible(perimeter))
}   #bboxToPoly
