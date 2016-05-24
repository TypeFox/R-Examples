f.bbox.information <- function( polygons, neighbours )
## purpose: find the largest bounding box of a polygon an its neighbours and find
##          the polygon with the smalest area
##
## arguments:
##           polygons = list of polyogons ("SpatialPolygons" "Spatial"    )
##           neighbours = list each element is a vector with the polygons indices of the neighbours 
##
##
##           wholes in polygons are not considered !!!!!!!!
##
##
## author: Ch. Hofer
## date: 29.1.2007
{
t.n.poly <- as.list( 1:length(polygons@polygons) )
#
t.bbox.list <- lapply(t.n.poly, 
    function(i, polygons, neighbours)
    {
	t.bbox.surrounding <- polygons[c(i, neighbours[[i]])]@bbox
	
	return( list(
		t.x.range.tmp <- diff(t.bbox.surrounding[1,]),
		t.y.range.tmp <- diff(t.bbox.surrounding[2,]),
		t.poly.area <- polygons@polygons[[ i ]]@Polygons[[1]]@area,
		t.bbox.area <- diff(t.bbox.surrounding[1,]) * diff(t.bbox.surrounding[2,])
	    )
	)
    },
    polygons,
    neighbours
)
#
t.bbox.matrix <- matrix( unlist( t.bbox.list), ncol = 4, byrow =T )
t.which.max.bbox <- apply(t.bbox.matrix, 2, which.max)
t.which.min.bbox <- apply(t.bbox.matrix, 2, which.min)
#
return(	list(
	largest.box.index = t.which.max.bbox[4], 
	min.area.index = t.which.min.bbox[3], 
	min.area = t.bbox.matrix[t.which.min.bbox[3],3], 
	t.xy.range = cbind(
	    t.bbox.matrix[t.which.max.bbox[1] ,1], 
	    t.bbox.matrix[t.which.max.bbox[2] ,2])
    )
)
} #end function
