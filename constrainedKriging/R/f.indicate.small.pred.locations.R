f.indicate.small.pred.locations <- function(
                                            polygons,
					    pixarea
                                            ) 
### purpose: function call returns a  logical vector; TRUE = area of polygons < t.cell.area
###          and FALSE = area of polygons >= t.cell.area
### arguments:
###            polygons, list of Polygons of Class "SpatialPolygons" "Spatial"
###            t.cell.area = area of a discretisation rectangle
###
### author: ch.hofer
###  date: 15.2.2007
{
sa.polygons <- lapply( polygons@polygons,
    		function( x, pixarea )
		{#
		    return( x@area < pixarea )
		},
		pixarea
	    )
#
return( as.vector( sa.polygons,  "logical" ) )
#
} ## end function

