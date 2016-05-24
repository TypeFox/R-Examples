f.gpc.poly <- function( polygons )

### purpose: function call returns at list where the elements are polygons of
###          the class gcp.poly
### arguments:
###            polygons, list of Polygons of Class "SpatialPolygons" "Spatial"
###
### author: ch.hofer
###  date: 15.2.2007
{


  t.gpc.polys<- lapply(polygons@polygons,
                function(x, t.cell.area){return(as(x@Polygons[[1]]@coords, "gpc.poly"))},
                )


return( t.gpc.polys )

} ## end function

