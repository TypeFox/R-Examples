f.p.in.poly <- function(polygons, pixcenter, rowwidth, colwidth)

###  purpose: function check if the point of t.grid lie in the polygons of the polygon list polygons
###          the function gives a matrix back (dim = n points x n polygons) with the values 0 = not a point
###        in the polygon and 1 = point in the polygon
###
### arguments:
###           polygons = Polygonobject Class = "SpatialPolygons" "Spatial"
###           pixcenter = n * 2 matrix with the coordinats of the points
###           rowwidth = mesh size of the pixel
###           colwidth = mesh size of the pixel
### author: Ch. Hofer
###  date: 20.2.2007
{
  

  t.n.poly <- length(polygons@polygons)
  
  #t.scatter <- (pi-3) / 10000
  
  tt <- lapply(polygons@polygons,
               function(t.polygon, pixcenter, rowwidth, colwidth){
                 
                 return(
                        point.in.polygon(
                                         pixcenter[,1], #+ rowwidth * t.scatter,  ## x-coords of the gird
                                         pixcenter[,2], #+ colwidth * t.scatter, ## y-coords of the grid
                                         t.polygon@Polygons[[1]]@coords[,1], ## x-coords of the polygon
                                         t.polygon@Polygons[[1]]@coords[,2]  ## y-coords of the polygon
                                         )
                        ) ## end return
                 
               },
               pixcenter,
               rowwidth,
               colwidth)
  
  t.point.in.polygon <- matrix(unlist(tt), ncol = t.n.poly)
  
  return(t.point.in.polygon)
}
