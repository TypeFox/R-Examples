f.pixelgrid <- function(
    		    polygons,
                    neighbours,
                    pixel.x.width,  ### width of the rectangle in x direction
                    pixel.y.width  ### width of the rectangle in y direction
                        )
 

## purpose: calculate the grid over the area of the target polygon and neighbours
##          in the the boundingbox
## arguments:
##
##
## author: Ch. Hofer
## date: 6. Feb 2007
{

## finden der maximalen boundingbox    
bboxInfo <- f.bbox.information( polygons = polygons, neighbours = neighbours )
t.x.range <- bboxInfo[[4]][1]    
t.y.range <- bboxInfo[[4]][2]     
    

#if(pixel.x.width != pixel.y.width){print("ATTENTION x.width and y.width are not equal")}
#if(min(c(pixel.x.width, pixel.y.width)) < 5){print("ATTENTION, mesh widht is smaller than 5 ->
#	cause a long calculation time")}

t.mesh.n.col <- ceiling( t.x.range / pixel.x.width )
t.mesh.n.row <- ceiling( t.y.range / pixel.y.width )

## t.info Erstellt das Gitter mit den Rechtecken
pixgridInfo <- precompute(nrows =  t.mesh.n.row,
                     ncols = t.mesh.n.col,
                     rowwidth =  pixel.y.width,
                     colwidth =  pixel.x.width,
                     rowsep = 0,
                     colsep = 0,
                     cat.level = 0
                     )


return( list(pixgridInfo = pixgridInfo, pixbboxInfo = bboxInfo ))

}
