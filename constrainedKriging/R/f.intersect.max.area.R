f.intersect.max.area <- function( polygons, gpc.pixel)

### purpose: function returns a binary vector 1 = max intersection area
###          0 else 
### arguments:
###            polygons = polygon list  
###            gpc.pixel = gpc.poly object of the pixel
###           
### author: ch.hofer
###  date: 27.2.2007
{
t.intersect.area.list <- lapply( gpc.pixel,
                                function( gpc.pixel, polygons){
                                  return( area.poly(intersect(gpc.pixel, polygons)))
                                },
                                polygons
			    ) 
			    t.ia <- as.vector( t.intersect.area.list, mode = "numeric")
			    if( sum( t.ia == max(t.ia) ) == 1)
			    {
				return( as.numeric( t.ia == max(t.ia) ) )
			    }else
			    {
				return( t.ia <- rep(0, length(t.ia)))
			    }
			    
}


