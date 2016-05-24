f.pixconfig <- function(
    			polygons,
    			neighbours,
			pixgrid,
			n = 1
		    )
### purpose: this function build the pixelconfiguration for all polygon configuration in polygons
###          which is used to calculate the  block-block covariance und point-block covariance
###          
###
### arguments:
###            polygons = List of Polygons ("SpatialPolygons" "Spatial")
###            neighbours =  list of polygon neighbours
###            pxigridInfo = information about basic pixel grid  (basic = pixel grid of the largest bbox)
###            polygon.centroids = polygons centroids (matrix 2 x n)
###            gpc.polygons= list of gpc.poly objects of the polygons
###            t.def.pred.loc = optional vector with the number of the polygon fow which the pixel configuration
###                             should be calculated
###            t.seed = number to set the random seed with set.seed(t.seed)
###
### author: Ch. Hofer
### date: 14.12.2006
{
#
pixgridInfo = pixgrid$pixgridInfo
#
# Kontrolle welche Polygone/Bl?cke habe eine Fl?che die kleiner als die Fl?che eines Pixel ist
# sa.poylgons = boolean vector of small area polygons
sa.polygons <- f.indicate.small.pred.locations( polygons= polygons,
    						pixarea= pixgridInfo$rowwidth * pixgridInfo$colwidth
					    )
# Schwerpunkte der Zielpolygone
polygon.centroids <- f.centroid.polygon( polygons )
#
# aus den sp-Polygonen wird ein gpc polygon gemacht
gpc.polygons<- f.gpc.poly( polygons )
#
### the left lower corner of the pixel grid is choosed randomly therefore
### one have to set a seed to get reproducable results
# # # if( n == 1 ){
# # #  print( "Pixelgrid is fix!!!\n" )
# # # }
# # # else
# # # {
# # # cat( "The mean of ", n, "random Pixelgrids are calculated  !!!\n" )
# # # }
#
polygons.config <- lapply( as.vector( 1:length( polygons@polygons ), mode = "list" ),
			function(x, neighbours)
			{
		  	t.poly.config <- c(x, neighbours[[x]])
		    	},
			neighbours
		    ) 
### Verschiebung des Pixelzentrums entweder n = 1 fix 
###
if( n > 1){
    delta.x.pix <-  runif( 1, -0.5*pixgridInfo$colwidth, pixgridInfo$colwidth )
    delta.y.pix <-  runif( 1, -0.5*pixgridInfo$rowwidth, pixgridInfo$rowwidth )
}
else
{
    ### einstellungen f?r die ersten Simulationen !!!!!!!
    ##t.p.x.null <-  t.grid.bbox[1,1]  + 0.48*t.rowwidth  ## f?r die Simulationsexperimente
    ##t.p.y.null <- t.grid.bbox[2,1] 
    delta.x.pix <-  0.5 * pixgridInfo$colwidth## + 3.5## f?r die Simulationsexperimente
    delta.y.pix <-  0.5 * pixgridInfo$rowwidth ##  - 0.7
}
  
    pixconfig <- lapply(
		    polygons.config,
		    function( ith.poly, polygons, sa.polygons, pixgridInfo, polygon.centroids,
			     t.gpc.pred.loc,  t.delta.x.pix,t.delta.y.pix){

			     pixel.config <- f.build.pixel(
				 			posindex = ith.poly,
							polygons = polygons[ ith.poly ],
							sa.polygons = sa.polygons[ ith.poly ],
							pixgridInfo  = pixgridInfo,
							polygon.centroids = matrix( polygon.centroids[ ith.poly, ], ncol = 2 ),
							gpc.polygons = gpc.polygons[ ith.poly ],
							delta.x.pix = delta.x.pix,
							delta.y.pix = delta.y.pix
					      )
			    
			    return( pixel.config )
			},
			    polygons,
			    sa.polygons,
			    pixgridInfo,
			    polygon.centroids,
			    gpc.polygons,
			    delta.x.pix,
			    delta.y.pix
			)
			
	return( pixconfig )
} ## end function

