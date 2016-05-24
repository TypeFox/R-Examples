f.bb.cov.one.configuration <- function(
    					pixconfig,
             pixcm,
             model,
              i
              )
### purpose: calculate the polygon-polygon covariance for  one given 
###          configuration pixconfig
###
### arguments:
###           pixconfig = list,
###                            	1. element = pixcenter = Koordinaten des Schwerpunktes des linken unteren Pixels
###                            	2. element = t.rowwidth = Pixelh?he
###                            	3. element = t.colwidth = Pixelbreite
###                            	4. element = t.nrows = # Zeilen
###                            	5. element = t.nrows = # Spalten
###                            	6. element = t.pixel.poly.dim = Vector of 2 nat?rliche Zahlen
###				   erste Zahl  = totale Anzahl Pixel, zweite Zahl = Anzahl Polygone
###                            	7. element = no.pix.in.poly= Vector mit m Elementen m = # Polygone
###				   Jedes Element beschreibt durch wieviele Pixel ein Polygon approximiert wird
###				   z.B c(8,3,2,10) = 1. Polygon durch 8, 2. Polygon durch 3, ...
###				8. element = sa.polygons = boolscher Vector mit m Elementen, m = # Polygone
###				    falls Element TRUE, Fl?che des Polygons kleiner als die Fl?che eines Pixels, FALSE sonst
###				9. element = t.centroids = m X 2 Matrix mit den Schwerpunkt-Koordinaten der m Polygone
###				10. element = posindex = m-Vector, mit dem Listen-Index der Polygon
###				11. element = t.which.poly = Liste, L?nge = totale Anzahl Pixel
###				    Beschreibt f?r jeden Pixelschwerpunkt in welchem Polygon er liegt.
###				    zb. 17 Element der Liste mit Wert 3 bedeutet, dass der 17 Pixelschwerpunkt in Polygon 3 liegt
###				    interger(0) bedeutet, dass der Pixelschwerpunkt in keinem Polygonliegt
###           pixcm =  p X p Kovarianzmatrix der Pixel p = Anzahl Pixel
###           t.grid.info = 
###           t.cov.spline = 
### author: Ch. Hofer
### date: 14.12.2006
{

  t.centroids <- pixconfig$polygon.centroids
  
  sa.polygons <-  pixconfig$sa.polygons
  
  ## pixcenter ist eine Matrix mit den Koordinaten der Mittelpunkte der Pixel
  ## pixcenter wird gebraucht um zu entscheiden ob ein Pixel in einem Polygon liegt
  
  t.col.range <- ( 2*i + 1 ):( 2*i + 2 )
  pixcenter  <- as.matrix( pixconfig$pixcenter[ , t.col.range], ncol = 2 )
  
  t.n.poly  <- length( pixconfig$posindex )
  
  t.col.range <- ( 1 + i * t.n.poly):( ( i+1 ) * t.n.poly )
  pix.in.poly <- as.matrix( pixconfig$pix.in.poly[ , t.col.range ], ncol = t.n.poly )
  
  ## no.pix.in.polyanzahl pixelcenter in den Polygonflaechen
  no.pix.in.poly <- apply( pix.in.poly, 2, sum)

 
  t.index <- vector()
  t.bl.bl.cov <- vector()
  
  ## Alle Kobinationen der unteren dreiecksmatrix plus der Diagonalelemente
  for( t.i in 1:t.n.poly ){ t.index <- c( t.index, t.i:t.n.poly ) }
  t.element.index <- cbind( t.index, rep( 1:t.n.poly, t.n.poly:1 ) )

  for(k in 1:dim( t.element.index )[ 1 ] )
  {
      if( no.pix.in.poly[ t.element.index[ k, 1 ] ] == 0 &&
	  no.pix.in.poly[ t.element.index[ k, 2 ] ] == 0 |
	  sa.polygons[ t.element.index[ k, 1 ] ] == T &&
	  sa.polygons[ t.element.index[ k, 2 ] ] == T |
          no.pix.in.poly[ t.element.index[ k, 1 ] ] == 0 &&
	  sa.polygons[ t.element.index[ k, 2 ] ] == T |
	  no.pix.in.poly[ t.element.index[ k, 2 ] ] == 0 &&
	  sa.polygons[ t.element.index[ k, 1 ] ] == T 
	  )
      {
	  t.dist <- f.row.dist( matrix( t.centroids[ t.element.index[ k, 1 ], ], ncol = 2 ),
	      			matrix( t.centroids[ t.element.index[ k, 2 ], ], ncol = 2 ) )
	  t.bl.bl.cov[k] <- f.pp.cov( t.dist = t.dist, model = model )
	  
      }
      if( no.pix.in.poly[ t.element.index[ k, 1 ] ] != 0 &&
	  no.pix.in.poly[ t.element.index[ k, 2 ] ] != 0 &&
	  sa.polygons[ t.element.index[ k, 1 ] ] == F &&
	  sa.polygons[ t.element.index[ k, 2 ] ] == F)
      {
	  t.tot.cov  <-  sum( pixcm[ pix.in.poly[, t.element.index[ k, 1 ] ] > 0, pix.in.poly[ , t.element.index[ k, 2 ] ] > 0])
	  ## t.tot.cov total covariance between to blocks
	  t.n.rect <- sum(pix.in.poly[, t.element.index[ k, 1 ] ] > 0 )  *  sum(pix.in.poly[, t.element.index[ k, 2 ] ] > 0)
	  ## number of rectangles which represent the two blocks
	  t.bl.bl.cov[ k ] <-    ( t.tot.cov )  / ( t.n.rect )
      }
      if( no.pix.in.poly[ t.element.index[ k,2 ] ] != 0 && 
	  sa.polygons[ t.element.index [k, 2 ] ] == F &&
	  sa.polygons[ t.element.index[ k, 1 ] ] == T 
	 |
          no.pix.in.poly[ t.element.index[ k, 1 ] ] == 0 &&
          sa.polygons[ t.element.index [k, 1 ] ] == F  &&
          sa.polygons[ t.element.index [k, 2 ] ] == F &&
	  no.pix.in.poly[ t.element.index[ k,2 ] ] != 0 )
      {
	  t.bl.bl.cov[k] <- f.point.block.cov( 	pixconfig,
	      					matrix(t.centroids[ t.element.index[ k, 1 ], ], ncol = 2 ),
						model = model
					    )[t.element.index[ k, 2 ] ]
	   print( t.bl.bl.cov[k] )
      }
      if( no.pix.in.poly[t.element.index[k,1]] != 0 &&
	  sa.polygons[t.element.index[k,1]] == F &&
	  sa.polygons[t.element.index[k,2]] == T |
      	  no.pix.in.poly[t.element.index[k,2]] == 0 &&
          sa.polygons[t.element.index[k,2]] == F  &&
          sa.polygons[ t.element.index [k, 1 ] ] == F &&
	  no.pix.in.poly[ t.element.index[ k,1 ] ] != 0)
      {
	  t.bl.bl.cov[k] <- f.point.block.cov(	pixconfig,
	      					matrix( t.centroids[ t.element.index[ k, 2 ], ], ncol = 2 ),
						model = model
					    )[ t.element.index[ k, 1 ] ]
	
      }
  } # end for (k in 1:dim ...
    
    
  t.bl.bl.cov.matrix <- matrix( 0, nrow = t.n.poly, ncol = t.n.poly )
  t.bl.bl.cov.matrix[ lower.tri( t.bl.bl.cov.matrix, diag = T ) ] <- t.bl.bl.cov
  t.bl.bl.cov.matrix <- t( t.bl.bl.cov.matrix )
  t.bl.bl.cov.matrix[ lower.tri( t.bl.bl.cov.matrix, diag = T ) ] <- t.bl.bl.cov
    

  return( list( t.bl.bl.cov.matrix ) )
              
           
} ### end of function
