f.point.block.cov <- function(	pixconfig,
                                locations,
                                model
			    )
### purpose: calculate the point-polygon covariance of one  polygon configuration
###          function returns n X m covariance matrix of n sample points and m polygons
###
### arguments:
###           pixconfig = list,
###                            	1. element = t.pixel.center = Koordinaten des Schwerpunktes des linken unteren Pixels
###                            	2. element = t.rowwidth = Pixelhoehe
###                            	3. element = t.colwidth = Pixelbreite
###                            	4. element = t.nrows = # Zeilen
###                            	5. element = t.nrows = # Spalten
###                            	6. element = t.pixel.poly.dim = Vector of 2 natuerliche Zahlen
###				   erste Zahl  = totale Anzahl Pixel, zweite Zahl = Anzahl Polygone
###                            	7. element = t.pixel.sum.in.poly = Vector mit m Elementen m = # Polygone
###				   Jedes Element beschreibt durch wieviele Pixel ein Polygon approximiert wird
###				   z.B c(8,3,2,10) = 1. Polygon durch 8, 2. Polygon durch 3, ...
###				8. element = sa.polygons = boolscher Vector mit m Elementen, m = # Polygone
###				    falls Element TRUE, Fl?che des Polygons kleiner als die Flaeche eines Pixels, FALSE sonst
###				9. element = polygon.centroids = m X 2 Matrix mit den Schwerpunkt-Koordinaten der m Polygone
###				10. element = posindex = m-Vector, mit dem Listen-Index der Polygon
###				11. element = t.which.poly = Liste, Laenge = totale Anzahl Pixel
###				    Beschreibt fuer jeden Pixelschwerpunkt in welchem Polygon er liegt.
###				    zb. 17 Element der Liste mit Wert 3 bedeutet, dass der 17 Pixelschwerpunkt in Polygon 3 liegt
###				    interger(0) bedeutet, dass der Pixelschwerpunkt in keinem Polygonliegt
###
###
###           locations = n x 2 matrix (coords of the sample)
###           model = list mit dem covairnace model
###
###
### output: C = n x m Kovarianzmatrix (n = Anzahl St?tzpunkte, m = Anzahl Polygone)
### author: Ch. Hofer
### date: 9-2-2010
{



## t.n.iter Anzahl Konfigurationen ?ber welche die Punkt Block Kovarianz gemittelt wird


t.n.iter <- as.list( 1:ceiling( dim(pixconfig$pixcenter)[2] / 2 ) )

t.pb.cov.matrix.list  <- lapply( t.n.iter,
    function( t.i, pixconfig, locations, model )
    {
	t.col.range <- ( 2 * ( t.i-1 ) + 1 ):( 2 * ( t.i-1 ) + 2 )
	t.pixel.center  <- matrix( pixconfig$pixcenter[ , t.col.range], ncol = 2 )

	t.n.poly  <- length( pixconfig$posindex )

	t.col.range <- ( 1 + ( t.i - 1 ) * t.n.poly):( t.i * t.n.poly )
#	pix.in.poly = n X m bin?re Zugeh?rigkeits-Matrix Element ith-Zeile, jth Spalte = 1 bedeutet: iter-Pixelschwerpunkt liegt in jtem Polygon, 0 sonst
	pix.in.poly <- as.matrix( pixconfig$pix.in.poly[ , t.col.range ], ncol = t.n.poly )

	sa.polygons <- pixconfig$sa.polygons[ t.col.range ]
	# Anzal Polygone mit einer kleineren Fl?che als die Fl?che eines Pixels
	t.poly.as.point <- sum( sa.polygons )
#
	polygon.centroids <- pixconfig$polygon.centroids

#	# St?tzpunkte-Polygon-Kovarianzmatrix mit NA gef?llt
	t.pb.cov.matrix <- matrix(nrow = nrow( locations ), ncol = t.n.poly )

#
# 	# Berechung der Punkt-Punkt Kovarianz falls die Polygon-Fl?che kleiner ist als die Pixelfl?che
#
	if( t.poly.as.point != 0 )
	{

#	# t.poly.as.point = number of polygons approximated by points
		t.dist <- as.vector(
   			f.row.dist(
                                  matrix( locations, ncol = 2),
                                  matrix( polygon.centroids[ sa.polygons == T, ], ncol = 2)
                      	)
         	)

  	t.pb.centroids.cov <- f.pp.cov( t.dist, model)

  	t.pb.cov.matrix[ , sa.polygons ] <- matrix(t.pb.centroids.cov, ncol = t.poly.as.point)

	} ## end if t.poly.as.point != 0)

#
# 	# Berechung der Punkt-Polygon Kovarianz falls die Polygon-Fl?che gr?sser ist als die Pixelfl?che
#
	if( t.poly.as.point < t.n.poly)
#	# f?r die Polygon welche durch Pixel approximiert werden
#	# t.poly.as.point = polygone die durch einen Punkt approximiert werden
#	# t.n.poly  = Polygonanzahl
	{

#	# t.poly.with.pixel Matrix # Zeilen = # Pixel des Pixelgrids, # Spalten = Polygone (bei CMCK)
#	# eine 1 bedeutet, dass des entsprechende Pixel (Zeilennummer) im entsprechenden Polygon liegt (Spaltennummer)

    	t.poly.with.pixel.matrix <- matrix(
 				pix.in.poly[,!sa.polygons],
				ncol = sum(!sa.polygons)
			    )


#	# Spalten der Matrix in Elemente einer Liste aufsplitten
	t.poly.with.pixel.list <- split(
  				t.poly.with.pixel.matrix,
				col( t.poly.with.pixel.matrix)
			    )

	t.pb.pixel.cov <- lapply(t.poly.with.pixel.list,
    			function(x, locations, t.pixel.center, pixconfig, model)
			{

			    t.locations.pix.cov <- f.point.rec.covmat(
							pxcoord  = locations[ ,1 ],
							pycoord  = locations[ ,2 ],
							crxcoord = t.pixel.center[ x == 1, 1 ],
							crycoord = t.pixel.center[ x == 1, 2 ],
							rxwidth  = pixconfig$colwidth,
   				 			rywidth  = pixconfig$rowwidth,
   			 	 			subdivisions = 1000,
   				 			param = model
						    )$result
			t.p.pixel.cov.matrix <- matrix(
				 			t.locations.pix.cov ,
    							ncol = nrow( matrix( t.pixel.center[ x == 1, ], ncol = 2 ) ),
							nrow = nrow( locations )
			    )

			    return( rowMeans( t.p.pixel.cov.matrix ) )


     			}, # end if t.poly.as.point < pixconfig$t.pixel.poly.dim[2]
			locations = locations,
       			t.pixel.center = t.pixel.center,
			pixconfig = pixconfig,
			model = model
		    ) ## end apply
	#cat("sum pb.cov:", sum( unlist( t.pb.pixel.cov ) ) , "\n")
	t.pb.cov.matrix[ , !sa.polygons ] <- unlist( t.pb.pixel.cov )

   	} ## end if t.poly.as.point < pixconfig$t.pixel.poly.dim[2]

 	return( t.pb.cov.matrix )
},
pixconfig = pixconfig,
locations = locations,
model = model
) # end lapply

t.mean.pb.cov  <- Reduce(f = "+", x = t.pb.cov.matrix.list ) / length( t.n.iter )

return( t.mean.pb.cov )

} ## end function




