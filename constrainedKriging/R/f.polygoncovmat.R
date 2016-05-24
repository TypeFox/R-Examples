f.polygoncovmat <- function(
                     pixconfig,
                     pixcm,
                     model,
		     n
		 )
### purpose: calculate the polygon-polygon covariance for every pixel 
###          configuration in the list pixconfig
###
###          functiion returns a list, each element is a cov matrix of
###          the appropiate polyon-polygon configuration the cov is calculated
###          in the function f.bb.cov.one.configuration.R
###
### arguments:
###           pixconfig = list,
###                            	1. element = t.pixel.center = Koordinaten des Schwerpunktes des linken unteren Pixels
###                            	2. element = t.rowwidth = Pixelh?he
###                            	3. element = t.colwidth = Pixelbreite
###                            	4. element = t.nrows = # Zeilen
###                            	5. element = t.nrows = # Spalten
###                            	6. element = t.pixel.poly.dim = Vector of 2 nat?rliche Zahlen
###				   erste Zahl  = totale Anzahl Pixel, zweite Zahl = Anzahl Polygone
###                            	7. element = t.pixel.sum.in.poly = Vector mit m Elementen m = # Polygone
###				   Jedes Element beschreibt durch wieviele Pixel ein Polygon approximiert wird
###				   z.B c(8,3,2,10) = 1. Polygon durch 8, 2. Polygon durch 3, ...
###				8. element = t.small.area.indicator = boolscher Vector mit m Elementen, m = # Polygone
###				    falls Element TRUE, Fl?che des Polygons kleiner als die Fl?che eines Pixels, FALSE sonst
###				9. element = t.centroids = m X 2 Matrix mit den Schwerpunkt-Koordinaten der m Polygone
###				10. element = t.poly.indices = m-Vector, mit dem Listen-Index der Polygon
###				11. element = t.which.poly = Liste, L?nge = totale Anzahl Pixel
###				    Beschreibt f?r jeden Pixelschwerpunkt in welchem Polygon er liegt.
###				    zb. 17 Element der Liste mit Wert 3 bedeutet, dass der 17 Pixelschwerpunkt in Polygon 3 liegt
###				    interger(0) bedeutet, dass der Pixelschwerpunkt in keinem Polygonliegt
###           pixcm =  p X p Kovarianzmatrix der Pixel p = totale Anzahl Pixel
###           t.grid.info = t.grid.info
###           t.cov.par = list with the covariance parameter
###           t.cov.spline =  spline of the point pixel cov for many distances 
### author: Ch. Hofer
### date: 28.2.2006
{
t.poly.index <- as.list( 1:length( pixconfig ) )
for( i in 0:( n - 1 ) )
{
	if(i == 0)
	{
	t.bb.cov.list <- lapply(
                   X = pixconfig,
                   FUN = f.bb.cov.one.configuration,
                   pixcm = pixcm,
                   model = model,
		   i = i
	       )
	t.poly.index <- as.list( 1:length( t.bb.cov.list) )
	}
	else
	{
	t.bb.cov.tmp <- lapply(
                   X = pixconfig,
                   FUN = f.bb.cov.one.configuration,
                   pixcm = pixcm,
                   model = model,
		   i = i
	       )

	t.bb.cov.list <- lapply( t.poly.index,
    		function( ith.poly, t.bb.cov.list, t.bb.cov.tmp, i)
    		{
	
		t.bb.cov.list[[ ith.poly ]][[ i + 1 ]] <- t.bb.cov.tmp[[ ith.poly ]][[ 1 ]] 
		return( t.bb.cov.list[[ ith.poly ]] ) 
    		},
    		t.bb.cov.list = t.bb.cov.list,
    		t.bb.cov.tmp = t.bb.cov.tmp,
		i = i
	    ) # end apply
	    rm( t.bb.cov.tmp )
	} # end else

} #end for
# Mittelwerte der Bolck-Kovarianzen    
t.bb.cov.mean.list <- lapply( t.bb.cov.list,
     function(x)
     {
 	return( Reduce(f = "+", x = x) / length(x) ) 
	

     }
 )
t.var.mean.bb.cov.list  <- lapply( t.bb.cov.list,
     function(t.cov.mat.list)
     {
 	t.m <- Reduce(f = "+", x = t.cov.mat.list) / length(t.cov.mat.list) 
	
	t.m.var <- Reduce(
	    		f = "+",
	    		x = lapply(t.cov.mat.list,
			    function( t.cov.mat, m )
			    { 
				( t.cov.mat - m )^2
			    },
			    t.m
			)
		    ) / (length(t.cov.mat.list) - 1) 
	return( sqrt( t.m.var ) )
     }
 )
# # Varianz der Bolck-Kovarianzen Mittelwerte   
# t.bb.cov.mean.var.list <- lapply( t.bb.cov.list,
#      function(x)
#      {
#  	return( Reduce(f = "var", x = x) / length(x) ) 
#      }
#  )


return( list( mean.bb.cov.mat = t.bb.cov.mean.list, var.mean.bb.cov.mat = t.var.mean.bb.cov.list)  )
} ## end function
