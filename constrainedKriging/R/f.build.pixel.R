   f.build.pixel  <- function(
                      posindex,
                      polygons, #list of class"SpatialPolygons" "Spatial"
                      sa.polygons,
                      pixgridInfo,
                      polygon.centroids,
                      gpc.polygons,
                      delta.x.pix,
		      delta.y.pix)
### purpose: this function build the pixel grid for a given polygon configuration polygons
###          if a pixel is in a polygon with the sa.polygons ==  TRUE then the pixel
###          will given to the polygon with the united largest area
### arguments:
###       polygons, #list of class"SpatialPolygons" "Spatial"
###       sa.polygons, = logical vector, T means polygon area is smaller than pixel
###       pixgridInfo =
###       polygon.centroids = matrix with coords of polygon centroids
###       gpc.polygons= gpc.poly object of the polygon list polygons
###       posindex = vector with the positioni numbers of the target polygon and
###                        its neighbours
###
###        t.seed = NULLL or number to set the random seed
###          the function returns the follwing list
###          t.pixel.config = list,
###                            1. element = coords of the lower left pixel
###                            2. element = vector with first number = number
###                               of pixel second number = number of polygons
###                            3. element = logical vector, TURE means the
###                               appropriate polyon area is smaller than the
###                               pixel area
###                            4. element = vector, numbers are the pixel centers
###                               in the appropiate polygon  e.g 0, 2, 1 means 0
###                               pixels in first poly  2 pixel in second poly
###                               and 1 pixel in third
###                             5. element = matrix, coords of the centrois of the
###                                polygons
###                             6. element = list, each element is a number, the number
###                                describe the polygon in which the appropriate
###                                pixel center is, e.g 13th element hast number
###                                3 that means that the 13th pixel center is
###                                in the 3th polygon
###                             7. element = vector, numbers represent which polyon are
###                                considered e.g. c(6,3,10) means its a polygon
###                                configuration with the 6th, 3th and 10th polygon
###
###
### author: Ch. Hofer
### date: 14.12.2006
{

  t.rowwidth <- pixgridInfo$rowwidth
  t.colwidth <- pixgridInfo$colwidth
  t.mesh.n.row <- pixgridInfo$nrows
  t.mesh.n.col <- pixgridInfo$ncol
  t.n.pred.loc <- length(polygons@polygons)


  #pixcenter <- as.matrix( expand.grid( t.y.pc, t.x.pc ) )
  #pixcenter <- pixcenter[,c(2,1)]
  # bounding box of PTNC
  t.grid.bbox <- bbox(polygons)

  ## print(t.grid.bbox)
  ## da das Zentrum des Pixels bestimmt ob ein Polygon approximiert wird
  ## wird der Startpunkt des Gitters nach links unten verschoben. Andernfalls
  ## h?tte die untere linke Ecke nicht die gleiche Wahrscheinlichkeit ein Pixel
  ##punkt zu besitzten

  t.grid.bbox[1,1] <- t.grid.bbox[1,1] - 0.5 * t.colwidth
  t.grid.bbox[2,1] <- t.grid.bbox[2,1] - 0.5 * t.rowwidth

  ## Zuf?lliger Startpunkt des Gitters im linken unteren Ecke falls ein seed gesetzt wurde
  ##  t.grid.points Gitterpunkte  (Zentrum des Pixels 0.5*colwidth und 0.5*rowwidth)
  ## a ist in Zeilen (row) Richtung und b ist in Spalten (col) Richtung

  t.grid.coords <- pixgridInfo$indices.preLimits

  ## pixcenter ist eine Matrix mit den Koordinaten der Mittelpunkte der Pixel
  ## pixcenter wird gebraucht um zu entscheiden ob ein Pixel in einem Polygon liegt

  pixcenter <- cbind( t.grid.bbox[1,1] + delta.x.pix + 0.5 * t.colwidth +  t.grid.coords$bx,
		      t.grid.bbox[2,1] + delta.y.pix + 0.5 * t.rowwidth  + t.grid.coords$ax)


  # sortierung der Pixelkoordinaten, dass sie von links unten spaltenweise geordnet sind
  if( dim(pixcenter)[1] > 1)
  {
      pixcenter <- pixcenter[ order( pixcenter[,1] ), ]
  }

  if( is.null(pixcenter) ){pixcenter <- matrix(pixcenter ,ncol = 2)}
  ## testen welche Pixelzentren in welchen Polygonen liegen

  pix.in.poly <- f.p.in.poly(
		       polygons = polygons,
		       pixcenter = pixcenter,
		       rowwidth = t.rowwidth,
		       colwidth = t.colwidth
		       )

  ##no.pix.in.poly anzahl pixelcenter in den Polygonfl?chen
  ##no.pix.in.poly <- apply(pix.in.poly,2,sum)

  no.pix.in.poly   <- unlist( lapply( split( pix.in.poly, col(pix.in.poly))  , sum) )

  ## diejenigen Polygone, welche keinen Pixelcenterpunkt enthalten
  ## werden auf TRUE  gesetzt, d.h. sie werden als Punkte behandelt
  sa.polygons[ no.pix.in.poly == 0 && sa.polygons == F] <- T


for(i in 1:t.n.pred.loc)
   {
    if( sa.polygons[i] == T &&
	no.pix.in.poly[i] > 0)
    {
	if( sum(sa.polygons) > 0 )
	{
	    t.gpc.pixel <- lapply(
	          split( pixcenter, row(pixcenter) ),
		function(x, cw,rw){
		    return(
			as(
			    rbind(
				cbind(t.l.x <- x[1] - 0.5*cw, t.l.y <- x[2] - 0.5*rw), ### lower left
				cbind(t.r.x <- x[1] + 0.5*cw, t.l.y), ## lower right
				cbind(t.r.x, t.r.y <- x[2] + 0.5*rw), ### upper right
				cbind(t.l.x, t.r.y) ## upper left
			    ),"gpc.poly"
			)
		    )
		},
		cw = t.colwidth,
		rw = t.rowwidth
	    )
	} # end if( sum(sa.polygons....
	#Liste mit den gemeinsamen Fl?chen
	t.inter.area.list <- lapply(
	    gpc.polygons,
	    FUN = f.intersect.area,
	    t.gpc.pixel[ pix.in.poly[,i] > 0 ]
	)

	t.inter.area <- matrix(
	    unlist(t.inter.area.list),
	    ncol = length(t.inter.area.list)
	)

	t.cut.poly <- t( apply(t.inter.area, 1,
		function(x)
		{
		    if( ( t.max <- which.max(x) ) != which.min(x) )
		    {
			x[1:length(x)] <- 0
			x[t.max] <- 1
		    }
		    else
		    {
			x[1:length(x)] <- 0
		    }
		    return(x)
		}
	    )#end apply
	)

	pix.in.poly[ pix.in.poly[,i] > 0, ] <- t.cut.poly

    } ## end if
  }# end for


  return( list(
	  pixcenter = pixcenter,
	  rowwidth = t.rowwidth,
	  colwidth = t.colwidth,
	  nrows = t.mesh.n.row,
	  ncols = t.mesh.n.col,
	  #pixel.poly.dim = dim(pix.in.poly),
	  no.pix.in.poly = no.pix.in.poly,
	  sa.polygons = sa.polygons,
	  polygon.centroids = matrix(polygon.centroids, ncol = 2),
	  posindex = posindex,
	  #which.poly = which.poly,
	  pix.in.poly = pix.in.poly
	  )
     )

  } ## end function

