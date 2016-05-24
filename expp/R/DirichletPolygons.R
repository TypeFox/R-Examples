setGeneric("DirichletPolygons", function(x, boundary, ...)   standardGeneric("DirichletPolygons") )


.DirichletPolygons <- function(x, boundary) {
		# x is a SpatialPointsBreeding, boundary is a SpatialPolygons*
            coords = coordinates(x)
            
			p  =  tile.list(deldir(coords[,1], coords[,2], suppressMsge = TRUE))
			p  = lapply(p, function(P) data.frame(x = P$x, y = P$y) )
			
			d  =  do.call(rbind, p)    
			d$id = rep(x@id, sapply(p, nrow) )
			
			P = split(d, d$id)
			P = lapply(P, function(x) { x = 
							x = rbind(x, x[1, ])  
							row.names(x) = paste(x$id, 1:nrow(x), sep = "_")
							x
							} )
			
			P = lapply(P,  function(x) SpatialPolygons(list( Polygons(list(Polygon(x[, c('x', 'y')])), ID = x$id[1] ) )) )
			P = lapply(P, function(pp) { proj4string(pp) = CRS(proj4string(x)); pp } )
      
			P = lapply(P,function(pj) gIntersection(boundary, pj, id = slot(slot(pj, 'polygons')[[1]], 'ID')))
			P = do.call(rbind, P) 
			P = SpatialPolygonsDataFrame(P, data = data.frame(ID = x@id, row.names = x@id))
			
		P
}

setMethod("DirichletPolygons",  
          signature  = c(x = "SpatialPointsBreeding", boundary = "missing"), 
          definition = function(x, ...) {
			if( !require(spatstat, quietly = TRUE) ) stop('spatstat package is not available.')
            coords = coordinates(x)
            ids = x@id
            rr = ripras(coords, shape = "convex", ...)
            rr = cbind(x = rr$bdry[[1]]$x, y = rr$bdry[[1]]$y)
            boundary =  SpatialPolygons(list( Polygons(list( Polygon(rbind(rr, rr[1, ] )) ) , 1) ) )
            proj4string(boundary) = proj4string(x)
            
            .DirichletPolygons(x, boundary)
            
          }
)

setMethod("DirichletPolygons",  
          signature  = c(x = "SpatialPointsBreeding", boundary = "integer"), 
          definition = function(x, boundary, width) {
            
                z = data.frame(coordinates(x), id = x@id )

                bb = data.frame(id = boundary, o = 1:length(boundary) )
                bb = merge(bb, z, by = 'id')
                bb = bb[order(bb$o), ]
                bb = rbind(bb, bb[1, ])

                P = readWKT( paste( "POLYGON((", paste(paste(bb$x, bb$y), collapse = ','), "))" ) )

                if( missing(width) ) {
                        # median distance between points
                        z12 =  cbind(bb[-nrow(bb), c('x', 'y')], bb[-1, c('x', 'y')])
                        width = mean(apply(z12, 1, function(x) spDists( as.matrix(t(x[1:2])), as.matrix(t(x[3:4])) ) ) )/2
                        }
                
                P = gBuffer(P, width = width)         
                
				.DirichletPolygons(x, P)
            
          }
)


setMethod("DirichletPolygons",  
          signature  = c(x = "SpatialPointsBreeding", boundary = "SpatialPolygons"), 
          definition = function(x, boundary) {
            
			.DirichletPolygons(x, boundary)
        
            
          }
)


























