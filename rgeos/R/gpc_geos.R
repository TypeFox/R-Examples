rgeos_SpatialPolygons2gpcpoly <- function(from) {
    
    if (!inherits(from,"SpatialPolygons"))
        stop("sp does not inherit from SpatialPolygons")
    
    gpcs = list()
    for(i in 1:length(from@polygons)) {
        polys = from@polygons[[i]]
        
        pts = list()
        for(j in 1:length(polys@Polygons)) {
            coords = polys@Polygons[[j]]@coords
            hole = polys@Polygons[[j]]@hole
            l=nrow(coords)
            
            pts[[j]] = list(x=coords[-l,1],y=coords[-l,2],hole=hole)
        }
        
        gpc = new("gpc.poly", pts = pts)
        gpcs[i] = gpc
    }
    
    if (length(gpcs) == 0)
        gpcs = NULL
    if (length(gpcs) == 1)
        gpcs = gpcs[[1]]

    return(gpcs)
}

rgeos_gpcpoly2SpatialPolygons <- function(from) {
    
    if (!is.list(from))
        from = list(from)
        
    res=lapply(from, function(gpc) {

        if (!inherits(gpc,"gpc.poly"))
            stop("object does not inherit from gpc.poly")
    
        if (length(gpc@pts) < 1)
            stop("must be at least one polygon")
    
        polylist = lapply( gpc@pts, function(pt) {
            x=pt$x
            y=pt$y
            l=length(x)
            
            if (x[1]!=x[l] | y[1]!=y[l]) {
                x = c(x,x[1])
                y = c(y,y[1])
            }
            
			hole = 
			
            return( Polygon(cbind(x,y),pt$hole) )
        })

        
        p = Polygons(polylist,"0")
        attr(p,"comment") = createPolygonsComment(p)
        
        return(p)
    })
    
    for (m in 1:length(from)) {
        res[[m]]@ID = as.character(m)
    }
    
    return( SpatialPolygons(res) )
}

setAs("SpatialPolygons", "gpc.poly", rgeos_SpatialPolygons2gpcpoly)
setAs("gpc.poly", "SpatialPolygons", rgeos_gpcpoly2SpatialPolygons)

setAs("SpatialPolygons", "gpc.poly.nohole", rgeos_SpatialPolygons2gpcpoly)
setAs("gpc.poly.nohole", "SpatialPolygons", rgeos_gpcpoly2SpatialPolygons)


areaGPC <- function(x.mat) {
    if(nrow(x.mat) < 3) 
        return(0);   
    x.segmat <- cbind(x.mat, rbind(x.mat[2:nrow(x.mat), ],
         x.mat[1, ]));
    abs(sum(x.segmat[,1] * x.segmat[,4] - x.segmat[,3]
        * x.segmat[,2])) / 2
}
