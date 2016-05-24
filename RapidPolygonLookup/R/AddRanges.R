AddRanges <- structure(function# Add xlim and ylim for each polygon
                       ### This function computes the bounding box for each polygon and adds this information
                       ### to the list. The bounding boxes can be used in various applications.
                       ### Our main motivation is for the massive PointsInPolygon search to exclude those 
                       ### polygons as candidates whose bounding box does not contain the current point.
                       (
                         poly.list ##<< polygon list with three elements: data, polys, and poly.centers
                       ){
                         poly.list$ranges <- list()
                         poly.list$ranges$x <- by(poly.list$polys[,c("X")],poly.list$polys[,c("PID")], FUN= range)
                         tmp <- names(poly.list$ranges$x)
                         poly.list$ranges$x <- do.call("rbind",poly.list$ranges$x) 
                         rownames(poly.list$ranges$x) <- tmp
                         
                         poly.list$ranges$y = by(poly.list$polys[,c("Y")],poly.list$polys[,c("PID")], FUN= range)
                         tmp <- names(poly.list$ranges$y)
                         poly.list$ranges$y <- do.call("rbind",poly.list$ranges$y) 
                         rownames(poly.list$ranges$y) <- tmp
                         return(poly.list)
                         ### Returns augmented polygon list with additional element -- "ranges" 
                       }, ex= function(){
                         
                         data(sf.polys, envir = environment())
                         sf.polys <- AddRanges(sf.polys)
                         str(sf.polys$ranges)
                         
                       }
)



