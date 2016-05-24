draw.shape <-
function(shape,type='poly',col=1,...){
# shape is a shapefile from the 'shapefiles' library
  if(type=='p') points(shape$shp$shp[,2:3],...) else
  lapply(shape$shp$shp,function(x){
    if(is.null(shape$dbf$dbf$col)==F) col=shape$dbf$dbf$col[x$record]
    start <- x$parts+1
    end <- c(x$parts[-1],x$num.points)
    i <- unlist(apply(cbind(start,end),1,function(x) c(x[1]:x[2],NA)))
    if(type=='poly') polygon(x$points[i,],col=col,...) 
    else lines(x$points[i,1],x$points[i,2],type,col=col,...) 
    })
  box()
}

