plot.ecespa.getis <-
function(x, type="k", dimyx=NULL, xy=NULL, eps=NULL, color=NULL, 
                    contour=TRUE , points=TRUE,...){
   
   lambda <- x$ppp$n/area.owin(x$ppp$window)
   if (type=="k") zg <- x$klocalgrid
   if (type== "l") zg <- sqrt(x$klocalgrid/pi)
   if (type== "n") zg <- x$klocalgrid*lambda
   if (type== "d") zg <- sqrt(x$klocalgrid/pi)-x$R
   
   lx <- length(x$x)
   n <- x$ppp$n
   # x carries first the x coords of the points, then the xcorrds of the grid
   xgrid <- sort(unique(x$x[(n+1):lx]))
   ygrid <- sort(unique(x$y[(n+1):lx]))
   
   # original map from getis
   map0 <- im(mat=t(matrix(x$klocalgrid[(n+1):lx],x$ny,x$nx)), xcol=xgrid,yrow=ygrid)
  
  
   if(!is.null(dimyx) |!is.null(xy) |!is.null(eps) )  map0 <- as.im(interp.im, W=x$ppp$window, Z=map0, dimyx=dimyx, xy=xy,eps=eps)
  # for irregular windows, restrict estimation to the polygonal boundary 
  map0 <- map0[x$ppp$window, drop=F]  
   
   if(is.null(color)) color <- topo.colors(64)
   plot(map0, col=color, main="",...)
    
    if(contour==TRUE) contour(map0, add=TRUE)
    if(points==TRUE) points(x$ppp, pch=16, cex=0.7)
}

