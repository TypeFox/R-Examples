nLTT.plot <- function(phy,xlab="Normalized Time",ylab="Normalized Lineages", ...) {

 xy <- ltt.plot.coords(phy,backward=TRUE,tol=1e-6);
 xy[,2] <- xy[,2]/max(xy[,2]); #normalize number lineages
 
 xy[,1] <- xy[,1] + abs(min(xy[,1])); #make sure time runs from 0..T
 xy[,1] <- xy[,1] / max(xy[,1]);      #normalize time
 
 plot.default(xy, xlab = xlab, ylab = ylab, xaxs = "r", yaxs = "r", 
					type = "S", ...)
}
		
nLTT.lines <- function(phy, ...) {

  xy <- ltt.plot.coords(phy, backward=TRUE, tol=1e-6)
  xy[,2] <- xy[,2]/max(xy[,2]); #normalize number lineages
 
  xy[,1] <- xy[,1] + abs(min(xy[,1])); #make sure time runs from 0..T
  xy[,1] <- xy[,1] / max(xy[,1]);      #normalize time
  lines(xy, type = "S", ...)
}
		