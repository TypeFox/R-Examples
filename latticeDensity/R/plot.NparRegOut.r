plot.NparRegOut <-
function(x,...){
#
  if(class(x) != "NparRegOut"){
    stop("This function only plots objects of type NparRegOut")
    }
  EW.locs <- as.vector(x$EW.locs)
  NS.locs <- as.vector(x$NS.locs)
  NparRegMap <- as.vector(x$NparRegMap)
  boundaryPoly <- as.matrix(x$boundaryPoly)
  PointPattern <- as.matrix(x$PointPattern)
  hole.list <- x$hole.list
  big.matrix <- matrix(nrow=length(EW.locs),ncol=length(NS.locs),NparRegMap)
  contour(x=EW.locs,y=NS.locs,z=big.matrix,...)
  lines(rbind(boundaryPoly,boundaryPoly[1,]),...)
  n.hole <- length(hole.list)
  for(i in 1:n.hole){
    lines(rbind(hole.list[[i]],hole.list[[i]][1,]),...)
    }
  }
  