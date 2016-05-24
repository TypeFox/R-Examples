colorbar3d <-
  function(zlim,mycolors,ncolor){
    
    # define unique vertices
    xy <- rbind(c(12,12),c(12,10),c(10,10),c(10,12))
    
    # define all vertices
    scales <- ncolor/(zlim[2]-zlim[1])
    itidx <- NULL
    newcoord <- NULL
    for (ii in 1:ncolor) {
      y <- zlim[1] + (ii-1)/scales
      newcoord <- rbind(newcoord,cbind(xy,y),cbind(xy,y+1/scales))
      itidx <- cbind(itidx,(cbind(c(1,2,6,5),c(3,4,8,7),c(1,4,8,5),
                                  c(2,3,7,6),c(1,2,3,4),c(5,6,7,8))+(ii-1)*8))
    }
    mybar <- qmesh3d(apply(newcoord,1,asHomogeneous),itidx)
    mybar$material <- list(color=rep(mycolors,each=24))
    mybar
    
  }