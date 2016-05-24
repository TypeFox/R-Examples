extract.graph <-
function (rl,rlist,nr)
  {
	if (class(rl)!="landscape") 
  {
  stop(paste(rl, " should be an object of class class 'landscape'.", sep=""), call. = FALSE)
  }
    mapsize <- rl$mapsize
    dist_m <- rl$minimum.distance
    disp <- rl$dispersal
    dist_m <-rl$minimum.distance
    rl1 <- rlist[[nr]]
    rl2 <- list(mapsize=mapsize, minimum.distance=dist_m, mean.area=mean(rl1$areas),
                     SD.area=sd(rl1$areas), number.patches=nrow(rl1),
                     dispersal=disp, nodes.characteristics=rl1)
	class(rl2) <- "landscape"
    if(nrow(rl1)>1)rl3 <- cluster.id(rl2)
    if(nrow(rl1)==1)rl3 <- rl2
    return(rl3)
  }
