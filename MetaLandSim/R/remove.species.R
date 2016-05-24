remove.species <-
function(sp)
  {
    if (class(sp)!="metapopulation") 
      {
        stop(paste(sp, " should be an object of class class 'metapopulation'.", sep=""), call. = FALSE)
      }
    l1 <- sp$mapsize
    l2 <- sp$minimum.distance 
    l3 <- sp$mean.area 
    l4 <- sp$SD.area
    l5 <- sp$number.patches
    l6 <- sp$dispersal 
    l7 <- sp$nodes.characteristics[,1:8]
    rland.out <- list(mapsize=l1, minimum.distance=l2, mean.area=l3,
                               SD.area=l4, number.patches=l5,dispersal=l6,nodes.characteristics=l7)
    class(rland.out) <- "landscape"
    return(rland.out)
  }
