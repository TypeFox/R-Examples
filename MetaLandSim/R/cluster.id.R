cluster.id <-
function (rl)
  {
	if (class(rl)!="landscape") 
  {
  stop(paste(rl, " should be an object of class class 'landscape'.", sep=""), call. = FALSE)
  }
    mapsize2 <- rl$mapsize
    dist_m2 <- rl$minimum.distance
    areaM2 <- rl$mean.area
    areaSD2 <- rl$SD.area
    Npatch2 <- rl$number.patches
    disp2 <- rl$dispersal
    rl2 <- rl$nodes.characteristics
    ID2 <- rl2$ID
    rl3 <- rl2[,c(1,2,3,4,8)]
    if(nrow(rl3) >= 2)
      {
        grouping <- hclust(dist(rl3[,1:2], method="euclidean"), "single")
        clusters <- cutree(grouping, h=disp2)
      } else clusters <- 1
    new_2 <- cbind(rl3, clusters)
    col1 <- rainbow(max(new_2[,6]))
    col2 <- as.data.frame(col1)      
    col2[,2] <- seq(1:nrow(col2))   
    col3 <- merge_order(new_2, col2, by.x="clusters",
                             by.y="V2", sort=FALSE, keep_order=TRUE) 
    col5 <- nndist (rl3[,1:2])      
    col4 <- data.frame(col3$x, col3$y, col3$areas, col3$radius, 
                       col3$clusters, col3$col1, col5,
                       as.numeric(col3$ID))
    names(col4)[names(col4)=="col3.x"] <- "x"
    names(col4)[names(col4)=="col3.y"] <- "y"
    names(col4)[names(col4)=="col3.areas"] <- "areas"
    names(col4)[names(col4)=="col3.radius"] <- "radius"
    names(col4)[names(col4)=="col3.clusters"] <- "cluster"
    names(col4)[names(col4)=="col3.col1"] <- "colour"
    names(col4)[names(col4)=="col5"] <- "nneighbour"
    names(col4)[names(col4)=="as.numeric.col3.ID."] <- "ID"
    rownames(col4) <- 1:nrow(col4)
    rland.out <- list(mapsize=mapsize2, minimum.distance=dist_m2, 
                      mean.area=mean(col4$areas), SD.area=sd(col4$areas), 
                      number.patches=nrow(col4), dispersal=disp2,
                      nodes.characteristics=col4)
    class(rland.out) <- "landscape" 
    return(rland.out)
  }
