removepoints <-
function (rl, nr)
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
    rl_0 <- rl$nodes.characteristics
    ID2 <- rl_0$ID
    nr_select <- nrow(rl_0)-nr
    rl_1 <- rl_0[sample(1:nrow(rl_0), nr_select, replace=FALSE), ]
    rl_2 <- rl_1[sort.list(as.numeric(rownames(rl_1))), ]
    names(rl_2)[names(rl_2) == "ID2"] <- "ID"
    rl_3 <- list(mapsize=mapsize2, minimum.distance=dist_m2,
                 mean.area=mean(rl_2$areas), SD.area=sd(rl_2$areas),
                 number.patches=nrow(rl_2), dispersal=disp2,
                 nodes.characteristics=rl_2)
	class(rl_3) <- "landscape"
    rl_4 <- cluster.id(rl_3)
    rownames(rl_4$nodes.characteristics) <- 1:nrow(rl_4$nodes.characteristics)
    class(rl_4) <- "landscape"
    return(rl_4)
  }
