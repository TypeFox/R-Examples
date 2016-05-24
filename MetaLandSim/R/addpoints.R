addpoints <-
function (rl, nr)
  {
if (class(rl)!="landscape") 
  {
   stop(paste(rl, " should be an object of class class 'landscape'.", sep=""), call. = FALSE)
  }	
    if(nr != 0)
      {
        mapsize2 <- rl$mapsize
        dist_m2 <- rl$minimum.distance
        areaM2 <- rl$mean.area
        areaSD2 <- rl$SD.area
        Npatch2 <- rl$number.patches
        disp2 <- rl$dispersal
        rl <- rl$nodes.characteristics
        ID2 <- rl$ID
        rl2 <- rl[,1:2]
        wind <- owin(xrange=c(0, mapsize2), yrange=c(0, mapsize2))
        suppressWarnings(pts_0 <- as.ppp(rl2, W = wind, fatal=TRUE))
        pts_1 <- rSSI(r=dist_m2, n = npoints(pts_0)+nr, win = wind, 
                      giveup = 1000, x.init=pts_0)
        df_pts0 <- as.data.frame(coords(pts_1))
        df_pts1 <- as.data.frame(coords(pts_0))
        df_pts2 <- df_pts0[!duplicated(rbind(df_pts1, df_pts0))[nrow(df_pts1) + 
                           1:nrow(df_pts0)],]
        nrow_0 <- nrow(rl)
        na_lines <- as.data.frame(matrix(NA, nrow=nr, ncol=ncol(rl)))
        colnames(na_lines) <- colnames(rl)
        rownames(na_lines) <- max(as.numeric(rownames(rl))) + 1:nrow(na_lines)
        rl <- rbind(rl, na_lines)
        rl[(nrow_0+1):nrow(rl), 1:2] <- df_pts2
        areas0 <- abs(rnorm(nr, mean = areaM2, sd = areaSD2))
        radius0 <- sqrt((areas0 * 10000) / pi)
        rl[(nrow_0 + 1):nrow(rl), "areas"] <- areas0
        rl[(nrow_0 + 1):nrow(rl), "radius"] <- radius0
        new_ID <- (max(ID2) + 1) : (max(ID2) + nr)
        rl[, 8] <- as.character(c(ID2, new_ID))
        grouping <- hclust(dist(rl[, 1:2], method = "euclidean"), "single")
        clusters <- as.data.frame(cutree(grouping, h=disp2))[, 1]
        rl[, "cluster"] <- clusters
        col1 <- rainbow(max(rl[, 5]))
        col2 <- as.data.frame(col1)
        col2[, 2] <- seq(1:max(rl[, 5]))
        col3 <- merge_order(rl, col2, by.x = "cluster",
                                 by.y = "V2", sort=FALSE,
                                 keep_order=TRUE)[, 9]
        rl[, "colour"] <- col3
        col4 <- nndist (rl[, 1:2])
        rl[, "nneighbour"] <- col4
        rland.out <- list(mapsize=mapsize2, minimum.distance=dist_m2,
                          mean.area=mean(rl$areas), SD.area=sd(rl$areas),
                          number.patches=nrow(rl), dispersal=disp2,
                          nodes.characteristics=rl)
        class(rland.out) <- "landscape"
        rland.out <- cluster.id(rland.out)
      }
    if(nr == 0)
      {
        rland.out <- rl
      }
    rownames(rland.out$nodes.characteristics) <- 1:nrow(rland.out$nodes.characteristics)
    class(rland.out) <- "landscape"
	return(rland.out)
  }
