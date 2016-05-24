convert.graph <-
function(dframe,mapsize,dispersal)
  {
    distxy <- pairdist(dframe[,2:3])
    diag(distxy) <- NA
    dist_m <- min(distxy, na.rm=TRUE)
    disp <- dispersal
    ma <- mean(dframe[, 4])
    SD.area <- sd(dframe[, 4])
    np <- nrow(dframe)
    distxy2 <- as.data.frame(distxy)
    names(distxy2) <- dframe[, 1]
    rownames(distxy2) <- dframe[, 1]
    dist_n <- distxy2
    A <- dframe[, 4]
    radius <- sqrt((A * 10000) / pi)
    grouping <- hclust(dist(dframe[, 2:3], method = "euclidean"), "single")
    clusters <- cutree(grouping, h=disp) 
    rg2 <- cbind(dframe, clusters)
    col1 <- rainbow(max(rg2[, "clusters"]))
    col2 <- as.data.frame(col1)
    col2[, 2] <- seq(1:max(rg2[, 6]))
    col3 <- merge_order(rg2, col2, by.x = "clusters", by.y = "V2", sort=FALSE,
                                                 keep_order=TRUE)
    col5 <- nndist (col3[, 3:4])
    species_data <- col3[,6:(ncol(col3)-1)]
    df1 <- data.frame(x = col3[, 3], y = col3[, 4], areas = col3[, 5], radius,
        cluster = col3[, 1], colour = col3[, 7], nneighbour = col5,
        ID = col3[, 2], species_data)
    name <- rep("species", ncol(df1) - 8)
    name[-1] <- paste(name[-1], 1:(length(name)-1), sep="")
    names(df1)[-1:-8] <- name
    result <- list(mapsize = mapsize, minimum.distance = dist_m, mean.area = ma, 
        SD.area = SD.area, number.patches = np,dispersal = dispersal,
        distance.to.neighbours = dist_n, nodes.characteristics = df1)
    class(result) <- "metapopulation"
    return(result)
  }
