edge.graph <-
function(rl)
  {
	if (class(rl)!="landscape") 
  {
  stop(paste(rl, " should be an object of class class 'landscape'.", sep=""), call. = FALSE)
  }
    dframe <- rl$nodes.characteristics
    mapsize <- rl$mapsize
    disp <- rl$dispersal
    distxy3 <- as.data.frame(pairdist(dframe[,1], dframe[,2]))
    names(distxy3) <- dframe$ID
    rownames(distxy3) <- dframe$ID
    comb1 <- as.data.frame(cbind(as.numeric(rownames(dframe)),
    as.numeric(as.vector(dframe[,8]))))
    names(comb1)[names(comb1)=="V1"] <- "rownames"
    names(comb1)[names(comb1)=="V2"] <- "ID"
    d <- distxy3 < disp
    ind <- which(d, arr.ind = TRUE, useNames = FALSE)
    rownames(ind) <- c(1:nrow(ind))
    ind <- as.data.frame(ind)
    names(ind) <- c("row","col")
    for(i in 1:nrow(ind))
      {
        z <- ind[i,1]
        ind[i,1] <- as.numeric(rownames(d)[z])
      }
    for (j in 1:nrow(ind))
      {
        x <- ind[j,2]
        ind[j,2] <- as.numeric(colnames(d)[x])
      }
    t1 <-ind[, 1]==ind[, 2]
    t2 <-cbind(ind,t1)
    t3 <- t2[t2[, 3]==0,]
    t4 <- as.data.frame(t3[, 1:2])
    for(i in 1:nrow(t4))
      {
        t4[i, 3] <- ifelse(t4[i, 1] > t4[i, 2],t4[i, 1],t4[i, 2])
        t4[i, 4] <- ifelse(t4[i, 1] > t4[i, 2],t4[i, 2],t4[i, 1])
      }
    t5 <- t4[,3:4]
    t6 <- unique(t5)
    names(t6)[names(t6) == "V3"] <- "node A"
    names(t6)[names(t6) == "V4"] <- "node B"
    for (x in 1:nrow(t6))
      {
        a <- as.character(t6[x, 1])
        b <- as.character(t6[x, 2])
        d2 <- distxy3[a, b]
        t6[x, 3] <- d2
      }
    names(t6)[names(t6) == "V3"] <- "distance"
    ID <- 1:nrow(t6)
    t7 <- cbind(t6, ID)
    nodeA <- vector(length=nrow(t7))
    nodeB <- vector(length=nrow(t7))
    for(z in 1:nrow(t7))
      {
        tA <- as.numeric(t7[z, 1])
        tB <- as.numeric(t7[z, 2])
        lineA <- which(dframe$ID == tA)
        lineB <- which(dframe$ID == tB)
        nodeA[z] <- dframe$areas[lineA]
        nodeB[z] <- dframe$areas[lineB]
      }
    t8 <- cbind(t7, nodeA, nodeB)
    names(t8)[names(t8) == "nodeA"] <- "area ndA"
    names(t8)[names(t8) == "nodeB"] <- "area ndB"
    XA <- vector(length=nrow(t8))
    YA <- vector(length=nrow(t8))
    XB <- vector(length=nrow(t8))
    YB <- vector(length=nrow(t8))
    for(g in 1:nrow(t8))
      {
        ndA <- as.numeric(t8[g, 1])
        ndB <- as.numeric(t8[g, 2])
        lineA <- which(dframe$ID == ndA)
        lineB <- which(dframe$ID == ndB)
        XA[g] <- dframe[lineA, "x"]
        YA[g] <- dframe[lineA, "y"]
        XB[g] <- dframe[lineB, "x"]
        YB[g] <- dframe[lineB, "y"]
      }
    t9 <- cbind(t8, XA, YA, XB, YB)
    t10 <- t9[c("ID", "node A", "node B", "area ndA", "area ndB", "XA",
                "YA", "XB", "YB", "distance")]
    rownames(t10) <- 1:nrow(t10)
    return(t10)
  }
