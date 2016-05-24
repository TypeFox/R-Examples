min_distance <-
function(rl)
  {
	if (class(rl)!="landscape") 
  {
  stop(paste(rl, " should be an object of class class 'landscape'.", sep=""), call. = FALSE)
  }
    df1 <- rl$nodes.characteristics
    disp <- rl$dispersal
    m <- df1[, c(8, 1, 2)]
    m2 <- as.matrix(m[,2:3])
    m3 <- pairdist(m2)
    m3[m3 > disp] <- 0
    m3[m3 == 0] <- NA
    m4 <- allShortestPaths(m3)
    nnodes <- ncol(m3)
    m5 <- matrix(NA, nnodes, nnodes)
    for(i in 1:nnodes)
      {
        for(j in 1:nnodes)
          {
            a <- extractPath(m4, i, j)
            m5[i, j] <- length(a) - 1
          }
      }
    m5[m5 == 1] <- 0
    m3[m3 > 0] <- 1
    m3[is.na(m3)] <- 0
    m6 <- m5 + m3
	m6[lower.tri(m6)] <- NA
    diag(m6) <- 0
    result <- as.data.frame(m6)
    names(result) <- df1$ID
    rownames(result) <- df1$ID
    result <- as.matrix(result)
    return (result)
  }
