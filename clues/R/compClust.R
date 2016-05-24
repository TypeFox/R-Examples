# v0.2.6 created on Feb. 10, 2009 by Weiliang Qiu
#
# compare different partitions for the same data set
# Purposes:
#  (1) get some senses if the cluster structure is overlapped, touched, separated, or
#      well-separated. If different partitions are quite similar based on
#      agreement indices, then the cluster structure is at least separated
#      and the obtained partitions are quite reliable. If different partitions
#      are quite different, then the cluster structure is probably overlapped or
#      not regular shape. In this, we can check the average silhouette index and/or
#      CH index to decide which partition get a relatively more separated cluster 
#      structures.  Of course cluster analysis is just an exploratory analysis.
#      Hence the function 'compClust' only provides some useful information to the
#      user/decision maker for their reference. 
#  (2) get some senses if some clustering methods performs similarly
#  (3) if true cluster membership is known (e.g. in simulation study), 
#      the function 'compClust' could check which cluster method performs the best
#      for a data set. 

# y -- n x p data matrix; rows are data points; columns are variables;
#        n -- number of data points; p -- number of variables
# memMat -- n x k membership matrix; k is the number of different partitions 
#           for the data matrix 'y'
# disMethod -- measure of the dissimilarity between data points
#
compClust <- function(y, memMat, disMethod = "Euclidean")
{
    disMethod <- match.arg(disMethod, c("Euclidean", "1-corr"))
 
    nc <- ncol(memMat)
 
    str <- colnames(memMat)
 
    avg.sVec <- rep(0, nc)
    names(avg.sVec) <- str
 
    CHVec <- rep(0, nc)
    names(CHVec) <- str
 
    RandMat <- matrix(1, nrow = nc, ncol = nc)
    rownames(RandMat) <- str
    colnames(RandMat) <- str
 
    HAMat <- matrix(1, nrow = nc, ncol = nc)
    rownames(HAMat) <- str
    colnames(HAMat) <- str
 
    MAMat <- matrix(1, nrow = nc, ncol = nc)
    rownames(MAMat) <- str
    colnames(MAMat) <- str
 
    FMMat <- matrix(1, nrow = nc, ncol = nc)
    rownames(FMMat) <- str
    colnames(FMMat) <- str
 
    JaccardMat <- matrix(1, nrow = nc, ncol = nc)
    rownames(JaccardMat) <- str
    colnames(JaccardMat) <- str
 
    # obtain the average silhouette index and CH index for each partition
    for(i in 1:nc)
    {
        CHVec[i] <- get_CH(y, as.numeric(memMat[, i]), disMethod)
        avg.sVec[i] <- get_Silhouette(y, as.numeric(memMat[, i]), disMethod)$avg.s
        # get agreement index matrix for each agreement index
        i1 <- i + 1
        if(i1 <= nc)
        { for(j in i1:nc)
            {
                RandMat[i, j] <- adjustedRand(as.numeric(memMat[, i]), 
                  as.numeric(memMat[, j]), randMethod = "Rand")
                RandMat[j, i] <- RandMat[i, j]
               
                HAMat[i,j] <- adjustedRand(as.numeric(memMat[, i]), 
                  as.numeric(memMat[, j]), randMethod = "HA")
                HAMat[j, i] <- HAMat[i, j]
               
                MAMat[i, j] <- adjustedRand(as.numeric(memMat[, i]), 
                  as.numeric(memMat[, j]), randMethod = "MA")
                MAMat[j, i] <- MAMat[i, j]
               
                FMMat[i, j] <- adjustedRand(as.numeric(memMat[, i]), 
                  as.numeric(memMat[, j]), randMethod = "FM")
                FMMat[j, i] <- FMMat[i, j]
               
                JaccardMat[i, j] <- adjustedRand(as.numeric(memMat[, i]), 
                  as.numeric(memMat[, j]), randMethod = "Jaccard")
                JaccardMat[j, i] <- JaccardMat[i, j]
            }
        }
    }
 
    res <- list(avg.s = avg.sVec, CH = CHVec, Rand = RandMat, 
      HA = HAMat, MA = MAMat, FM = FMMat, Jaccard = JaccardMat)
 
    return(res)
}

