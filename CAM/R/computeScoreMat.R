computeScoreMat <-
function(X, scoreName, numParents, output, numCores, selMat, parsScore, intervMat, intervData)
{
    
    # numParents indicates how many parents we consider. If numParents = 1 (default), then the 
    # score matrix is of dimension (p-1) x p. If numParents = 2, then the  
    # score matrix is of dimension (p-1)(p-2) x p and so on...
    #
    # scoreMat[i,j] equals the GAIN in score if we consider i being a parent of j. 
    # it should therefore be positive.
    
    p <- dim(X)[2]
    n <- dim(X)[1]
    rowParents <- t(combn(p,numParents))
    
    tt <- expand.grid(1:dim(rowParents)[1], 1:p)
    allNode2 <- tt[,2]
    allI <- tt[,1]
    if(numCores == 1)
    {
        scoreMat <- mapply(computeScoreMatParallel, MoreArgs = list(rowParents = rowParents, selMat = selMat, scoreName = scoreName, X = X, output = output, parsScore = parsScore, intervMat = intervMat, intervData = intervData), node2 = allNode2, i = allI)
    } else
    {
        scoreMat <- mcmapply(computeScoreMatParallel, MoreArgs = list(rowParents = rowParents, selMat = selMat, scoreName = scoreName, X = X, output = output, parsScore = parsScore, intervMat = intervMat, intervData = intervData), node2 = allNode2, i = allI, mc.cores = numCores)
    }
    
    scoreMat <- matrix(scoreMat,dim(rowParents)[1],p)
    # initScore[i] equals the variance of variable j. 
    initScore <- rep(NA,p)
    for(i in 1:p)
    {
        if(intervData)
        {
            X2 <- X[!intervMat[,i],]
        } else
        {
            X2 <- X
        }
        vartmp <- var(X2[,i])
        initScore[i] <- -log(vartmp)
        # scoreMat[i,j] equals the GAIN in score if we consider i being a parent of j. 
        scoreMat[,i] <- scoreMat[,i] - initScore[i]
    }
    return(list(scoreMat = scoreMat, rowParents = rowParents, scoreEmptyNodes = initScore))
}
