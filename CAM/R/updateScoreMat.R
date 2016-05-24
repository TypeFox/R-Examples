updateScoreMat <-
function(scoreMat, X, scoreName, i, j, scoreNodes, Adj, output, numCores, maxNumParents, parsScore, intervMat, intervData)
    # new edge: from i to j
{
    p <- dim(X)[2]
    existingParOfJ <- which(Adj[,j] == 1)
    notAllowedParOfJ <- setdiff(which(scoreMat[,j] == -Inf), c(existingParOfJ,j))
    
    # if there is something left that we need to update
    if(length(existingParOfJ) + length(notAllowedParOfJ) < p-1)
    {
        # update column for j
        rowParents <- matrix(c(existingParOfJ,NA), p, length(existingParOfJ)+1, byrow = TRUE)
        rowParents[,length(existingParOfJ)+1] <- 1:p
        toUpdate <- setdiff(1:p,c(j,existingParOfJ,notAllowedParOfJ))
        if(length(existingParOfJ)< maxNumParents)
        {
            if(numCores == 1)
            {
                scoreUpdate <- mapply(computeScoreMatParallel,MoreArgs = list(rowParents = rowParents, selMat = matrix(TRUE,p,p), scoreName = scoreName, X = X, output = output, node2 = j, parsScore = parsScore, intervMat = intervMat, intervData = intervData), i = toUpdate)
            } else
            {
                scoreUpdate <- mcmapply(computeScoreMatParallel,MoreArgs = list(rowParents = rowParents, selMat = matrix(TRUE,p,p), scoreName = scoreName, X = X, output = output, node2 = j, parsScore = parsScore, intervMat = intervMat, intervData = intervData), i = toUpdate, mc.cores = numCores)
            }
        } else
        {
            scoreUpdate <- -Inf
        }
        scoreMat[toUpdate,j] <- scoreUpdate - scoreNodes[j]
    }    
    return(scoreMat)
}
