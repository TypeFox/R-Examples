selLasso <-
function(X,pars = list(),output = FALSE,k)
{
    if(output)
    {
        cat("Performing variable selection for variable", k, ": \n")
    }
    result <- list()
    p <- dim(as.matrix(X))
    if(p[2] > 1)
    {
        selVec <- rep(FALSE, p[2])
        modfitGam <- train_lasso(X[,-k],X[,k],pars)
        selVecTmp <- rep(FALSE, p[2] - 1)
        selVecTmp[which(modfitGam$model$beta != 0)] <- T
        selVec[-k] <- selVecTmp 
    } else
    {
        selVec <- list()
    }
    return(selVec)
}
