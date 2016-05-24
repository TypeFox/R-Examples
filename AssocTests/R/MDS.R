## multidimensional scaling
MDS <- function(S, n, TopK, SignEigenPoint)
{   
    H <- diag(n) - matrix(data=1, nrow=n, ncol=n)/n
    B <- H %*% S %*% H

    CL <- eigen(B)
    Val <- CL$values
    Vec <- CL$vectors

    # Significante Eigenvalues
    if (is.null(TopK))
    {
#        TopK <- FindTopK (Val, n, SignEigenPoint)
        TopK <- tw(Val, n, SignEigenPoint)$SigntEigenL
	      if (TopK < 1) TopK <- 1
	      
	      TopK <- min(TopK,10)
    }
    
    res <- matrix(data=0, nrow=n, ncol=TopK)
    for (i in 1:TopK)
    {
        res[,i] <- Vec[, i]*sqrt(Val[i])
    }

    res
}
