multscatter <- function(scatterlist, X, toshape=TRUE)
    {
    k <- length(scatterlist)
    p <- ncol(X)
    res <- array(NA,dim=c(p,p,k))
    
    if (toshape==FALSE){
          for (i in 1:k)  res[,,i] <- do.call(scatterlist[i],list(X))
    } else {
          for (i in 1:k)  {
                Mi <- do.call(scatterlist[i],list(X))
                res[,,i] <- Mi/det(Mi)^(1/p)
                }
    }
    
    res
    }
