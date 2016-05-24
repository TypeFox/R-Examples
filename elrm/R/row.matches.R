`row.matches` <-
function(y, X)
{
    i <- seq(nrow(X))
    j <- 0
    
    while(length(i) && (j <- j + 1) <= ncol(X))
    {
        i <- i[X[i, j] == y[j]];
    }
    
    return(i);
}

