suff.stat <-
function(ModelMatrix, Table)
{
    nr <- nrow(ModelMatrix)
    S <- rep(0,nr)  
    for( h in 1:nr)
    {
       S[h] <- (Table) %*% ModelMatrix[h,]  
    }
    return(S)
}
