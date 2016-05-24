SIR <-
function(S, S.hat)
    {
    
    na.fail(S)
    na.fail(S.hat)
    
    if (any(dim(S)!=dim(S.hat))) stop("'S' and 'S.hat' must have the same dimensions")
    
    if(!all(sapply(S, is.numeric))) stop("'S' must be numeric")
    if(!all(sapply(S.hat, is.numeric))) stop("'S.hat' must be numeric")
    
    S <- scale(S) 
    S.hat <- scale(S.hat)
    DM <- abs(cor(S, S.hat))
    p <- dim(DM)[1]
    
    maxcors <- numeric(p)
    
    for (i in 1:(p-1))
        {
        ind.i <- which(DM==max(DM), arr.ind=TRUE)
        maxcors[i] <- DM[ind.i[1,1], ind.i[1,2]]
        DM <- DM[-ind.i[1,1], -ind.i[1,2]]
        }
    maxcors[p] <- as.numeric(DM)
    res.in.db <- -10*log10(2-2*maxcors)
    res <- mean(res.in.db)
    return(res)
    }
