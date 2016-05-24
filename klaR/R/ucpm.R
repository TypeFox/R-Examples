ucpm <- function(m, tc, ec = NULL)
{
# membership values
# true classes
# estimated classes
    klassen <- function(y){
        n <- length(y)
        m <- length(table(y))
        leer <- numeric(n*m)
        Y <- matrix(leer, nrow=n, ncol=m)
        for(i in 1:n) Y[i,y[i]] <- 1
        return(Y)
    }
    if(is.null(colnames(m))) colnames(m) <- levels(tc)
    membercheck(m)
    G <- ncol(m)
    N <- nrow(m)
    dummy <- matrix(0, nrow=N, ncol=G)
    if (is.null(ec)) 
        ec <- factor(max.col(m), levels = seq(along = colnames(m)), 
            labels = colnames(m))
    CR <- mean(ec == tc)
    e.ec <- klassen(ec)
    e.c <- klassen(tc)
    dummy <- apply((e.ec - m)^2, 1, sum)
    AS <- 1 - sqrt(G) / (N * sqrt(G-1)) * sum(sqrt(dummy))
    dummy <- apply((e.c - m)^2, 1, sum)
    AC <- 1 - sqrt(G) / (N * sqrt(G-1)) * sum(sqrt(dummy))
    CF <- mean(apply(m,1,max))
    CFvec <- apply(m, 1, max)
    CFvec <- unlist(lapply(split(CFvec,tc), mean))
    return(list(CR=CR, AC=AC, AS=AS, CF=CF, CFvec=CFvec))
}
