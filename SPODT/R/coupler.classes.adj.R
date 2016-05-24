coupler.classes.adj <-
function(mat)
{
    n <- nrow(mat)
    v <- 1:n
    mat <- mat - diag(1, nrow=n, ncol=n)
    
    l <- lapply(v, coupler.classes.adj.i, mat=mat)
    lg <- lapply(l, length)
    nm <- rownames(mat)

    nm <- rep(nm, lg)  
    l <-  unlist(l)
        
    couples <- matrix(c(nm, l), nrow=length(l), ncol=2)
    
    nm <- as.integer(nm)   
    l <-  as.integer(l)
    
    ind <- which(nm < l)
    
    couples <- couples[ind,]

    return(couples)
}
