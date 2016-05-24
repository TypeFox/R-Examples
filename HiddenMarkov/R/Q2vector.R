Q2vector <- function(Q){
    #   map mmpp Q rate matrix to vector
    #   use log like constraints
    m <- nrow(Q)
    p <- log(as.vector(Q[as.logical(1-diag(m))]))
    return(p)
}
