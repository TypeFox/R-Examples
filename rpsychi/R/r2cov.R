r2cov <- function(sd, R){
    Dsqrt <- matrix(0, ncol=length(sd), nrow=length(sd))
    diag(Dsqrt) <- sd
    output <- Dsqrt %*% R %*% Dsqrt
    return(output)
}
