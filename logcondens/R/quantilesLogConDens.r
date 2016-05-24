quantilesLogConDens <- function(ps, res){

    if (any(ps < 0 | ps > 1)){stop("All entries of the argument ps given to quantilesLogConDens() must be in [0, 1]!\n")}

    x <- res$x
    phi <- res$phi
    Fhat <- res$Fhat
    n <- length(x)
    res <- matrix(NA, ncol = 2, nrow = length(ps))
    
    for (i in 1:length(ps)){

        p0 <- ps[i]
    
        if (p0 == 0){q <- -Inf}
        if (p0 == 1){q <- x[n]}
    
        if ((p0 > 0) && (p0 < 1)){
            n <- length(x)
            xj <- max(x[Fhat <= p0])
            j <- length(x[x <= xj])
            q <- xj + (x[j + 1] - x[j]) * qloglin((p0 - Fhat[j]) / (Fhat[j + 1] - Fhat[j]), (x[j + 1] - x[j]) * (phi[j + 1] - phi[j]))
        }
    
    res[i, ] <- c(p0, as.numeric(q))
    }
    
    colnames(res) <- c("ps", "quantile")
    return(res)
}
