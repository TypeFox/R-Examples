GelmanRubin <- function(theta){
    ## theta is a matrix of outputs from various chains

    if(!is.matrix(theta)){
        stop("theta must be a matrix")
    }

    nObs <- nrow(theta)
    nCols <- ncol(theta)
    n1 <- floor(nObs*0.5)
    n2 <- nObs - n1

    if(nObs<100)
        stop("There must be at least 100 observations from each chain")

    if(nCols<2)
        stop("There must be at least two chains")

    theta <- theta[-(1:n1),] # take only the second half of the data

    vars <- apply(theta, 2, var)
    means  <- apply(theta, 2, mean)
    mBar <- mean(means)

    B <- n2*sum((means-mBar)^2)/(nCols-1)
    W <- sum(vars)/nCols
    sigmaSq <- ((n2-1)*W+B)/(n2)
    vHat <- sigmaSq+B/(n2*nCols)
    df <- n2
    R <- sqrt(vHat/W*(df/(df-2)))

    results.df <- data.frame(n = n2,B,W,vHat,R)
    cat(paste(R,"\n"))
    invisible(results.df)
}

GR<-function(theta){
    return(GelmanRubin(theta))
}
