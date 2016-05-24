calcBias <- function(x) {
    
    totmean <- mean(x)
    devusermeans <- rowMeans(x) - totmean
    devitemmeans <- colMeans(x) - totmean
    biasmat <- matrix(0, nrow = nrow(x), ncol = ncol(x))
    biasmat <- t(apply(biasmat, 1, function(temp) temp + devitemmeans))
    biasmat <- apply(biasmat, 2, function(temp) temp + devusermeans)
    biasmat <- biasmat + totmean
    
    list(devusermeans = devusermeans, devitemmeans = devitemmeans, totmean = totmean, biasmat = biasmat)
} 
