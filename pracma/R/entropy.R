##
##  e n t r o p y . R  (Fast) Approximate Entropy
##


approx_entropy <- function(ts, edim = 2, r = 0.2*sd(ts), elag = 1) {

    N <- length(ts)
    result <- numeric(2)

    for (j in 1:2) {
        m <- edim + j - 1
        phi <- zeros(1, N-m+1)
        dataMat <- zeros(m, N-m+1)
        for (i in 1:m)
            dataMat[i, ] <- ts[i:(N-m+i)]

        for (i in 1:(N-m+1)) {
            tempMat <- abs(dataMat - repmat(dataMat[, i, drop = FALSE], 1, N-m+1))
            boolMat <- apply(tempMat > r, 2, max)
            phi[i]  <- sum(!boolMat)/(N-m+1)
        }
        result[j] <- sum(phi)/(N-m+1)
    }

    apen <- log(result[1]/result[2])
    return(apen)
}


sample_entropy <- function(ts, edim = 2, r = 0.2*sd(ts), tau = 1) {
    stopifnot(is.numeric(ts), is.numeric(edim))

    if (tau > 1) {
        s <- seq(1, length(ts), by = tau)
        ts <- ts[s]
    }

    N <- length(ts)
    correl  <- numeric(2)
    datamat <- zeros(edim + 1, N - edim)
    for (i in 1:(edim+1))
        datamat[i, ] <- ts[i:(N - edim + i - 1)]

    for (m in edim:(edim+1)) {
        count <- zeros(1, N-edim)
        tempmat <- datamat[1:m, ]

        for (i in 1:(N-m-1)) {
            # calculate Chebyshev distance
            X <- abs(tempmat[, (i+1):(N-edim)] - 
                            repmat(tempmat[, i, drop=FALSE], 1, N-edim-i))
            dst <- apply(X, 2, max)

            # calculate Heaviside function
            d <- (dst < r)

            count[i] <- sum(d) / (N - edim)
        }

        correl[m - edim + 1] <- sum(count) / (N - edim)
    }

   return(log(correl[1]/correl[2]))
}

