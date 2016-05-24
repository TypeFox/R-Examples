whitening <- function(data, channel, method = c("PCA", "ZCA"), k = 4, r = 1, data.name) {
    if (missing(data)) 
        stop("'data' argument is not specified")
    if (!is.emg(data)) 
        stop("an object of class 'emg' is required")
    if (missing(channel)) {
        if (missing(data.name)) 
            data <- extractchannel(data) else data <- extractchannel(data, data.name = data.name)
    } else {
        if (missing(data.name)) 
            data <- extractchannel(data, channel) else data <- extractchannel(data, channel, data.name)
    }
    method <- match.arg(method)
    
    y <- data$values - mean(data$values)
    Y <- head(y, -k)
    for (i in 1:(k - 1)) Y <- cbind(Y, head(tail(y, -i), -(k - i)))
    Y <- cbind(Y, tail(y, -k))
    
    Cy <- cov(Y)
    
    res1 <- eigen(Cy)
    D <- diag(1/res1$values)
    V <- res1$vectors
    
    Y2 <- matrix(rep(tail(y, k), dim(Y)[2]), ncol = dim(Y)[2])
    
    if (method == "PCA") 
        R <- rbind(Y, Y2) %*% V %*% D^(1/2) %*% t(V) else R <- t(D^(1/2) %*% t(V) %*% t(rbind(Y, Y2)))
    
    object <- emg(R[, r], data$samplingrate, data$units, data$data.name)
    return(object)
} 
