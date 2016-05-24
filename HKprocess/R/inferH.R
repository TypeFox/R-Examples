inferH <- function(data,n,add = 0.001,minu = 0.001,maxu = 0.999) {
    # Set natural logarithm of eq. 10 in Tyralis and Koutsoyiannis (2014)
    logphxfunction <- function(H,x) {
        size = length(x)
        maxlag <- size - 1
        a1 <- ltza(acfHKp(H,maxlag),x)
        
        -0.5*a1[4]-0.5*(size-1)*log(a1[3]*a1[1]-(a1[2])^2)+(0.5*size-1)*
        log(a1[3])
    }
    
    # Set a function to find the maximum of logphxfunction
    optimprice <- function(x) {
        f <- function(H) {logphxfunction(H,x = x)}

        as.vector(unlist(optimize(f = f,interval = c(0.00001,0.99999),
        maximum = TRUE)[2]))
    }
    logM <- optimprice(data) # Maximum value of logphxfunction
    # An Accept-Reject algorithm to simulate from H
    accrej <- function(x,logM,add,minu = minu,maxu = maxu) {
        dist = 1/(maxu - minu)
        f <- function(H) {logphxfunction(H,x = x)}
        logM1 <- logM - log(dist) + add
        u <- runif(1,min=0,max=1)
        logu <- log(u) + logM1
        y <- runif(1,min = minu,max = maxu)
        while (logu > f(y) - log(dist))
        {u <- runif(1,min = 0,max = 1)
         logu <- log(u) + logM1
         y <- runif(1,min = minu,max = maxu)
        }
        y
    }
    y <- c()
    for (i in 1:n) {
        y[i] <- accrej(x = data,logM,add,minu = minu,maxu = maxu)
    }
    return(y)
}