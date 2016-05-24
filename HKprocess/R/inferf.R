inferf <- function(data,n,add=0.001,minu = -0.999,maxu = 0.999) {
    # Set natural logarithm of eq. 10 in Tyralis and Koutsoyiannis (2014)
    logpfxfunction <- function(f1,x) {
        size = length(x)
        a1 <- ltza(f1^(0:(size-1)),x)

        -0.5*a1[4]-0.5*(size-1)*log(a1[3]*a1[1]-(a1[2])^2)+(0.5*size-1)*
        log(a1[3])
    }
    
    # Set a function to find the maximum of logpfxfunction
    optimprice <- function(x) {
        f <- function(f1) {logpfxfunction(f1,x = x)}
        
        as.vector(unlist(optimize(f = f,interval = c(-0.99999,0.99999),
        maximum = TRUE)[2]))
    }
    logM <- optimprice(data) # Maximum value of logpfxfunction
    # An Accept-Reject algorithm to simulate from phi
    accrej <- function(x,logM,add,minu = minu,maxu = maxu) {
        dist = 1/(maxu - minu)
        f <- function(f1) {logpfxfunction(f1,x = x)}
        logM1 <- logM - log(dist) + add
        u <- runif(1,min = 0,max = 1)
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