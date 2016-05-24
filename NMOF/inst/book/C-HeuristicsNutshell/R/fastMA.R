# fastMA.R -- version 2010-12-11
# N = length of times series
# O = order of moving average
# trials = number of runs for speed test
N <- 1000L; O <- 50L; trials <- 100L
y <- rnorm(N)

# variant 1 -- slow
MA1 <- numeric(N)
system.time({
    for (i in 1L:trials) {
        for (t in O:N) MA1[t] <- mean(y[(t-O+1L):t])
    }
})

# ... continues fastMA.R
# variant 2 -- compute mean
MA2 <- numeric(N)
system.time({
    for(i in 1L:trials) {
        for(t in O:N) MA2[t] <- sum(y[(t-O+1L):t])/O
    }
})
stopifnot(all.equal(MA1[-(1L:(O-1L))],MA2[-(1L:(O-1L))]))

# ... continues fastMA.R
# variant 3 -- updating
MA3 <- numeric(N)
system.time({
    for(i in 1L:trials) {
        MA3[O] <- sum(y[1L:O])/O
        for(t in (O+1L):N) 
            MA3[t] <- MA3[t-1L]+y[t]/O - y[t-O]/O
    }
})
stopifnot(all.equal(MA1[-(1L:(O-1L))],MA3[-(1L:(O-1L))]))

# ... continues fastMA.R
## variant 4 -- use filter
MA4 <- numeric(N)
system.time({
    for(i in 1L:trials) MA4 <- filter(y, rep(1/O,O), sides = 1L)
})
stopifnot(all.equal(MA1[-(1L:(O-1L))],MA4[-(1L:(O-1L))]))

# variant 5 -- use cumsum
MA5 <- numeric(N)
system.time({
    for(i in 1L:trials) {
        MA5 <- cumsum(y) / O
        MA5[O:N] <- MA5[O:N] - c( 0,MA5[1L:(N-O)] )
    }
})
stopifnot(all.equal(MA1[-(1L:(O-1L))],MA5[-(1L:(O-1L))]))