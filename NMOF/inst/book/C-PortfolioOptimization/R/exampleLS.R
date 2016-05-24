# exampleLS.R -- version 2010-12-21
require("NMOF")

# create random data
na <- 500L
C <- array(0.6, dim = c(na,na)); diag(C) <- 1
minVol <- 0.20; maxVol <- 0.40
Vols <- (maxVol - minVol) * runif(na) + minVol
Sigma <- outer(Vols,Vols) * C

# objective function
OF <- function(x, data) {
    xx <- as.logical(x)
    w <- x/sum(x)
    res <- crossprod(w[xx],data$Sigma[xx,xx])
    res <- tcrossprod(w[xx],res)
    res
}

# neighborhood function
neighbour <- function(xc, data) {
    xn <- xc
    p <- sample.int(data$na, data$nn, replace = FALSE)
    xn[p] <- abs(xn[p]-1L)
    # reject infeasible solution
    if( (sum(xn)>data$Ksup) || (sum(xn)<data$Kinf) ) {
        return(xc)} else return(xn)
}

# data 
data <- list(Sigma = Sigma, 
              Kinf = 30L, 
              Ksup = 60L, 
                na = na, 
                nn = 1L)

# random solution
card0 <- sample(data$Kinf:data$Ksup, 1L, replace = FALSE) 
assets <- sample.int(na, card0, replace = FALSE)
x0 <- numeric(na)
x0[assets] <- 1L

# settings
algo <- list(x0 = x0, 
      neighbour = neighbour, 
             nS = 5000L)

system.time(sol1 <- LSopt(OF, algo, data))
# result
sqrt(sol1$OFvalue)
# plot best solution over time
par(ylog = TRUE); plot(sqrt(sol1$Fmat[ ,2L]), type = "l")

# run more trials
trials <- 100L
allRes <- restartOpt(LSopt, n = trials, OF, algo = algo, data = data)
allResOF <- numeric(trials)
for (i in 1L:trials) allResOF[i] <- sqrt(allRes[[i]]$OFvalue)

# a simpler objective function
OF2 <- function(x, data) {
    xx <- as.logical(x); w <- 1/sum(x)
    res <- sum(w * w * data$Sigma[xx,xx])
    res
}