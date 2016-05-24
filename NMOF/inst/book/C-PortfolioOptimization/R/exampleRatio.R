# exampleRatio.R -- version 2011-01-12
require(NMOF)

# objective function
OFcmR <- function(sol,data) {
    Rw <- sol$Rw
    losses <- Rw - abs(Rw)
    gains <- Rw + abs(Rw)
    nL <- sum(losses < 0)
    nG <- sum(gains  > 0)
    vG <- sum(gains^data$eG)
    vL <- sum(abs(losses)^data$eL)
    (vL/nL) / (vG/nG)
}

neighbourUK <- function(sol, data){
    wn <- sol$w
    J <- wn > 0; K <- sum(J)
    eps <- data$eps * runif(1)
     if (K > data$Kinf && K < data$Ksup) {
        toSell <- wn > 0
        toBuy  <- wn < data$wsup
    } else {
        if (K == data$Ksup) {
            toSell <- wn > 0
            toBuy  <- J & (wn < data$wsup)
        } else {
            toSell <- wn > eps
            toBuy  <- wn < data$wsup
        }
    }
    i <- resample(which(toSell),1) 
    j <- resample(which(toBuy),1) 
    eps <- min(wn[i], data$wsup - wn[j], eps)
    wn[i] <- wn[i] - eps
    wn[j] <- wn[j] + eps
    Rw <- sol$Rw + data$R[,c(i,j)] %*% c(-eps,eps)
    list(w = wn, Rw = Rw)
}

# prepare data
na <- dim(fundData)[2L]
ns <- dim(fundData)[1L]

data <- list(R = fundData,
        na = na, ns = ns, 
        eps = 0.5/100,
        wsup = 0.1,
        eG = 2, eL = 2,
        Kinf = 10L, Ksup = 50L)

# initial solution
card0 <- sample(data$Kinf:data$Ksup, 1) 
assets <- sample.int(data$na, card0, replace = FALSE)
w0 <- numeric(data$na); w0[assets] <- 1/card0
sol0 <- list(w = w0, Rw = fundData %*% w0)

algo <- list(x0 = sol0, neighbour = neighbourUK, 
             nS = 1000L, nT = 10L,
             nD = 10000L, q = 0.9)
system.time(res <- TAopt(OFcmR,algo,data))
plot(res$Fmat[,1], type = 'l')
res$OFvalue; sum(res$xbest$w <= 1e-8); sum(res$xbest$w > 1e-8)