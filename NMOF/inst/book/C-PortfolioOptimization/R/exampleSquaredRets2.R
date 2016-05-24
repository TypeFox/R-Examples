# exampleSquaredReturns2.R -- version 2011-01-06
require(NMOF)

neighbourU <- function(sol, data){
    wn <- sol$w
    toSell <- wn > data$winf
    toBuy  <- wn < data$wsup
    i <- resample(which(toSell), size = 1L)
    j <- resample(which(toBuy), size = 1L)
    eps <- runif(1) * data$eps
    eps <- min(wn[i] - data$winf, data$wsup - wn[j], eps)
    wn[i] <- wn[i] - eps
    wn[j] <- wn[j] + eps
    Rw <- sol$Rw + data$R[,c(i,j)] %*% c(-eps,eps)
    return(list(w = wn, Rw = Rw))
}

OF <- function(sol, data) crossprod(sol$Rw)

# prepare data
na <- dim(fundData)[2L]; ns <- dim(fundData)[1L]
winf <- 0.0; wsup <- 0.05
data <- list(R = fundData,
            na = na, 
            ns = ns, 
           eps = 0.5/100,
          winf = winf, 
          wsup = wsup)

# random solution
w0 <- runif(data$na); w0 <- w0/sum(w0)
x0 <- list(w = w0, Rw = fundData %*% w0)
algo <- list(x0 = x0, 
      neighbour = neighbourU, 
             nS = 2000L, 
             nT = 10L,
             nD = 5000L, 
              q = 0.20)
system.time(res2 <- TAopt(OF,algo,data))
100*sqrt(crossprod(fundData %*% res2$xbest$w)/ns)