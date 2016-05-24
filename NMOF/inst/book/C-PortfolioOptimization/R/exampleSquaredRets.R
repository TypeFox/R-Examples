# exampleSquaredReturns.R -- version 2011-04-01
require(NMOF)
resample <- function(x, ...) x[sample.int(length(x), ...)]

na <- dim(fundData)[2L]
ns <- dim(fundData)[1L]
winf <- 0.0
wsup <- 0.05
data <- list(R = t(fundData),
            RR = crossprod(fundData),
            na = na, 
            ns = ns, 
           eps = 0.5/100,
          winf = winf, 
          wsup = wsup)

neighbour <- function(w, data){
    eps <- runif(1) * data$eps
    toSell <- w > data$winf
    toBuy  <- w < data$wsup
    i <- resample(which(toSell), size = 1L)
    j <- resample(which(toBuy), size = 1L)
    eps <- min(w[i] - data$winf, data$wsup - w[j], eps)
    w[i] <- w[i] - eps
    w[j] <- w[j] + eps
    w
}

OF1 <- function(w, data) {
    Rw <- crossprod(data$R,w)  
    crossprod(Rw)
}
OF2 <- function(w, data) {
    aux <- crossprod(data$RR,w) 
    crossprod(w,aux)
}

# random solution
w0 <- runif(na); w0 <- w0/sum(w0)

algo <- list(x0 = w0, 
      neighbour = neighbour, 
             nS = 2000L, 
             nT = 10L,
             nD = 5000L, 
              q = 0.20)
system.time(res <- TAopt(OF1,algo,data))
100*sqrt(crossprod(fundData %*% res$xbest)/ns)
system.time(res <- TAopt(OF2,algo,data))
100*sqrt(crossprod(fundData %*% res$xbest)/ns)

# benchmark
require(quadprog)
covMatrix <- crossprod(fundData)
A <- rep(1, na); a <- 1
B <- rbind(-diag(na),
            diag(na))
b <- rbind(array(-data$wsup, dim = c(na,1)),
           array( data$winf, dim = c(na,1)))
system.time({
    result <- solve.QP(Dmat = covMatrix, 
                       dvec = rep(0,na),
                       Amat = t(rbind(A,B)),
                       bvec = rbind(a,b),
                        meq = 1)
})
wqp <- result$solution 
# compare results
100 * sqrt( crossprod(fundData %*% wqp)/ns )
100 * sqrt( crossprod(fundData %*% res$xbest)/ns )
# check constraints
min(res$xbest); max(res$xbest); sum(res$xbest)  # TA
min(wqp); max(wqp); sum(wqp)  # QP 