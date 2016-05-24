fisherinfo <- function(beta, X, risksetlist, event){
    n <- length(event)
    P <- length(beta)
    f <- as.numeric(X %*% beta)
    ef <- exp(f)
    info <- matrix(nrow = P, ncol = P, 0)
    index <- which(event == 1)
    for(p in 1:P){
        for(q in 1:P){
            part1 <- part2 <- rep(0, n)
            for(i in index) {
                j <- risksetlist[[i]]
                ef.j <- ef[j]
                risk <- sum(ef.j)
                X.j.p <- X[j, p]
                X.j.q <- X[j, q]
                part1[i] <- sum(ef.j * X.j.p * X.j.q)/risk
                part2[i] <- sum(ef.j * X.j.p) * sum(ef.j * X.j.q)/(risk * risk)
            }
            info[p, q] <- sum(event * part1) - sum(event * part2)
        }
    }
    return(info)}