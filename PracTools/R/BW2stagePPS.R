
BW2stagePPS <- function(X, pp, psuID){
    M <- length(unique(psuID))
    Ni <- table(psuID)
    cl.tots <- by(X, INDICES = psuID, FUN = sum)
    cl.vars <- by(X, INDICES = psuID, FUN = var)

    tU <- sum(cl.tots)
    S2U1 <- sum(pp * (cl.tots/pp - tU)^2)
    B2 <- S2U1 / tU^2

    ybarU <- mean(X)
    W2 <- sum(Ni^2 * cl.vars / pp) / tU^2
    S2U <- var(X)
    k <- (B2 + W2) / (S2U/ybarU^2)

    c("B2"=B2,
      "W2"=W2,
      "unit relvar" = S2U/ybarU^2,
      "B2+W2"=B2 + W2,
      "k"=k,
      "delta"=B2 / (B2 + W2))
}
