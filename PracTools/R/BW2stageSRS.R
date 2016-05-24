
BW2stageSRS <- function(X, psuID){
    M <- length(unique(psuID))
    Ni <- table(psuID)
    Nbar <- length(X)/M

    ti <- by(X, INDICES = psuID, FUN = sum)
    S2Ui <- by(X, INDICES = psuID, FUN = var)

    tbarU <- mean(ti)
    tU <- M*tbarU
    S2U1 <- var(ti)
    B2 <- S2U1 / tbarU^2

    ybarU <- mean(X)
    S2U <- var(X)

    W2 <- M * sum(Ni^2 * S2Ui) / tU^2

    c("B2"=B2,
      "W2"=W2,
      "unit relvar"=S2U/ybarU^2,
      "B2+W2"=B2 + W2,
      "k"=(B2 + W2)/(S2U/ybarU^2),
      "delta" = B2 / (B2 + W2)
      )
}
