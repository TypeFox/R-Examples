
BW3stagePPS <- function(X, pp, psuID, ssuID){
    tUi <- by(X, INDICES = psuID, FUN = sum)
    tU <- sum(tUi)

    S2U1pwr <- sum(pp * (tUi/pp - tU)^2)
    B <- S2U1pwr / tU^2

        # extract first ssu in each psu and count ssu's per psu
        # identify psu IDs for 1st instance of each ssuID
    first <- do.call("rbind",list(by(psuID,ssuID,head,1)))
        # identify  first row for each ssu
        xx <- do.call("rbind",list(by(1:length(X),ssuID,head,1)))

    Ni <- table(first)
    N <- sum(Ni)
    M <- length(unique(psuID))
    Nbar <- mean(Ni)

        # ssu totals and their vars within psu
    tij <- by(X, INDICES = ssuID, FUN = sum)
    S2U2i <- by(as.vector(tij), INDICES = psuID[xx], FUN = var)
        # identify any PSUs with only 1 SSU
    not.na <- !is.na(S2U2i)
    W2 <- sum(Ni[not.na]^2 * S2U2i[not.na]/pp[not.na]) / tU^2

    Q <- length(X)
    Qij <- table(ssuID)
    Qbbar <- mean(Qij)
    Qi <- table(psuID)
    Qbar <- mean(Qi)
    S2U3i <- by(X, INDICES = psuID, FUN = var)
    S2U3ij <- by(X, INDICES = ssuID, FUN = var)

        # compute W needed for delta1
    W <- sum(Qi^2 * S2U3i/pp) / (tU)^2

    ppr <- rep(pp, Ni)
    Nir <- rep(Ni, Ni)
    W3 <- sum(Nir * Qij^2 * S2U3ij/ppr) / tU^2

    V <- var(X)/mean(X)^2

    delta1 <- B / (B + W)
    delta2 <- W2 / (W2 + W3)

    c("B"=B, "W"=W, "W2"=W2, "W3"=W3,
        "unit relvar"=V,
        "k1"=(B+W)/V, "k2"=(W2+W3)/V,
        "delta1"=delta1, "delta2"=delta2)
}
