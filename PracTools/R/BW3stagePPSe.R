
BW3stagePPSe <- function(dat, v, Ni, Qi, Qij, m){
    y <- dat[, v]
        # 3rd stage conditional wts
    wk.ij <- dat$w / dat$w2ij
        # 3rd stage sample counts
    qij <- table(dat$ssuID)
    qbbar <- mean(qij)

        # extract first row of dat for each PSU
    xx.psu <- do.call("rbind",list(by(1:nrow(dat),dat$psuID,head,1)))
    w1i.psu <- dat[xx.psu, "w1i"]
        # compute PSU 1-draw probs
        # assume that w1i = 1/(m*pp)
    pp <- 1/(m * dat[xx.psu,]$w1i)
        # extract first row of dat for each SSU
    xx.ssu <- do.call("rbind",list(by(1:nrow(dat),dat$ssuID,head,1)))
    ni <- table(dat[xx.ssu, "psuID"])
    f2i <- ni/Ni
    nbar <- mean(ni)
    w2ij.ssu <- dat[xx.ssu, "w2ij"]

        # SSU sel probs conditional on PSU selected
    w2ijC.ssu <- w2ij.ssu/dat[xx.ssu, "w1i"]
        # estimated no. of PSUs
    M.hat <- sum(w1i.psu)

    Qbbar <- sum(w2ij.ssu*Qij) / sum(w2ij.ssu)    # estimated mean no. of elements per SSU
    Qbar <- sum(w1i.psu*Qi) / M.hat               # estimated mean no. of elements per PSU

    tij <- by(wk.ij*y, dat$ssuID, sum)
    ti <- by(as.vector(w2ijC.ssu*tij), dat[xx.ssu,]$psuID, sum)
    t.pwr <- sum(dat$w * y)

    S2ai <- by(as.vector(tij), INDICES = dat[xx.ssu,]$psuID, var)
    S3ij <- by(y, INDICES = dat$ssuID, var)
    V3ij <- Qij * (Qij/qij - 1) * S3ij
    V3ijb <- Qij^2 * S3ij
    S2bi <- by(as.vector(V3ij), INDICES = dat[xx.ssu,]$psuID, sum)/ni
    sV3i <- by(as.vector(V3ij), INDICES = dat[xx.ssu,]$psuID, sum)
    sV3ib <- by(as.vector(V3ijb), INDICES = dat[xx.ssu,]$psuID, sum)

    S1a <- sum((ti/pp - t.pwr)^2)/(m-1)
    S1b <- sum(Ni^2/ni/pp^2*((1-f2i)*S2ai + f2i*S2bi))/m
    V3i <- by(y, INDICES = dat$psuID, wtdvar, w = dat$w)

    Vtsu <- sum(Ni^2/ni^2 * w1i.psu^2 * sV3i)
    Vssu <- sum(Ni^2/ni/(m*pp)^2*(1-f2i)*(S2ai - S2bi))
    Vpsu <- (S1a - S1b)/m

    B <-  (S1a - S1b) / t.pwr^2
    W <- sum(Qi^2 * V3i / m /pp^2) / t.pwr^2

    W2 <- sum(Ni^2/m/pp^2*(S2ai - S2bi))
    W2 <- W2 / t.pwr^2
    W3 <- sum(Ni^2/ni/m/pp^2 * sV3ib)
    W3 <- W3 / t.pwr^2

    delta1 <- B / (B + W)
    delta2 <- W2 / (W2 + W3)

    c(Vpsu=Vpsu, Vssu=Vssu, Vtsu=Vtsu,
      B=B, W=W, W2=W2, W3=W3,
      delta1=delta1, delta2=delta2)
}
