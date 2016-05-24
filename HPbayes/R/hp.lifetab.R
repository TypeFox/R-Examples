hp.lifetab <-
function(hpp, nax, age = seq(0, 85, 1), l0 = 1e+05, with.CI=FALSE, CI=95) 
{
    lt <- NULL
    lt.lo <- NULL
    lt.hi <- NULL
    
    H.new.hat <- hp.nqx(H.out = hpp, age = age)
    nqx <- rep(NA, length(age[-1]))
    nage <- length(age)
    for (i in 1:length(nqx)) {
        nqx[i] <- median(H.new.hat[, i])
    }
    nqx <- c(nqx, 1)
    n <- c(diff(age), 999)
    npx <- 1 - nqx
    lx <- round(cumprod(c(l0, 1 - nqx)))
    ndx <- -diff(lx)
    lxpn <- lx[-1]
    nLx <- n * lxpn + ndx * nax
    Tx <- rev(cumsum(rev(nLx)))
    lx <- lx[1:length(age)]
    ex <- Tx/lx
    lt <- cbind(Age = age, nax = nax, nqx = round(nqx, 4), npx = round(npx, 
        4), ndx = ndx, lx = lx, nLx = round(nLx), Tx = round(Tx), 
        ex = round(ex, 2))

    if(with.CI==TRUE){
    loCI <- ((100-CI)/2)/100
	hiCI <- 1-(((100-CI)/2)/100)
	
    nqx.lo <- rep(NA, length(age[-1]))
    nage <- length(age)
    for (i in 1:length(nqx.lo)) {
        nqx.lo[i] <- quantile(H.new.hat[, i], probs=loCI)
    }
    nqx.lo <- c(nqx.lo, 1)
    n <- c(diff(age), 999)
    npx.lo <- 1 - nqx.lo
    lx.lo <- round(cumprod(c(l0, 1 - nqx.lo)))
    ndx.lo <- -diff(lx.lo)
    lxpn.lo <- lx.lo[-1]
    nLx.lo <- n * lxpn.lo + ndx.lo * nax
    Tx.lo <- rev(cumsum(rev(nLx.lo)))
    lx.lo <- lx.lo[1:length(age)]
    ex.lo <- Tx.lo/lx.lo
    lt.lo <- cbind(Age = age, nax = nax, nqx.lo = round(nqx.lo, 4), npx.lo = round(npx.lo, 
        4), ndx.lo = ndx.lo, lx.lo = lx.lo, nLx.lo = round(nLx.lo), Tx.lo = round(Tx.lo), 
        ex.lo = round(ex.lo, 2))

    nqx.hi <- rep(NA, length(age[-1]))
    nage <- length(age)
    for (i in 1:length(nqx.hi)) {
        nqx.hi[i] <- quantile(H.new.hat[, i], probs=hiCI)
    }
    nqx.hi <- c(nqx.hi, 1)
    n <- c(diff(age), 999)
    npx.hi <- 1 - nqx.hi
    lx.hi <- round(cumprod(c(l0, 1 - nqx.hi)))
    ndx.hi <- -diff(lx.hi)
    lxpn.hi <- lx.hi[-1]
    nLx.hi <- n * lxpn.hi + ndx.hi * nax
    Tx.hi <- rev(cumsum(rev(nLx.hi)))
    lx.hi <- lx.hi[1:length(age)]
    ex.hi <- Tx.hi/lx.hi
    lt.hi <- cbind(Age = age, nax = nax, nqx.hi = round(nqx.hi, 4), npx.hi = round(npx.hi, 
        4), ndx.hi = ndx.hi, lx.hi = lx.hi, nLx.hi = round(nLx.hi), Tx.hi = round(Tx.hi), 
        ex.hi = round(ex.hi, 2))

    }
    return(list(lt = lt, lt.lo = lt.lo, lt.hi = lt.hi))
}

