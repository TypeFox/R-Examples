pbounds <- function(h=0, pT, iD=list(T=0)) {
    K<-length(pT$t)
    if(iD$T >= K) {
        print("Error: Cannot compute conditional bounds for T>=K")
        rep(NaN, K)
    } else {
        -h * sqrt(pT$t[(iD$T+1):K] * pT$Imax - ifelse(iD$T == 0, 0, pT$t[iD$T] * pT$Imax)) +
            (pT$b[(iD$T+1):K] * sqrt(pT$t[(iD$T+1):K] * pT$Imax) - ifelse(iD$T == 0, 0, iD$z * sqrt(pT$t[iD$T] * pT$Imax))) /
                sqrt(pT$t[(iD$T+1):K] * pT$Imax - ifelse(iD$T == 0, 0, pT$t[iD$T] * pT$Imax))
    }
}

