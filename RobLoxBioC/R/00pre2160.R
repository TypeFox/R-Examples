## due to a change to .C and .Call in 2.16.0
.getA1.locsc <- RobLox:::.getA1.locsc
.getA2.locsc <- RobLox:::.getA2.locsc
.geta.locsc <- RobLox:::.geta.locsc
.getb.locsc <- RobLox:::.getb.locsc

.getlsInterval <- function (r, rlo, rup){
    if (r > 10) {
        b <- 1.618128043
        const <- 1.263094656
        A2 <- b^2 * (1 + r^2)/(1 + const)
        A1 <- const * A2
    }
    else {
        A1 <- .getA1.locsc(r)
        A2 <- .getA2.locsc(r)
        b <- .getb.locsc(r)
    }
    if (rlo == 0) {
        efflo <- (A1 + A2 - b^2 * r^2)/1.5
    }
    else {
        A1lo <- .getA1.locsc(rlo)
        A2lo <- .getA2.locsc(rlo)
        efflo <- (A1 + A2 - b^2 * (r^2 - rlo^2))/(A1lo + A2lo)
    }
    if (rup > 10) {
        bup <- 1.618128043
        const.up <- 1.263094656
        A2up <- bup^2 * (1 + rup^2)/(1 + const.up)
        A1up <- const.up * A2up
        effup <- (A1 + A2 - b^2 * (r^2 - rup^2))/(A1up + A2up)
    }
    else {
        A1up <- .getA1.locsc(rup)
        A2up <- .getA2.locsc(rup)
        effup <- (A1 + A2 - b^2 * (r^2 - rup^2))/(A1up + A2up)
    }
    return(effup - efflo)
}

.onestep.locsc.matrix <- RobLox:::.onestep.locsc.matrix
.kstep.locsc.matrix <- RobLox:::.kstep.locsc.matrix