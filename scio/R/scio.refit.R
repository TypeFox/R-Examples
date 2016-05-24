scio.refit <- function(S, Omega, thr = 1e-04, ...) {
    if(!require(QUIC))  {
        warning("Refitting did not run because package QUIC is not available!")
        return(list(w=NULL))
    }
    HUGE <- 1e8
    ss <- dim(Omega)

    if (length(ss)>2) {
        w <- 0*Omega 
        for (jj in 1:ss[3]) {
            rho <- HUGE*(abs(Omega[,,jj])<thr) + thr
            w[,,jj] <-  QUIC(S, rho, ...)$X
        }
    } else {
        rho <- HUGE*(abs(Omega)<thr) + thr
        w <- QUIC(S, rho, ...)$X
    }
    return(list(w=w))
}
