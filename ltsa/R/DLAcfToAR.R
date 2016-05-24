`DLAcfToAR` <-
function(r, useC=TRUE, PDSequenceTestQ=FALSE)
{
    if (abs(r[1]) >= 1)
        stop(paste("error - lag one autocorrelation, r[1] =",r[1],", must be <1"))
    R <- c(1,r)
    np1 <- length(R)
    n <- np1-1
    EPS <- .Machine$double.eps # 1+EPS==1, machine epsilon
    if (useC)    {
        out<-.C("DLar", as.double(R), phi=as.double(numeric(length(R))), zta=as.double(numeric(length(R))),
                      sigmasq=as.double(numeric(length(R))), as.integer(length(R)), 
                      as.double(EPS), fault = as.integer(1),
                      PACKAGE="ltsa" )
        fault<-out$fault
        if (fault == 1)
                if  (PDSequenceTestQ)
                    return(FALSE)
                else
                    stop("error: r is not a p.d. sequence")
        phi<-(out$phi)[1:n]
        zta<-(out$zta)[1:n]
        sigmasq<-(out$sigmasq)[-1]
        ans <- matrix( c(phi,zta,sigmasq), ncol=3, nrow=n)
    }
    else {
        zta <- sigmasq <- numeric(np1)
        sigmasq[1] <- R[1]
        zta[2]<- phi <- R[2]/R[1]
        sigmasqkm1 <- R[1] * (1 - phi^2)
        sigmasq[2] <- sigmasqkm1
        if (sigmasqkm1<EPS)
            if (PDSequenceTestQ)
                return(FALSE)
            else
                stop("error: r is not a p.d. sequence")
        if (n > 1) 
            for(k in 2:n) {
                zta[k+1] <- phikk <- (R[k + 1] - phi %*% rev(R[2:k]))/sigmasqkm1
                sigmasqk <- sigmasqkm1 * (1 - phikk^2)
                phinew <- phi - phikk * rev(phi)
                phi <- c(phinew, phikk)
                sigmasqkm1 <- sigmasqk
                sigmasq[k + 1] <- sigmasqk
                if (sigmasqk < EPS) {
                    if (PDSequenceTestQ)
                        return(FALSE)
                else
                    stop("error: r is not a p.d. sequence")
                }
            }
        ans <- matrix( c(phi,zta[-1],sigmasq[-1]), ncol=3)
        }
    if (PDSequenceTestQ)
                    ans<-TRUE
                else 
                    dimnames(ans) <- list(1:length(phi), c("phi", "phikk", "sigsqk"))
    ans
}

