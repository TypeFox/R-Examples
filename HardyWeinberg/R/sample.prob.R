sample.prob <- function (n, nm, mA, nA, nfAB, verbose=FALSE) 
{
    # conditioning on the number of A alleles and the number of males.
    nf <- n - nm
    nt <- nm + 2*nf
    nB <- nt - nA
    mB <- nm - mA
    nfAA <- 0.5 * (nA - mA - nfAB)
    nfBB <- 0.5 * (nB - mB - nfAB)
    if(nfBB < 0) {
      cat("begin\n")
      cat(n,nm,nf,"\n")
      cat(nA,nB,nt,"\n")
      cat(mA,mB,"\n")
      cat(nfAA,nfAB,nfBB,"\n")
      stop("here")
      
    }
    if(verbose) {
        cat(n,nm,nf,"\n")
        cat(nt,nA,nB,"\n")
        cat(nm,mA,mB,"\n")
        cat(nfAA,nfAB,nfBB,"\n")
    }
    numer <- lfactorial(nm) + lfactorial(nf) + lfactorial(nA) + lfactorial(nB) + nfAB * log(2) 
    denom <- lfactorial(mA) + lfactorial(mB) + lfactorial(nfAA) + lfactorial(nfAB) + lfactorial(nfBB) + lfactorial(nt)
    prob <- exp(numer - denom)
    return(prob)
}
