wandafromx <-
function (x) 
{
## note: this is a temporary copy of the wandafromx function from EbayesThresh
# modified to give a list output instead of multiargument returns (accepted in
# old versions of R). [cf R CMD check issues in binhf examples].

    thi <- sqrt(2 * log(length(x)))
    lo <- c(0, 0.04)
    hi <- c(thi, 3)
    startpar <- c(1, 0.5)
    if (exists("optim")) {
        uu <- optim(startpar, negloglik.laplace, method = "L-BFGS-B", 
            lower = lo, upper = hi, xx = x)
        uu <- uu$par
    }
    else {
        uu <- nlminb(startpar, negloglik.laplace, lower = lo, 
            upper = hi, xx = x)
        uu <- uu$parameters
    }
    a <- uu[2]
    w <- wfromt(uu[1], a = a)
    return(list(w=w, a=a))
}
