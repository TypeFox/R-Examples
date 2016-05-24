cudaMultireg.slice <-
function (slicedata, ymaskdata, R=3000, keep=5, nu.e=3, fsave=NA, zprior=FALSE, rng=0) 
{
    X <- slicedata$X
    nvar <- slicedata$nvar
    nobs <- slicedata$nobs
    yn <- ymaskdata$yn # y data, masked and standardized
    nreg <- ymaskdata$nreg
    stopifnot(nobs == nrow(yn))
    #-----------------------------
    # Mcmc regression stuff
    if (!zprior) # Default prior for Z: simplest case
        Z <- matrix(rep(1, nreg), ncol = 1)
    else  # Z prior as CSF/GRY/WHT segmented regions
		    Z <- read.Zsegslice(slicedata = slicedata, ymaskdata = ymaskdata)
    nz <- ncol(Z)
    # Allocate space for the draws and set initial values of Vbeta and Delta
    Deltabar <- matrix(0, nz, nvar)
    Vbetadraw = matrix(double(floor(R/keep) * nvar * nvar), ncol = nvar * 
        nvar)
    Deltadraw = matrix(double(floor(R/keep) * nz * nvar), ncol = nz * 
        nvar)
    taudraw = matrix(double(floor(R/keep) * nreg), ncol = nreg)
    betadraw = array(double(floor(R/keep) * nreg * nvar), dim = c(nreg, 
        nvar, floor(R/keep)))
    #-----------------------------
    runif(1)
    cat("\nBegin of Cuda call \n")
    outcuda <- .C("cudaMultireg", as.single(yn), as.single(X), 
        as.single(Z), as.single(Deltabar), as.integer(nz), as.integer(nu.e), 
        as.integer(nreg), as.integer(nobs), as.integer(nvar),
        as.integer(R), as.integer(keep), as.integer(rng),
        Vbetadraw = as.single(Vbetadraw), Deltadraw = as.single(Deltadraw),
        betadraw = as.single(betadraw), taudraw = as.single(taudraw))
    cat("\nEND of Cuda call \n")
    #-----------------------------
    Vbetadraw <- matrix(outcuda$Vbetadraw, ncol = nvar * nvar, 
        byrow = T)
    Deltadraw <- matrix(outcuda$Deltadraw, ncol = nz * nvar, 
        byrow = T)
    betadraw <- array(outcuda$betadraw, dim = c(nreg, nvar, floor(R/keep)))
    taudraw <- matrix(outcuda$taudraw, ncol = nreg, byrow = T)
    #-----------------------------
		# Apply "bayesm" methods used for summary compatibility purposes
    attributes(taudraw)$class = c("bayesm.mat", "mcmc")
    attributes(taudraw)$mcpar = c(1, R, keep)
    attributes(Deltadraw)$class = c("bayesm.mat", "mcmc")
    attributes(Deltadraw)$mcpar = c(1, R, keep)
    attributes(Vbetadraw)$class = c("bayesm.var", "bayesm.mat", 
        "mcmc")
    attributes(Vbetadraw)$mcpar = c(1, R, keep)
    attributes(betadraw)$class = c("hcoef.post") # attributes for new plot method 
    out <- list(Vbetadraw = Vbetadraw, Deltadraw = Deltadraw, 
        betadraw = betadraw, taudraw = taudraw)
    if(!is.na(fsave)) {
        cat("saving simulation ", fsave, "...")
        save(out, file = fsave)
        cat("\n")
    }
    invisible(out)
}
