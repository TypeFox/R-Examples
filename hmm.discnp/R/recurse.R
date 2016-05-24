recurse <- function(fy,tpm,ispd,lns)
{
#
# Function recurse to calculate the ``recursive probabilities'',
# given the parameters theta, and the observations y.
#

# Set a bunch of constants:
K  <- nrow(tpm)
K2 <- K*K
L  <- ncol(fy)
M  <- K*L
nreps <- length(lns)
nxi   <- L - nreps
N <- K*K*nxi
epsilon <- sqrt(.Machine$double.eps)
if(is.matrix(ispd)) {
    if(ncol(ispd) != nreps)
        stop("Number of columns of \"ispd\" must equal \"nreps\".\n")
    cis <- FALSE
    nis <- nreps
} else {
    cis <- TRUE
    nis <- 1
}

# Recursive probabilities:

        rp <- .Fortran(
                'recurse',
		fy=as.double(fy),
                xispd=as.double(ispd),
                tpm=as.double(tpm),
		nreps=as.integer(nreps),
                epsilon=as.double(epsilon),
		lns=as.integer(lns),
                nstate=as.integer(K),
                nis=as.integer(nis),
                cis=as.logical(cis),
                wrk=double(K2),
		xlc=double(L),
                ntot=as.integer(L),
                nxi=as.integer(nxi),
                alpha=double(M),
                beta=double(M),
                gamma=double(M),
                xi=double(N),
		xisum=double(K2),
                PACKAGE="hmm.discnp"
        )

list(gamma=matrix(rp$gamma,nrow=K),xi=matrix(rp$xisum,nrow=K),llc=rp$xlc)
}
