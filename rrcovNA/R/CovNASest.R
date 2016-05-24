CovNASest <- function(x,
                    bdp=0.5,
                    arp=0.1,
                    eps=1e-5,
                    maxiter=120,
                    nsamp=500,
                    impMeth=c("norm" , "seq", "rseq"),
                    seed=NULL,
                    trace=FALSE,
                    tolSolve=10e-14,
                    scalefn,
                    method=c("sfast", "surreal", "bisquare", "rocke", "suser", "sdet"),
                    control,
                    t0,
                    S0,
                    initcontrol
                )
{
    impMeth <- match.arg(impMeth)
    method <- match.arg(method)
    if(!missing(control)){
        defcontrol <- CovControlSest()     # default control
        if(bdp == defcontrol@bdp)           bdp <- control@bdp          # for s-fast and surreal
        if(arp == defcontrol@arp)           arp <- control@arp          # for rocke type
        if(eps == defcontrol@eps)           eps <- control@eps          # for bisquare and rocke
        if(maxiter == defcontrol@maxiter)   maxiter <- control@maxiter  # for bisquare and rocke

        if(nsamp == defcontrol@nsamp)       nsamp <- control@nsamp
        if(is.null(seed) || seed == defcontrol@seed)         seed <- control@seed
        if(trace == defcontrol@trace)       trace <- control@trace
        if(tolSolve == defcontrol@tolSolve) tolSolve <- control@tolSolve
        if(method == defcontrol@method)     method <- control@method
    }

    if(length(seed) > 0) {
        if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))  {
            seed.keep <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir=.GlobalEnv))
        }
        assign(".Random.seed", seed, envir=.GlobalEnv)
    }

    if(is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
                    dimnames = list(names(x), deparse(substitute(x))))

    xcall <- match.call()
    n <- nrow(x)
    p <- ncol(x)

    ## drop all rows which contain only missings
    na.x <- rowSums(ifelse(is.na(x),1,0)) == ncol(x)
    ok <- !na.x
    x <- x[ok, , drop = FALSE]

    dimn <- dimnames(x)
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]
    if(p < 2)
        stop("Need at least 2 columns ")

        ximp <- .imputation(x, impMeth = impMeth)
        sss <- CovSest(ximp,
                        bdp=bdp,
                        arp=arp,
                        eps=eps,
                        maxiter=maxiter,
                        nsamp=nsamp,
                        seed=seed,
                        trace=trace,
                        tolSolve=tolSolve,
                        scalefn=scalefn,
                        method=method,
                        control=control,
                        t0=t0,
                        S0=S0,
                        initcontrol=initcontrol)

        ans <- new("CovNASest",
               call = xcall,
               iter=sss@iter,
               crit=sss@crit,
               cov=sss@cov,
               center=sss@center,
               n.obs=n,
               X = as.matrix(x),
               method=sss@method)
    ans
}
