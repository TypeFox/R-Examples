## Function which performs importance sampling using the
## adaptive mixture of Student-t densities as the importance density
## __input__
## N       : [integer>0] number of draws used in the importance sampling (default: 1e5)
## KERNEL  : [function] which computes the kernel !! must be vectorized !!
## G       : [function] used in the importance sampling !! must be vectorized !! (default: NULL, i.e. posterior mean)
## mit     : [list] with mixture information
## ...     : additional parameters used by 'KERNEL' or/and by 'G'
## __output__
## [list] with the following components:
## $ghat   : [vector] importance sampling estimator
## $NSE    : [vector] numerical standard error estimated as in Geweke (1989, p.1321)
## $RNE    : [vector] relative numerical efficiency estimated as in Geweke (1989, p.1321)
## __20080521__
'AdMitIS' <- function(N=1e5, KERNEL, G=function(theta){theta}, mit=list(), ...)
  {
    if (N<1)
      stop ("'N' should be larger than 1")
    if (missing(KERNEL))
      stop ("'KERNEL' is missing in 'AdMitIS'")
    KERNEL <- match.fun(KERNEL)
    if (!any(names(formals(KERNEL))=="log"))
      stop ("'KERNEL' MUST have the logical argument 'log' which returns the (natural) logarithm of 'KERNEL' if 'log=TRUE'")      
    G <- match.fun(G)
    theta <- rMit(N, mit)
    ## arguments in list(...)
    args <- list(...)
    nargs <- names(args)
    ## arguments for KERNEL
    argsK <- formals(KERNEL)
    nargsK <- names(argsK)
    ## arguments for do.call('fn.w')
    argsw <- list(theta)
    argsw <- c(argsw, list(KERNEL=KERNEL, mit=mit, log=FALSE))
    nargsargsK <- nargsK[charmatch(nargs, nargsK, nomatch=0)]
    argsw <- c(argsw, args[nargsargsK])
    w <- do.call('fn.w', argsw)
    ## arguments for G
    argsG <- formals(G)
    nargsG <- names(argsG)
    ## arguments for do.call('G')
    argsg <- list(theta)
    nargsargsG <- nargsG[charmatch(nargs, nargsG, nomatch=0)]
    argsg <- c(argsg, args[nargsargsG])
    g <- do.call('G', argsg)    

    k <- ncol(g)
    w <- w/sum(w)
    ghat <- apply(w*g, 2, sum)
    tmp <- g-matrix(ghat, N, k, byrow=TRUE)
    NSE <- sqrt( apply(w^2*tmp^2, 2, sum) )
    RNE <- (apply( w*tmp^2, 2, sum) / N) / NSE^2 ## naive variance / NSE^2
    
    list(ghat=as.vector(ghat),
         NSE=as.vector(NSE),
         RNE=as.vector(RNE))
  }
