"dccm.nma" <-
  function(x, nmodes=NULL, ncore=NULL, ...) {
    nma <- x
    if (missing(nma))
      stop("dccm.nma: must supply a 'nma' object, i.e. from 'nma'")
    if(!"nma" %in% class(nma))
      stop("dccm.nma: must supply 'nma' object, i.e. from 'nma'")

    ## Check for multiple cores
    ncore <- setup.ncore(ncore, bigmem = FALSE)
   
    if(ncore > 1) {
      mcparallel <- get("mcparallel", envir = getNamespace("parallel")) 
      mccollect <- get("mccollect", envir = getNamespace("parallel")) 
    }

    ## Inner product between all pairs of residues
    cross.inner.prod <- function(a, b) {
      mat <- apply(a, 1, "%*%", t(b))
      return(mat)
    }

    ## Calc initial correlations for a subset of modes
    corrmats <- function(r.inds, core.id, nma, corr.mat, freqs) {
      for ( i in r.inds ) {
        mode <- matrix(nma$U[,i], ncol=3, byrow=TRUE)
        corr.mat <- corr.mat + (cross.inner.prod(mode, mode) / (freqs[i]**2))
        
        if(core.id==1)
          setTxtProgressBar(pb, i)
      }
      return(corr.mat)
    }
    
    if(!is.null(nma$frequencies)) {
      freqs <- nma$frequencies
    }
    else {
      freqs <- nma$force.constants
    }
    
    if(is.null(nmodes))
      nmodes <- length(nma$L)
    else {
      nmodes <- nmodes + nma$triv.modes
      if(nmodes>length(nma$L)) {
        warning("'nmodes' larger than the number of modes")
        nmodes <- length(nma$L)
      }
    }

    ## Initialize progress bar
    ##ptm <- proc.time()
    pb <- txtProgressBar(min=(nma$triv.modes+1), max=nmodes + nma$natoms, style=3)

    ## Allocate the correl matrix
    corr.mat <- matrix(0, nma$natoms, nma$natoms)
    
    ## Which modes to use for calculation
    mode.inds <- (nma$triv.modes+1):nmodes
    core.ids  <- rep(1:ncore, length.out=length( mode.inds ))

    if(ncore>1) 
      jobs <- list()
    
    for ( i in 1:ncore ) {
      rinds <- mode.inds[ which(core.ids==i) ]

      if(ncore>1) {
        q <- mcparallel(corrmats(rinds, i, nma, corr.mat, freqs))
        jobs[[i]] <- q
      }
      else
        corr.mat <- corrmats(rinds, i, nma, corr.mat, freqs)
    }

    ## Collect all jobs, and sum matrices
    if(ncore>1) {
      res <- mccollect(jobs, wait=TRUE)
      for ( job in res ) {
        corr.mat <- corr.mat + job
      }
    }
  
    ## Basis for normalization
    a <- vector('numeric', length=nrow(corr.mat))
    k <- length(mode.inds) ## for ProgressBar !
    inds <- rep(1:nrow(corr.mat), each=3)
    
    for ( j in (nma$triv.modes+1):nmodes )  {
      v <- nma$U[, j] * nma$U[, j]
      a <- a + ( tapply( v, inds, sum) / (freqs[j]**2))
      
      k <- k+1
      setTxtProgressBar(pb, k)
    }
    
    close(pb)
    
    a <- sqrt(a)
    bn <- a%o%a
    
    ## Normalized correlation matrix
    corr.mat <- corr.mat / bn
    class(corr.mat) <- c("dccm", "matrix")
    
    ##t <- proc.time() - ptm
    ##cat(" Done in", t[[3]], "seconds.\n")
    
    return(corr.mat)
  }
