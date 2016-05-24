## =======================================================
## (Adaptive) Metropolis Sampler
##
## Implementation of the RAM (robust adaptive Metropolis)
## sampler of
## Vihola, M. (2011) Robust adaptive Metropolis algorithm with
## coerced acceptance rate. Statistics and Computing.
## [online] http://www.springerlink.com/content/672270222w79h431/
## (Accessed December 8, 2011).

## Version 1.1

## June 22, 2012 -- Andreas Scheidegger
## =======================================================


MCMC <- function(p, n, init, scale=rep(1, length(init)),
                 adapt=!is.null(acc.rate), acc.rate=NULL, gamma=0.5, list=TRUE, n.start=0, ...) {

  ## checks
  if(adapt & !is.numeric(acc.rate)) stop('Argument "acc.rate" is missing!')
  if(gamma<0.5 | gamma>1) stop('Argument "gamma" must be between 0.5 and 1!')

  ## number of adaption steps
  if(is.numeric(adapt)) n.adapt <- adapt
  if(adapt==TRUE) n.adapt <- Inf
  if(adapt==FALSE) n.adapt <- 0

  ## number of parameter
  d <- length(init)

  ## matrix to store MC chain
  X <- matrix(NA, ncol=d, nrow=n)
  colnames(X) <- names(init)
  X[1,] <- init

  ## vector to store log densities p(x)
  p.val <- rep(NA, n)
  p.val[1] <- p(X[1,], ...)

  ## initial S
  if(length(scale)==d) {
    M <- diag(scale)
  } else {
    M <- scale
  }
  ## check
  if(ncol(M) != length(init)) stop("Length or dimension of 'init' and 'scale' do not match!")
  if(length(init)==1) stop('One-dimensional sampling is not possible!')

  S <-  t(chol(M))

  ## initialize progress bar
  cat('  generate', n, 'samples \n')
  pb <- txtProgressBar(min=0, max=n, style=3)
  update.step <- max(5, floor(n/100))

  k <- 0
  for(i in 2:n){

    if(i %% update.step == 0) {
      setTxtProgressBar(pb, i)
    }

    ## proposal value
    U <- rnorm(d)
    X.prop <- c( X[i-1,] + S %*% U )
    names(X.prop) <- names(init)

    ## calculate density at X.prop
    p.val.prop <- p(X.prop, ...)

    ## acceptance probability
    alpha <- min(1, exp( p.val.prop - p.val[i-1] )) # for log density

    if(!is.finite(alpha)) alpha <- 0    # if zero divided by zero

    ## accept with P=alpha
    if(runif(1)<alpha) {
      X[i,] <- X.prop                   # accept
      p.val[i] <- p.val.prop
      k <- k+1
    } else {
      X[i,] <- X[i-1,]                  # or not
      p.val[i] <- p.val[i-1]
    }

    ## compute new S
    ii <- i+n.start
    if(ii < n.adapt) {
      adapt.rate <-  min(5, d*ii^(-gamma))
      M <- S %*% (diag(d) + adapt.rate*(alpha - acc.rate) * U%*%t(U)/sum(U^2)) %*% t(S)

      ## check if M is positive definite. If not, use nearPD().
      eig <- eigen(M, only.values = TRUE)$values
      tol <- ncol(M)*max(abs(eig))*.Machine$double.eps

      if( !isSymmetric(M) | is.complex(eig) | !all(Re(eig)>tol) ){
        ## nearPD() computes the 'nearest' positive definite matrix
        M <- as.matrix(Matrix::nearPD(M)$mat)
      }

      S <- t(chol(M))

    }
  }

  close(pb)                             #close progress bar

  ## calculate accpetance rate
  acceptance.rate <- round(k/(n-1), 3)

  if(list) {
    return(list(samples=X,
                log.p=p.val,
                cov.jump=M,
                n.sample=n,
                acceptance.rate=acceptance.rate,
                adaption=adapt,
                sampling.parameters=list(sample.density=p,
                  acc.rate=acc.rate,
                  gamma=gamma)
                ) )
  } else {
    cat("Acceptance rate:", acceptance.rate, "\n")
    return(X)
  }
}


## ----------------------
## Adds more samples to an existing chain
## does work with objet generated with MCMC.parallel()
## but updating is not computed  parallel.

MCMC.add.samples <- function(MCMC.object, n.update, ...) {

    ## if single chain
    if(!is.null(names(MCMC.object))) {
        if(is.matrix(MCMC.object)) stop("Only MCMC objects generated with option 'list=TRUE' can be updated!")

        ## get values from last sampling
        p <- MCMC.object$sampling.parameters$sample.density
        init <- MCMC.object$samples[nrow(MCMC.object$samples),]
        scale <- MCMC.object$cov.jump
        acc.rate <- MCMC.object$sampling.parameters$acc.rate
        gamma <- MCMC.object$sampling.parameters$gamma
        n.before <- MCMC.object$n.sample    # number of existing samples
        adapt <- MCMC.object$adaption

        ## generate new samples
        samp.update <- MCMC(p=p, n=n.update, init=init, scale=scale,  adapt=adapt, acc.rate=acc.rate,
                            gamma=gamma, list=TRUE, n.start=n.before, ...)

        ## update old sampling object
        MCMC.object$cov.jump <- samp.update$cov.jump
        m <- c(MCMC.object$n.sample, samp.update$n.sample)
        MCMC.object$acceptance.rate <-  1/sum(m)*(m[1]*MCMC.object$acceptance.rate + m[2]*samp.update$acceptance.rate)
        MCMC.object$n.sample <- MCMC.object$n.sample + n.update

        MCMC.object$samples <- rbind(MCMC.object$samples, samp.update$samples)
        MCMC.object$log.p <- c(MCMC.object$log.p, samp.update$log.p)

        ## return the updated list
        return(MCMC.object)
    }

    ## if list of chains
    if(is.null(names(MCMC.object))) {
        ## recursive call of MCMC.add.samples() to update single chains
        MCMC.object <- lapply(MCMC.object, function(x) MCMC.add.samples(x, n.update=n.update, ...))
        return(MCMC.object)
    }
}

## ----------------------
## Wrapper for parallel calculation of independent chains
MCMC.parallel <- function(p, n, init, n.chain=4, n.cpu, packages=NULL, dyn.libs=NULL,
                          scale=rep(1, length(init)), adapt=!is.null(acc.rate),
                          acc.rate=NULL, gamma=0.55, list=TRUE, ...) {

  require(parallel)

  ## initialisation of (local) parallel computing
  cl <- makeCluster(min(n.cpu, detectCores()))

  ## stop parallel computing on exit
  on.exit({ stopCluster(cl); print("Cluster stopped.")})

  ## export complete work space of master
  varlist <- unique(c(ls(), ls(envir=.GlobalEnv), ls(envir=parent.env(environment()))))
  clusterExport(cl, varlist=varlist, envir=environment())

  ## init random generators
  clusterSetRNGStream(cl)

  ## export 'packages', 'dyn.libs', and current working directory
  wd <- getwd()
  clusterExport(cl, varlist=c("packages", "dyn.libs", "wd"), envir=environment())

  ## wrapper function to be called in parallel
  MCMC.wrap <- function(x, ...) {
      require(Matrix)
      if(!is.null(packages)) sapply(packages, function(x) require(x, character.only=TRUE))
      ## load and unload dynamic libraries
      if (!is.null(dyn.libs)) {
          sapply(dyn.libs, function(x) dyn.load(paste(wd, x, sep = "/")))
          on.exit( sapply(dyn.libs, function(x) dyn.unload(paste(wd, x, sep = "/"))) )
      }
      MCMC(...)
  }


  ## sample n chains in parallel
  result <- clusterApply(cl, 1:n.chain, MCMC.wrap, p=p, n=n, init=init,
                        scale=scale,  adapt=adapt, acc.rate=acc.rate,
                        gamma=gamma, list=list, ...)

  return(result)

}


## ----------------------
## converts a sample into coda object

convert.to.coda <- function(sample) {
    ## if single chain
    if(!is.null(names(sample))) {
        if(is.matrix(sample)) {
            obj <- coda::mcmc(sample)
        }
        if(is.list(sample)) {
            obj <- coda::mcmc(sample$samples)
        }
        return(obj)
    } else {

        ## if sample is a list of chains
        if(is.matrix(sample[[1]])) {
            obj <- as.mcmc.list(lapply(sample, coda::mcmc))
        }
        if(is.list(sample[[1]])) {
            obj <- as.mcmc.list(lapply(sample, function(x) {coda::mcmc(x$samples)}))
        }
        return(obj)
    }
}


