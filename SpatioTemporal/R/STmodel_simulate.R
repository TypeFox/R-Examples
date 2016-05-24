#############################################
## S3-METHOD THAT SIMULATES FROM A STmodel ##
#############################################
##Functions in this file:
## simulate.STmodel     EX:ok

##' Data is simulated for the space-time locations in \code{object} using the
##' parameters in \code{x}.
##' 
##' @title Simulate Data from the Spatio-Temporal Model
##' 
##' @param object A \code{STmodel} object to perform unconditional
##'   simulation from.
##' @param nsim Number of replicates to simulate.
##' @param seed if !=NULL used in a call to \code{\link[base:set.seed]{set.seed}},
##'   allowing for replicatable simulation studies.
##' @param x Parameters to use when simulating the data; both regression and
##'   covariance parameters must be given, see \code{\link{loglikeSTgetPars}}.
##' @param nugget.unobs Value of nugget at unonserved locations, either a scalar
##'   or a vector with one element per unobserved site.
##' @param ... Additional parameters for \code{\link[base:set.seed]{set.seed}}
##' 
##' @return A list containing:
##'   \item{param}{Parameters used in the simulation, i.e. \code{x}.}
##'   \item{B}{The simulated beta fields in a (number of locations) -
##'            by - (number of temporal trends) - by -
##'            (number of replicates) array.}
##'   \item{X}{The simulated spatio-temporal fields in a
##'            (number of timepoints) - by - (number of locations) - by -
##'            (number of replicates) array. Row and column names indicate
##'            the time and locations for each point.}
##'   \item{obs}{A list with one element per replicate, containing the simulated
##'              observations extracted at space-time locations matching those in
##'              \code{object$obs}. To replace the observations with the i:th
##'              simulated values do:\cr
##'              \code{object$obs <- res$obs[[i]]}.}
##'
##' @example Rd_examples/Ex_simulate_STmodel.R
##' 
##' @author Johan Lindström
##' 
##' @family STmodel methods
##' @importFrom stats simulate
##' @importFrom MASS mvrnorm
##' @method simulate STmodel
##' @export
simulate.STmodel <- function(object, nsim=1, seed=NULL, x, nugget.unobs=0, ...){
  ##check class belonging
  stCheckClass(object, "STmodel", name="object")

  ##set random seed
  if( !is.null(seed) ){
    set.seed(seed, ...)
  }
  
  ##first figure out a bunch of dimensions
  dimensions <- loglikeSTdim(object)

  ##extract and check parameters
  if( length(x)==dimensions$nparam.cov ){
    ##compute alpha and gamma as cond. exp. given data
    tmp <- predict.STmodel(object, x, only.pars=TRUE, type="p")$pars
    x <- c(tmp$gamma.E,tmp$alpha.E, x)
  }
  if(length(x)!=dimensions$nparam){
    stop( sprintf("'x' should have %d or %d elements",
                  dimensions$nparam.cov, dimensions$nparam) )
  }
  ##add names to the parameters for future reference
  names(x) <- loglikeSTnames(object, TRUE)

  ##extract parameters from x
  tmp <- loglikeSTgetPars(x, object)
  gamma <- tmp$gamma
  alpha <- tmp$alpha
  cov.pars.beta <- tmp$cov.beta
  cov.pars.nu <- tmp$cov.nu
  
  ##nugget fo unobserved sites
  nugget.unobs <- internalFixNuggetUnobs(nugget.unobs, object,
                                         cov.pars.nu$nugget)
  ##get timepoints
  T <- object$trend$date
  
  ##create new time vectors F, one for each observation time.
  ##construct a matrix holding the temporal basis for each observation
  ##(including the constant)
  F <- matrix(1, length(T), dimensions$m)
  if( dimensions$m>1 ){
    I <- -which( names(object$trend)=="date" )
    F[, 2:dimensions$m] <- as.matrix(object$trend[, I, drop=FALSE])
  }
  colnames(F) <- rep("const", dim(F)[2])
  if( dimensions$m>1 ){
    colnames(F)[2:dimensions$m] <- names(object$trend)[I]
  }

  ##compute the new full distance matrices.
  D.nu <- crossDist( object$locations[, c("x.nu","y.nu"), drop=FALSE] )
  colnames(D.nu) <- rownames(D.nu) <- object$locations$ID
  D.beta <- crossDist( object$locations[, c("x.beta","y.beta"), drop=FALSE] )
  colnames(D.beta) <- rownames(D.beta) <- object$locations$ID

  ##calculate the land use regression for the temporal trends, as column
  mu.B <- matrix(calc.mu.B(object$LUR.all, alpha), ncol=1)
  
  ##create covariance matrices, beta-field
  sigma.B <- makeSigmaB(cov.pars.beta$pars, dist = D.beta,
                        type = object$cov.beta$covf,
                        nugget = cov.pars.beta$nugget)
  ##and nu-field
  ##just create residual matrix for all locations but one time point,
  ##Due to independence we just sample repeatedly from this.
  sigma.nu <- makeSigmaNu(cov.pars.nu$pars, dist = D.nu,
                          type = object$cov.nu$covf,
                          nugget = nugget.unobs,
                          random.effect = cov.pars.nu$random.effect,
                          blocks1 = dimensions$n, ind1 = 1:dimensions$n)
  
  ##calculate block cholesky factor of the matrices
  Rsigma.B <- try( makeCholBlock(sigma.B, n.blocks=dimensions$m), silent=TRUE)
  Rsigma.nu <- try( makeCholBlock(sigma.nu, n.blocks=1), silent=TRUE)

  ##array of simulated beta-fields
  B <- array(0, c(dimensions$n, dimensions$m, nsim))
  dimnames(B)[[1]] <- as.character( object$locations$ID )
  dimnames(B)[[2]] <- names( object$LUR.all )
  dimnames(B)[[3]] <- paste("nsim", as.character(1:nsim), sep=".")
  ##simulate data from B
  for(j in 1:dimensions$m){
    Ind <- (1:dimensions$n) + (j-1)*dimensions$n
    if( class(Rsigma.B)=="try-error" ){
      B[,j,] <-  t( mvrnorm(n=nsim, mu=mu.B[Ind], Sigma=sigma.B[Ind,Ind]) )
    }else{
      R <- t(Rsigma.B[Ind,Ind])
      B[,j,] <- matrix(mu.B[Ind], dimensions$n, nsim) + 
        R %*% matrix(rnorm(dim(R)[1]*nsim), dim(R)[1], nsim)
    }
  }
  ##create a time points by observations points matrix
  ##to hold the simulated data.
  X <- array(0, c(dim(F)[1], dimensions$n, nsim))
  dimnames(X)[[1]] <- as.character( T )
  dimnames(X)[[2]] <- object$locations$ID
  dimnames(X)[[3]] <- dimnames(B)[[3]]
  ##if needed add spatio-temporal covariate(s)
  if( dimensions$L!=0 ){
    ##extract relevant spatio-temporal covariates
    ##since this is STmodel we should be assured of consistency.
    for(i in 1:dim(object$ST.all)[3]){
      X[,,1] <- X[,,1] + gamma[i]*object$ST.all[,,i]
    }
    ##and replicate for the other X:s
    if(nsim>1){
      for(j in 2:nsim){
        X[,,j] <- X[,,1]
      }
    }
  }
    
  ##now combine simulated B, with the temporal trends and simulated
  ##data from the residuals field
  for(k in 1:nsim){
    X[,,k] <- X[,,k] + F %*% t(B[,,k])
  }
  ##simulate data from the nu fields
  if( class(Rsigma.nu)=="try-error" ){
    e <- t( mvrnorm(n=dim(F)[1]*nsim, mu=rep(0,dim(X)[2]), Sigma=sigma.nu) )
  }else{
    ##simulate from the residuals
    e <- rnorm( prod(dim(X)) )
    dim(e) <- c(dim(X)[2], dim(F)[1]*nsim)
    e <- t(Rsigma.nu) %*% e
  }
  e <- array(e, dim(X)[c(2,1,3)])
  e <- aperm(e, c(2,1,3))
  X <- X + e


  ##extract a vector matching the observed locations (i.e. matching the
  ##obs vector in data.obs)
  if( length(object$obs$obs)==0 ){
    obs = NULL
  }else{
    I <- ( (match(object$obs$ID, dimnames(X)[[2]])-1) * dim(X)[1] + 
          match(object$obs$date, convertCharToDate(dimnames(X)[[1]])) )
    obs <- list()
    for(i in 1:nsim){
      obs[[i]] <- object$obs
      tmp <- X[,,i]
      obs[[i]]$obs <- tmp[I]
      ##sort the observations (should be done since object$obs
      ##are sorted, but better safe than sorry)
      IND <- order(obs[[i]]$date, obs[[i]]$idx)
      obs[[i]] <- obs[[i]][IND,,drop=FALSE]
    }
  }
  out <- list(param=x, B=B, X=X, obs=obs)
  return( out )
}##simulate.STmodel
