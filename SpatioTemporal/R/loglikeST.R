#################################################
## FILE CONTAINING THE LOGLIKELIHOOD FUNCTIONS ##
#################################################
##Functions in this file:
## loglikeST        EX:ok
## loglikeSTnaive   EX:with loglikeST

##' Computes the log-likelihood for the spatio-temporal model.  \code{loglikeST}
##' uses an optimised version of the log-likelihood, while \code{loglikeSTnaive}
##' uses the naive (slow) version and is included mainly for testing and speed
##' checks.
##' 
##' @title Compute the Log-likelihood for the Spatio-Temporal Model
##' @param x Point at which to compute the log-likelihood, should be only
##'   \emph{log}-covariance parameters if \code{type=c("p","r")} and
##'   regression parameters followed by \emph{log}-covariance parameters if
##'   \code{type="f"}. If \code{x=NULL} the function acts as an alias for
##'   \code{\link{loglikeSTnames}} returning the expected names of the
##'   input parameters.
##' @param STmodel \code{STmodel} object with the model for which to compute
##'   the log-likelihood.
##' @param type A single character indicating the type of log-likelihood to
##'   compute. Valid options are "f", "p", and "r", for \emph{full},
##'   \emph{profile} or \emph{restricted maximum likelihood} (REML).
##' @param x.fixed Parameters to keep fixed, \code{NA} values in this vector is
##'   replaced by values from \code{x} and the result is used as \code{x}, ie. \cr
##'   \code{ x.fixed[ is.na(x.fixed) ] <- x} \cr \code{ x <- x.fixed }.
##' 
##' @return Returns the log-likelihood of the spatio temporal model. 
##' 
##' @section Warning: \code{loglikeSTnaive} may take long to run. However for
##'   some problems with many locations and short time series
##'   \code{loglikeSTnaive} could be faster than \code{loglikeST}.
##' 
##' @example Rd_examples/Ex_loglikeST.R
##'
##' @author Johan Lindström
##' 
##' @family STmodel functions
##' @family likelihood functions
##' @family estimation functions
##' @export
loglikeST <- function(x=NULL, STmodel, type="p", x.fixed=NULL){
  ##check class belonging
  stCheckClass(STmodel, "STmodel", name="STmodel")
  ##first ensure that type is lower case
  type <- tolower(type)
  ##check if type is valid and if x.fixed should be expanded
  x <- stCheckLoglikeIn(x, x.fixed, type)
    
  ##if x is null or contains NA
  if( is.null(x) || any(is.na(x)) ){
    ##return the expected variable names
    return( loglikeSTnames(STmodel, all=(type=="f")) )
  }
  ##else, calculate loglikelihood
  
  ##first figure out a bunch of dimensions
  dimensions <- loglikeSTdim(STmodel)
  ##check parameter sizes
  if((type=="f" && length(x)!=dimensions$nparam) ||
     (type!="f" && length(x)!=dimensions$nparam.cov)){
    stop("Number of parameters, length(x), incorrect.")
  }
  
  ##extract parameters from x
  tmp <- loglikeSTgetPars(x, STmodel)
  if(type=="f"){
    gamma <- tmp$gamma
    alpha <- tmp$alpha
  }
  cov.pars.beta <- tmp$cov.beta
  cov.pars.nu <- tmp$cov.nu

  ##extract the observations
  Y <- STmodel$obs$obs
  ##extract sparse matrices
  F <- expandF(STmodel$F, STmodel$obs$idx, n.loc=dimensions$n.obs)
  if( type=="f" ){
    ##calculate the land use regression for the temporal trends, as column
    mu.B <- matrix(calc.mu.B(STmodel$LUR, alpha), ncol=1)
    if( dimensions$L!=0 ){ ##we have spatio temporal covariates, subtract these
      Y <- Y - STmodel$ST %*% gamma
    }
  }

  ##create covariance matrices, beta-field
  sigma.B <- makeSigmaB(cov.pars.beta$pars, dist = STmodel$D.beta,
                        type = STmodel$cov.beta$covf,
                        nugget = cov.pars.beta$nugget)
  ##and nu-field
  sigma.nu <- makeSigmaNu(cov.pars.nu$pars, dist = STmodel$D.nu,
                          type = STmodel$cov.nu$covf,
                          nugget = cov.pars.nu$nugget,
                          random.effect = cov.pars.nu$random.effect,
                          blocks1 = STmodel$nt, ind1 = STmodel$obs$idx,
                          sparse=TRUE)
  ##calculate block cholesky factor of the matrices
  ##storing in-place to conserve memory
  sigma.B <- try( makeCholBlock(sigma.B, n.blocks=dimensions$m),
                 silent=TRUE)
  if(class(sigma.B)=="try-error"){
    return(-.Machine$double.xmax)
  }
  sigma.nu <- try( chol(sigma.nu), silent=TRUE)
  if(class(sigma.nu)=="try-error"){
    return(-.Machine$double.xmax)
  }
  ##loglikelihood calculations follow:
  ##-log(det(sigma.B)^.5)
  l <- -sumLogDiag( sigma.B )
  ##-log(det(sigma.nu)^.5)
  l <- l -sumLog(diag( sigma.nu ))

  ##invert the matrices
  ##calculate inverse of sigma.B
  sigma.B <- invCholBlock(sigma.B, n.blocks=dimensions$m)
  ##calculate inverse of sigma.nu
  sigma.nu <- chol2inv(sigma.nu)
  ##inv(sigma.nu) * Y
  i.sR.Y <- as.matrix(sigma.nu %*% Y)
  
  if(type=="f"){
    ##calculate inv(matrix)*mu
    ##inv(sigma.B) * mu.B
    i.sB.mu.B <- blockMult(sigma.B, mu.B, n.blocks=dimensions$m)
    ##-1/2 mu.B' * inv(sigma.B) * mu.B
    l <- l - dotProd(i.sB.mu.B, mu.B)/2
    ##-1/2 Y' * inv(sigma.nu) * Y
    l <- l - dotProd(i.sR.Y, Y)/2
  }else{
    ##inv(sigma.B) * X
    iS.X <- calc.iS.X(STmodel$LUR, sigma.B)
    ##inv(sigma.nu) * M
    if(dimensions$L!=0){
      i.sR.M <- as.matrix(sigma.nu %*% STmodel$ST)
    }
    ##F'*inv(sigma.nu)*Y
    F.i.sR.Y <- calc.tFX(STmodel$F, i.sR.Y, STmodel$obs$idx,
                         n.loc=dimensions$n.obs)
    ##F'*inv(sigma.nu)*M
    if(dimensions$L!=0){
      F.i.sR.M <- calc.tFX(STmodel$F, i.sR.M, STmodel$obs$idx,
                           n.loc=dimensions$n.obs)
    }
  }##if(type=="f"){...}else{...}
  
  ##inv(sigma.B|Y) = F'*inv(sigma.nu)*F + inv(sigma.B)
  sigma.B.Y <- as.matrix( t(F) %*% sigma.nu %*% F) + sigma.B
  ##calculate cholesky factor of inv(sigma.B|Y)
  sigma.B.Y <- try( chol(sigma.B.Y), silent=TRUE)
  if( class(sigma.B.Y)=="try-error" ){
    return(-.Machine$double.xmax)
  }

  ##-log(det( inv(sigma.B|Y) )^.5)
  l <- l - sumLogDiag( sigma.B.Y )

  if(type=="f"){
    ##mu.B.Y = inv(sigma.B)*mu.B + F'*inv(sigma.nu)*Y
    mu.B.Y <- i.sB.mu.B + calc.tFX(STmodel$F, i.sR.Y, STmodel$obs$idx,
                                   n.loc=dimensions$n.obs)
    ##calculate inv(R)'*mu.B.Y
    mu.B.Y <- solveTriBlock(sigma.B.Y, mu.B.Y, transpose=TRUE)
    ##+1/2 mu.B.Y * inv(i.sigma.B.Y) * mu.B.Y
    l <- l + norm2(mu.B.Y)/2
  }else{
    ##chol(inv(sigma.B.Y))^-T * F' * inv(sigma.nu) * Y
    Y.hat <- solveTriBlock(sigma.B.Y, F.i.sR.Y, transpose=TRUE)
    ##chol(inv(sigma.B.Y))^-T * F' * inv(sigma.nu) * M
    if(dimensions$L!=0){
      M.hat <- solveTriBlock(sigma.B.Y, F.i.sR.M, transpose=TRUE)
    }
    ##calculate inverse of inv(sigma.B|Y) (keep in place)
    sigma.B.Y <- invCholBlock(sigma.B.Y)
    ##sigma.B.Y * inv(sigma.B) * X
    sigma.B.Y.iS.X <- sigma.B.Y %*% iS.X
    ##inv(sigma.alpha|Y) = X'*inv(sigma.B)*X -
    ##       X'*inv(sigma.B)*sigma.B|Y*inv(sigma.B)*X
    i.sigma.alpha.Y <- -t(iS.X) %*% sigma.B.Y.iS.X +
      calc.X.iS.X(STmodel$LUR, iS.X)
    ##calculate cholesky factor of inv(sigma.alpha|Y)
    i.sigma.alpha.Y <- try( chol(i.sigma.alpha.Y), silent=TRUE)
    if(class(i.sigma.alpha.Y)=="try-error"){
      return(-.Machine$double.xmax)
    }
    if(type=="r"){
      l <- l - sumLogDiag( i.sigma.alpha.Y )
    }
    ## X' * inv(sigma.B) * sigma.B.Y * F' * inv(sigma.nu) * Y
    Y.hat.2 <- t(sigma.B.Y.iS.X) %*% F.i.sR.Y
    ## chol(inv(sigma.alpha.Y))^-T * Y.hat.2
    Y.hat.2 <- solveTriBlock(i.sigma.alpha.Y, Y.hat.2, transpose=TRUE)
    ##Y'*sigma_hat*Y
    Y.sigma.hat.Y <- dotProd(Y,i.sR.Y) - norm2(Y.hat) - norm2(Y.hat.2)
    l <- l - Y.sigma.hat.Y/2
  
    ##parts relevant only with spatio-temporal model
    if(dimensions$L!=0){
      ## X' * inv(sigma.B) * sigma.B.Y * F' * inv(sigma.nu) * M
      M.hat.2 <- t(sigma.B.Y.iS.X) %*% F.i.sR.M
      ## chol(inv(sigma.alpha.Y))^-T * M.hat.2
      M.hat.2 <- solveTriBlock(i.sigma.alpha.Y, M.hat.2, transpose=TRUE)
    
      Y.sigma.hat.M <- (t(Y) %*% i.sR.M - t(Y.hat) %*% M.hat -
                        t(Y.hat.2) %*% M.hat.2)
      Y.sigma.hat.M <- t(Y.sigma.hat.M)
      M.sigma.hat.M <- (t(STmodel$ST) %*% i.sR.M - t(M.hat) %*% M.hat -
                        t(M.hat.2) %*% M.hat.2)
      ##calculate cholesky factor of M.sigma.hat.M
      M.sigma.hat.M <- try( chol(M.sigma.hat.M), silent=TRUE)
      if(class(M.sigma.hat.M)=="try-error"){
        return(-.Machine$double.xmax)
      }
      ##-log(det( M.sigma.hat.M )^.5)
      if(type=="r"){ 
        l <- l - sumLogDiag( M.sigma.hat.M )
      }

      ## chol(inv(sigma.alpha.Y))^-T * Y.hat.2
      Y.sigma.hat.M <- solveTriBlock(M.sigma.hat.M, Y.sigma.hat.M, transpose=TRUE)
      l <- l + norm2(Y.sigma.hat.M)/2
    }##if(dimensions$L!=0)
  }##if(type=="f") ... else ...

  ##ensure that l is treated as a double value (and not a sparse matrix)
  l <- as.double(l)
  ##add safe guard (infinite or nan results almost always due to
  ##matrix inprecision, lets return a very small value)
  if( !is.finite(l) )
    l <- -.Machine$double.xmax
  return(l)
}##function loglikeST


###############################################
## Log-likelihood using the full formulation ##
###############################################

##' @rdname loglikeST
##' @export
loglikeSTnaive <- function(x=NULL, STmodel, type="p", x.fixed=NULL){
  ##check class belonging
  stCheckClass(STmodel, "STmodel", name="STmodel")
  ##first ensure that type is lower case
  type <- tolower(type)
  ##check if type is valid and if x.fixed should be expanded
  x <- stCheckLoglikeIn(x, x.fixed, type)
  
  ##if x is null or contains NA
  if( is.null(x) || any(is.na(x)) ){
    ##return the expected variable names
    return( loglikeSTnames(STmodel, all=(type=="f")) )
  }
  ##else, calculate loglikelihood
  
  ##first figure out a bunch of dimensions
  dimensions <- loglikeSTdim(STmodel)
  ##check parameter sizes
  if((type=="f" && length(x)!=dimensions$nparam) ||
     (type!="f" && length(x)!=dimensions$nparam.cov)){
    stop("Number of parameters, length(x), incorrect.")
  }
    
  ##extract parameters from x
  tmp <- loglikeSTgetPars(x, STmodel)
  if(type=="f"){
    gamma <- tmp$gamma
    alpha <- tmp$alpha
  }
  cov.pars.beta <- tmp$cov.beta
  cov.pars.nu <- tmp$cov.nu

  ##extract the observations
  Y <- STmodel$obs$obs
  ##extract sparse matrices
  LUR <- Matrix::bdiag(STmodel$LUR)
  F <- expandF(STmodel$F, STmodel$obs$idx, n.loc=dimensions$n.obs)
  if( type=="f" ){
    ##subtract mean value from observations
    Y <- as.matrix( Y - (F %*% (LUR %*% unlist(alpha))) )
    if(dimensions$L!=0){
      ##also subtract spatio-termporal covariate
      Y <- Y - STmodel$ST %*% gamma
    }
  }##if( type=="f" )
  
  ##create covariance matrices, beta-field
  sigma.B.full <- makeSigmaB(cov.pars.beta$pars, dist = STmodel$D.beta,
                             type = STmodel$cov.beta$covf,
                             nugget = cov.pars.beta$nugget, sparse=TRUE)
  sigma.B.full <- F %*% Matrix(sigma.B.full %*% t(F))

  ##and nu-field (do NOT use sparse, due to the full matrix below)
  sigma.nu <- makeSigmaNu(cov.pars.nu$pars, dist = STmodel$D.nu,
                          type = STmodel$cov.nu$covf,
                          nugget = cov.pars.nu$nugget,
                          random.effect = cov.pars.nu$random.effect,
                          blocks1 = STmodel$nt, ind1 = STmodel$obs$idx,
                          sparse=FALSE)

  
  ##Total covariance matrix
  sigma.nu <- sigma.nu + as.matrix(sigma.B.full)
  ##calculate (block) cholesky factor of the matrices
  ##storing in-place to conserve memory
  sigma.nu <- try( chol(sigma.nu), silent=TRUE)
  if(class(sigma.nu)=="try-error"){
    return(-.Machine$double.xmax)
  }

  ##loglikelihood calculations follow:
  ##-log(det(sigma.nu)^.5)
  l <- -sumLogDiag( sigma.nu )
  ##calculate if(type=="f") inv(R)'*(Y-mean.val) else inv(R)'*Y
  Y <- solveTriBlock(sigma.nu, Y, transpose=TRUE)
  ##if(type=="f")  -1/2 (Y-mean)' * inv(sigma.nu) * (Y-mean)
  ##     else      -1/2 Y' * inv(sigma.nu) * Y
  l <- l - norm2(Y)/2
  if(type!="f"){
    ##Create the Ftmp = [FX M] matrix
    Ftmp <- as.matrix( F %*% LUR )
    ##Add the spatio-temporal covariate (if it exists)
    if( dimensions$L!=0 ){
      Ftmp <- cbind(Ftmp, STmodel$ST)
    }
  
    ##calculate inv(R)'*Ftmp
    Ftmp <- solveTriBlock(sigma.nu, Ftmp, transpose=TRUE)
    ##calculate Ftmp'*inv(Sigma)*Y
    FY <- t(Ftmp) %*% Y
    ##calculate [FX M]*invSigma*[FX M]'
    sigma.alt <- t(Ftmp) %*% Ftmp
    ##calculate cholesky factor (storing in-place to conserve memory)
    sigma.alt <- try( chol(sigma.alt), silent=TRUE)
    if( class(sigma.alt)=="try-error" ){
      return(-.Machine$double.xmax)
    }
    ##-log(det(sigma.alt)^.5)
    if(type=="r"){
      l <- l - sumLogDiag( sigma.alt )
    }
    ##calculate inv(R)'*Y
    FY <- solveTriBlock(sigma.alt, FY, transpose=TRUE)
    ##+1/2 FY' * inv(sigma.alt) * FY
    l <- l + norm2(FY)/2
  }##if(type!="f")
  
  ##ensure that l is treated as a double value (and not a sparse matrix)
  l <- as.double(l)
  
  ##add safe guard (infinite or nan results almost always due to
  ##matrix inprecision, lets return a very small value)
  if( !is.finite(l) )
    l <- -.Machine$double.xmax
  return(l)
}##function loglikeSTnaive
