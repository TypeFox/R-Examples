HL_process.args <- function (...) { ## from gsl package
  a <- list(...)
  attr <- attributes(a[[which.max(unlist(lapply(a, length)))]])
  a <- lapply(a, as.vector)
  out <- do.call(rbind, a)
  out <- list(`1`=out[1,],`2`=out[2,]) ## == out <- split(out, row(out)) but profiling shows the latter is slow
  names(out) <- paste("arg", 1:length(a), sep = "")
  return(c(out, attr = list(attr)))
}

HL_strictify <- function (val, status) { ## from gsl package
    val[status < 0] <- NaN ##FR inverted the test // gsl package !!
    return(val)
}

bessel_lnKnu <- function (nu, x, give = FALSE, strict = TRUE) { ## from bessel_lnKnu in gsl package
    jj <- HL_process.args(nu, x)
    nu.vec <- jj$arg1
    x.vec <- jj$arg2
    attr <- jj$attr
    jj <- .C("bessel_lnKnu_e", as.double(nu.vec), as.double(x.vec), 
        as.integer(length(x.vec)), val = as.double(x.vec), err = as.double(x.vec), 
        status = as.integer(0 * x.vec), PACKAGE = "spaMM")
    val <- jj$val
    err <- jj$err
    status <- jj$status
    attributes(val) <- attr
    attributes(err) <- attr
    attributes(status) <- attr
    if (strict) {
        val <- HL_strictify(val, status)
    }
    if (give) {
        return(list(val = val, err = err, status = status))
    }
    else {
        return(val)
    }
}

"MaternCorr" <- function(d, rho=1, smoothness, nu=smoothness, Nugget=0L) UseMethod("MaternCorr") 

Matern.corr <- MaternCorr ## for back compat as it is in spMMjob.R

MaternCorr.ff <- function (d, rho=1, smoothness, nu=smoothness, Nugget=0L) { ## rho is alpha in fields
  ## ideally (but not necess) on a 'dist' so the diagonal is not  manipulated 
    if (any(d[] < 0)) 
        stop("distance argument must be nonnegative")
    dscal <- ff::ff(vmode="double",dim=dim(d))
    dscal[] <- d[] * rho
    isd0 <- d[] == 0L ## regular matrix
    dscal[][isd0] <- 1e-10 ## avoids errors on distance =0; but the resulting corrvals can be visibly < 1 for small nu ## FR->FR make value dependent on rho, nu ?
    logcon <- (nu - 1)*log(2)+ lgamma(nu) 
    corrvals <- ff::ff(vmode="double",dim=dim(d))
    ## bessel_lnKnu -> various local copies -> memory issues
    for (it in seq(nrow(corrvals))) corrvals[it,] <- - logcon + nu*log(dscal[it,])+ bessel_lnKnu(x=dscal[it,], nu=nu) ## 
    corrvals[] <- exp(corrvals[]) 
    corrvals[][!isd0] <- (1-Nugget)* corrvals[][!isd0]
    corrvals[][isd0] <- 1 ## 
    corrvals[][corrvals[] < 1e-16] <- 0L ## an attempt to deal with problem in chol/ldl/svd which don't like 'nearly-identity' matrices
    attr(corrvals,"corr.model") <- "Matern"
    return(corrvals)
}

MaternCorr.default <- function (d, rho=1, smoothness, nu=smoothness, Nugget=0L) { ## rho is alpha in fields
  ## ideally (but not necess) on a 'dist' so the diagonal is not  manipulated 
  if (any(d < 0)) 
    stop("distance argument must be nonnegative")
  dscal <- d * rho
  isd0 <- d == 0L
  dscal[isd0] <- 1e-10 ## avoids errors on distance =0; but the resulting corrvals can be visibly < 1 for small nu ## FR->FR make value dependent on rho, nu ?
  logcon <- (nu - 1)*log(2)+ lgamma(nu) 
  corrvals <- - logcon + nu*log(dscal)+ bessel_lnKnu(x=dscal, nu=nu) ## 
  ##    corrvals <- - logcon + nu*log(dscal)+ log(besselK(x=dscal, nu=nu)) ## function from package gsl
  corrvals <- exp(corrvals) 
  corrvals[!isd0] <- (1-Nugget)* corrvals[!isd0]
  corrvals[isd0] <- 1 ## 
  corrvals[corrvals < 1e-16] <- 0L ## an attempt to deal with problem in chol/ldl/svd which don't like 'nearly-identity' matrices
  attr(corrvals,"corr.model") <- "Matern"
  return(corrvals)
}


## ess <- function(nu,d) {exp(-(d/(2*sqrt(nu)))^2)} ...


#### demo
if (F) {
 mym <- matrix(c(3,2,1,2,3,2,1,2,3),ncol=3)
 cL <- t(chol(mym))
 cL %*% t(cL)
 LDL <- eigen(mym,symmetric=T)
 eL <- LDL$vectors %*% diag(sqrt(LDL$values))    
 eL %*% t(eL)
 seL <- LDL$vectors %*% diag(sqrt(LDL$values)) %*% t(LDL$vectors)   ## symmetric
 seL %*% t(seL)
}

## FR->FR we also use this function once on d2hdv2 in HLfit...
`designL.from.Corr` <- function(m=NULL,symSVD=NULL,try.chol=TRUE,try.eigen=FALSE,threshold=1e-06,debug=FALSE,SVDfix=1/10) {
  ## cf return value: the code must compute 'L', and if the type of L is not chol, also 'corr d' and 'u'
  type <- NULL
  if ( ! is.null(symSVD)) {
    u <- symSVD$u ## local copy needed for processing attributes at the end of the function
    d <- symSVD$d ## of corr matrix !
    L <- ZWZt(u,sqrt(d))  
    type <- "symsvd"      
  } else if (try.chol) {
    if (.spaMM.data$options$USEEIGEN) {
      if (inherits(m,"blockDiag")) { 
        stop("designL.from.Corr code should be allowed again to handle blockDiag objects")
        #trychol <- RcppChol.blockDiag(m) ## cf pb RcppEigen / generic
      } else trychol <- RcppChol(m)
      if (trychol$Status==TRUE) { ## if chol computation successful
        L <- trychol$L
        type <- "chol"  
      } else {
        mreg <- m *(1-1e-8)
        diag(mreg) <- diag(mreg) + 1e-8 
        trychol <- RcppChol(mreg)
        if (trychol$Status==TRUE) {
          L <- trychol$L
          type <- "chol"  
        } ## else type remains NULL      
      } 
    } else { ## chol de R
      cholR <- try(chol(m),silent=TRUE) ## dim -> unique geo coordinates
      if (inherits(cholR,"try-error")) { ## attempt at regularization 050213
        mreg <- m *(1-1e-8)
        diag(mreg) <- diag(mreg) + 1e-8 
        cholR <- try(chol(mreg),silent=TRUE) 
        if ( ! inherits(cholR,"try-error")) type <- "chol"
      } else type <- "chol"
      if (!is.null(type)) L<-t(cholR) ## such that L %*% t(L) = the correlation matrix in both cases        
    }    
  }
  if ( is.null(type) ) { ## hence if none of the chol algos (nor symSVD input) has been used 
    ## slower by more stable. Not triangular but should not be a pb
    if (try.eigen) {
      LDL <- try(eigen(m,symmetric=TRUE),silent=TRUE) ## may _hang_ in R2.15.2 on nearly-I matrices
      if ( ! inherits(LDL,"try-error")) d <- LDL$values 
    }
    if ( (! try.eigen) || inherits(LDL,"try-error") || any(d < -1e-08)) {
      if (.spaMM.data$options$USEEIGEN) { ## see package irlba for SVD of sparse matrices
        symSVD <- selfAdjointSolverCpp(m) ## such that v= t(u)without any sign issue  
        u <- symSVD$u
        d <- symSVD$d  
        type <- "symsvd"          
      } else {
        SVD <- try(svd(m)) ## "the SVD implementation of Eigen (...) is not a particularly fast SVD method." (RcppEigen vignette)
        if (inherits(SVD,"try-error")) {
          print("spaMM retries SVD following 'svd' failure.") 
          ## numerically different but otherwise equivalent computation
          m <- diag(rep(1-SVDfix,ncol(m))) + m*SVDfix 
          SVD <- try(svd(m)) 
          if (! inherits(SVD,"try-error") ) SVD$d <- 1+ (SVD$d-1)/SVDfix 
        } 
        if (inherits(SVD,"try-error")) {
          ## typically m is I + a few large elements
          ## Most informative post: http://r.789695.n4.nabble.com/Observations-on-SVD-linpack-errors-and-a-workaround-td837282.html
          ####################### aide sur le .Rdata :  zut <- as.list(attr(resu,"condition")$call) est un HLCor call
          print("Singular value decomposition failed.") 
          print(" See documentation of the 'SVDfix' argument of 'designL.from.Corr'")
          print("   for ways to handle this.")
          return(try(stop(),silent=TRUE)) ## passes control to calling function
        } else {
          d <- SVD$d
          ## must be valid for sym (semi) PD matrices using U, V being eigenvectors of m %*% t(m)
          ## symm matrix => $u, $v match left and right eigenvectors of original Corr matrix
          ## FR->FR but they can be of opposite sign with negative $d...
          u <- SVD$u
          type <- "svd"          
        }         
      } ## svd by R
    } else { ## we have an LDL decomp
      u <- LDL$vectors
      d <- LDL$values
      type <- "eigen"
    }
    if (any(d< -1e-08)) {
      ## could occur for two reasons: wrong input matrix; or problem with R's svd which $d are the eigenvalues up to the sign, 
      ##   and thus $d can be <0 for eigenvalues >0 (very rare, observed in a 2x2 matrix) 
      mess <- pastefrom("correlation matrix has suspiciously large negative eigenvalue(s).")
      print(mess)
      return(try(stop(),silent=T)) ## passes control to calling function
    } else { ## we have a not-too-suspect decomp
      # d[d< threshold]<- threshold ## wrong for corrmats, would be OK for d2hdv2 computation which uses this function 
      if (any(d<threshold)) d <- threshold + (1-threshold) * d ## 17/05/2014
      L <- ZWZt(u,sqrt(d))
    }
  } 
  if ( ! is.null(m)) colnames(L) <- rownames(L) <- colnames(m) ## for checks in HLfit ## currently typically missing from symSVD   
  attr(L,"type") <- type
  if ( ! is.null(symSVD)) {
    attr(L,"corr.model") <- symSVD$corr.model
  } else attr(L,"corr.model") <- attr(m,"corr.model")
  ## add the 'sanitized' matrix decomp as attribute of the L matrix
  if (type != "chol") {
    decomp <- list(u=u,d=d)
    if ( ! is.null(symSVD)) decomp$adjd <- symSVD$adjd ## useful for SEM CAR; otherwise may be NULL
    attr(L,type) <- decomp
  }
  return(L)
} 

## we want to compute (1) [2e ligne de ChanK97 eq 11] the conditional covariance Cov - Cov inv(Cov+I) Cov  (we return a Cholesky factor of it)
## and (2) the conditional partial regression coeffs for Lv: Cov inv(Cov+I)
##  It turns out that the two are identical : 
## If Corr= LDL', cov= lam LDL', we want  L [lam D - lam^2 D^2/(lam D+I)] L'
## =  L [lam D/(lam D +I)] L' 
CondNormfn <- function(decomp,lambda) {
  diago <- decomp$d/(decomp$d+1/lambda)
  sqrtCondCovLv <- sweep(decomp$u,2,sqrt(diago),`*`); ## so that cond Corr = this.t(this)
  condLvReg <- tcrossprodCpp(sqrtCondCovLv) ## conditional regr = cond Corr
  # Un probleme est que la repres sqrtCondCovLv est loin d'être "numerically unique". 
  # Donc même si on a des distrib equivalentes pour differents sqrtCondCovLv 
  # (en particulier des condLvReg equivalents)
  # on va avoir sqrtCondCovLv %*% rnorm nettement divergents sous linux vs Windows 
  return(list(sqrtCondCovLv=sqrtCondCovLv,condLvReg=condLvReg))
} 

`make_scaled_dist` <- function(uniqueGeo,uniqueGeo2=NULL,distMatrix,rho,rho.mapping=seq_len(length(rho)),
                               dist.method="Euclidean",return_matrix=FALSE) {
  if (length(rho)>1L && dist.method!="Euclidean") { 
    mess <- pastefrom("'rho' length>1 not allowed for non-Euclidian distance.",prefix="(!) From ")
    stop(mess)
  }
  if ( missing(distMatrix) ) { ## ## then provide it (most of this fn's code)
    if ( missing(uniqueGeo) ) {
      mess <- pastefrom("missing(distMatrix) && missing(uniqueGeo).",prefix="(!) From ")
      stop(mess)
    } 
    distnrow <- nrow(uniqueGeo)
    distncol <- NROW(uniqueGeo2)
    scaled.dist <- NULL
    if ( distnrow*distncol > spaMM.getOption("ff_threshold") ) {
      if (requireNamespace("ff",quietly=TRUE)) {
        scaled.dist <- ff::ff( vmode ="double", dim=c(distnrow,distncol))
      } else message("Package 'ff' for large matrices may be needed but is not installed.")
    }
    if (dist.method=="Euclidean") {
      if (length(rho)==1L) {
        uniqueScal <- uniqueGeo * rho 
      } else if (ncol(uniqueGeo)==length(rho.mapping)) {
        uniqueScal <- t(t(uniqueGeo) * rho[rho.mapping]) ## valid for vectorial rho...
      } else {
        mess <- pastefrom("invalid length(rho[rho.mapping]).",prefix="(!) From ")
        print(mess)
        mess  <- paste("Length should be either 1 or",ncol(uniqueGeo))
        stop(mess)
      }
      if (! is.null(uniqueGeo2)) {
        if (length(rho)==1L) {
          uniqueScal2 <- uniqueGeo2 * rho  
        } else uniqueScal2 <- t(t(uniqueGeo2) * rho[rho.mapping]) 
      } else {
        uniqueScal2 <- NULL
      }
      if (inherits(scaled.dist,"ff_matrix")) {
        scaled.dist[] <- proxy::dist(x=uniqueScal,y=uniqueScal2, method=dist.method) 
      } else scaled.dist <- proxy::dist(x=uniqueScal,y=uniqueScal2, method=dist.method) 
    } else { ## not Euclidean
      if (inherits(scaled.dist,"ff_matrix")) {
        scaled.dist[] <- rho * proxy::dist(uniqueGeo,y=uniqueGeo2,method=dist.method) 
      } else scaled.dist <- rho * proxy::dist(uniqueGeo,y=uniqueGeo2,method=dist.method)  
    }
  } else { ## distMatrix provided
    scaled.dist <- rho * distMatrix
  }
  ## Here proxy::dist always returns a dist object, maybe dist(0) 
  ## but rho * dist(0) is numeric(0); we standardise it:
  if ( identical(scaled.dist, numeric(0))) scaled.dist <- dist(0)
  if (return_matrix) {
    if (inherits(scaled.dist,"dist")) {
      scaled.dist <- as.matrix(scaled.dist)
    } else if (inherits(scaled.dist,"crossdist")) scaled.dist <- scaled.dist[] ## []: same effect as what oen would expect from non-existent as.matrix.crossdist()
  }
  return(scaled.dist)
}

getDistMat <- function(object,scaled=FALSE) {
  if (! is.null(msd.arglist <- attr(object,"msd.arglist"))) {
    if ( ! scaled)  {
      msd.arglist$rho <- 1 
      msd.arglist$`rho.mapping` <- NULL 
    }
    return(do.call(make_scaled_dist,msd.arglist))
  } else {
    message("no Matern-correlated random effects")
    return(NULL)
  }
}

