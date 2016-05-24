#' Computes various types of residuals
#'
#' @description Computes several types of residuals for objects of class \code{clme}.
#'
#' @rdname clme_resids
#'
#' @param formula a formula expression. The constrained effect(s) must come before any unconstrained covariates on the right-hand side of the expression. The first \code{ncon} terms will be assumed to be constrained.
#' @param data data frame containing the variables in the model. 
#' @param gfix optional vector of group levels for residual variances. Data should be sorted by this value.
#' @param ncon the number of variables in \code{formula} that are constrained.
#'
#' @details 
#' For fixed-effects models \eqn{Y = X\beta + \epsilon}{Y = X*b + e}, residuals are given as \eqn{\hat{e} = Y - X\hat{\beta}}{ ehat = Y - X*betahat}.
#' For mixed-effects models \eqn{Y = X\beta + + U\xi + \epsilon}{Y = X*b + U*xi + e}, three types of residuals are available.
#' \eqn{PA = Y - X\hat{\beta}}{ PA = Y - X*betahat}\\
#' \eqn{SS = U\hat{\xi}}{ SS = U*xihat}\\
#' \eqn{FM = Y - X\hat{\beta} - U\hat{\xi}}{ FM = Y - X*betahat - U*xihat}
#' 
#' @return
#' List containing the elements \code{PA}, \code{SS}, \code{FM}, \code{cov.theta}, \code{xi}, \code{ssq}, \code{tsq}.
#' \code{PA}, \code{SS}, \code{FM} are defined above (for fixed-effects models, the residuals are only \code{PA}). Then \code{cov.theta} is the unconstrained covariance matrix of the fixed-effects coefficients, \code{xi} is the vector of random effect estimates, and \code{ssq} and \code{tsq} are unconstrained estimates of the variance components.
#' 
#' @note
#' There are few error catches in these functions. If only the EM estimates are desired, 
#' users are recommended to run \code{\link{clme}} setting \code{nsim=0}.
#' 
#' By default, homogeneous variances are assumed for the residuals and (if included) 
#' random effects. Heterogeneity can be induced using the arguments \code{Nks} and \code{Qs},
#'  which refer to the vectors \eqn{ (n_{1}, n_{2}, \ldots, n_{k}) }{(n1, n2 ,... , nk)} and
#'   \eqn{ (c_{1}, c_{2}, \ldots, c_{q}) }{(c1, c2 ,... , cq)}, respectively. See 
#'   \code{\link{CLME-package}} for further explanation the model and these values.
#' 
#' See \code{\link{w.stat}} and \code{\link{lrt.stat}} for more details on using custom 
#' test statistics.
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' \dontrun{
#' data( rat.blood )
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' 
#' clme.out <- clme_resids(mcv ~ time + temp + sex + (1|id), data = rat.blood, ncon=1 )
#' }
#' 
#' @importFrom MASS ginv
#' @export
#' 
#'
clme_resids <-
function( formula, data, gfix=NULL, ncon=1 ){
  
  ##
  ## I should consider a way to condense this so it doesn't need to 
  ## be copied between here, clme() and resid_boot().
  ##
  suppressMessages( mmat <- model_terms_clme( formula, data, ncon ) )
  formula2 <- mmat$formula
  Y  <- mmat$Y
  P1 <- mmat$P1
  X1 <- mmat$X1
  X2 <- mmat$X2
  U  <- mmat$U
  
  if( is.null(gfix) ){
    gfix <- rep("Residual", nrow(X1)) 
  } else{
    data <- with( data, data[order(gfix),])
  }
  Nks <- table(gfix)
  
  if( !is.null(U) ){
    if( !is.null(mmat$REidx) ){
      Qs        <- table( mmat$REidx )
      names(Qs) <- mmat$REnames
    } else{
      Qs  <- table( rep("tsq", ncol(U)) )
    }    
  } else{
    Qs <- NULL
  }
  
  #############################################
  
  N  <- sum(Nks)
  N1 <- 1 + cumsum(Nks) - Nks
  N2 <- cumsum(Nks)
  
  Q  <- length(Qs)
  Q1 <- 1 + cumsum(Qs) - Qs
  Q2 <- cumsum(Qs)
  
  X  <- as.matrix( cbind(X1, X2) )
  K  <- length(Nks)
  
  
  # Initial values
  theta <- ginv( t(X)%*%X )%*%( t(X)%*%Y )
  
    ssq   <- vector()
    for( k in 1:K ){
      Yk <- Y[ N1[k]:N2[k] ]
      Xk <- X[ N1[k]:N2[k],]  
      ssq[k] <- sum( (Yk - Xk%*%theta)^2 ) / (Nks[k])
    }
  
  ## Obtain the estimates of epsilon and delta
  ssqvec  <- rep(ssq,Nks)    
  XSiX <- t(X) %*% (X/ssqvec)
  XSiY <- t(X) %*% (Y/ssqvec)
  
  if( Q > 0 ){
    mq.phi <- minque( Y=Y, X1=X1, X2=X2, U=U, Nks=Nks, Qs=Qs )[1:Q]
    tsq    <- mq.phi
    
    tsqvec  <- rep(tsq,Qs)    
    C       <- U * tsqvec
    
    U1         <- apply( U , 2 , FUN=function(x,sq){x*sq} , 1/sqrt(ssqvec) )
    tusu       <- t(U1) %*% U1
    diag(tusu) <- diag(tusu) + 1/tsqvec
    tusui      <- solve(tusu)   
    XSiU <- t(X) %*% (U/ssqvec)
    USiY <- t(U) %*% (Y/ssqvec)
    XPiX <- XSiX - XSiU%*%tusui%*%t(XSiU)
    XPiY <- XSiY - XSiU%*%(tusui%*%USiY)
  } else{
    XPiX <- XSiX
    XPiY <- XSiY
  }
  
  # H <- X%*%ginv( XPiX )%*%(t(X)%*%PsiI)
  #eps  <- c( Y - H%*%Y )  
  Yhat <- X %*% ginv( XPiX )%*%XPiY
  PA   <- c(Y - Yhat)
  
  if( Q > 0 ){
    # xi    <- c( t(C) %*% PsiI %*% eps )  
    USiR <- t(U) %*% (PA/ssqvec)
    USiU <- t(U) %*% (U/ssqvec)
    xi   <- tsqvec * ( USiR - USiU%*%(tusui%*%USiR) )
    SS <- U %*% xi
    FM <- PA - SS

    resid.out <- list( PA=c(PA), SS=c(SS), FM=c(FM), cov.theta=solve(XPiX),
                       xi=c(xi), ssq=ssq, tsq=tsq )
  } else{
    resid.out <- list( PA=c(PA), cov.theta=XPiX, ssq=ssq )
  }
  
  return( resid.out )
    
}












