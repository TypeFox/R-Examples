#' @export
rotate_ranef <- function(.mod, ...){
  UseMethod("rotate_ranef", .mod)
}

#' @export
#' @rdname rotate_ranef.mer
#' @method rotate_ranef default
rotate_ranef.default <- function(.mod, ...){
  stop(paste("there is no rotate_ranef() method for objects of class",
             paste(class(.mod), collapse=", ")))
}

#' Calculate s-dimensional rotated random effects
#' 
#' This function calculates reduced dimensional rotated random effects. 
#' The rotation reduces the influence of the residuals from other levels
#' of the model so that distributional assessment of the resulting
#' random effects is possible.
#' 
#' @export
#' @method rotate_ranef mer
#' @S3method rotate_ranef mer
#' @aliases rotate_ranef
#' @param .mod an object of class \code{mer} or \code{lmerMod}.
#' @param .L a matrix defining which combination of random effects are of interest.
#' @param s the dimension of the subspace of interest.
#' @param .varimax if \code{.varimax = TRUE} than the raw varimax rotation 
#'   will be applied to the resulting rotation.
#' @param ... do not use
#' @author Adam Loy \email{loyad01@@gmail.com}
#' @references
#' Loy, A. & Hofmann, H. (in press). Are you Normal? 
#' The Problem of Confounded Residual Structures in Hierarchical Linear Models.
#' \emph{Journal of Computational and Graphical Statistics}.
rotate_ranef.mer <- function(.mod, .L, s = NULL, .varimax = FALSE, ...) {
  y <- .mod@y
  X <- lme4::getME(.mod, "X")
  Z <- BlockZ(.mod)
  
  n <- nrow(X)
  p <- ncol(X)
  ngrps <- unname( summary(.mod)@ngrps )
  
  vc <- lme4::VarCorr(.mod)
  Di <- bdiag( lme4::VarCorr(.mod) ) / (unname(attr(vc, "sc")))^2
  D  <- kronecker( Diagonal(ngrps), Di )
  
  Aslot <- .mod@A # ZDZ'
  zdzt <- crossprod( .mod@A )
  V  <- Diagonal( n ) + zdzt
  V.chol <- chol( V )
  Vinv  <- chol2inv( V.chol ) 
  
  XVXinv <- solve( t(X) %*% Vinv %*% X )
  VinvX  <- Vinv %*% X
  M      <- VinvX %*% XVXinv %*% t(VinvX)
  P      <- .Call("cxxmatsub", BB = as.matrix(Vinv), CC = as.matrix(M), 
                  PACKAGE = "HLMdiag")
  
  betahat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
  mr <- y - X %*% betahat
  
  bvec <- D %*% t(Z) %*% Vinv %*% mr
    
  pzdl <- P %*% Z %*% D %*% .L
  A <- crossprod( pzdl )
  B <- t(.L) %*% D %*% t(Z) %*% P %*% Z %*% D %*% .L ## diagnostic se
  W <- try( mcrotate(A, B, s) )
  if( class(W) == "try-error") {W <- NA} else {W <- as.matrix(W)}
    
  if( .varimax == TRUE) {
    W <- try( varimax(W, normalize = FALSE)$loadings )
    if( class(W) == "try-error" ) W <- NA 
  }
    
  return( as.numeric( t(W) %*% as.numeric( t(.L) %*% bvec ) ) )
}


#' @export
#' @rdname rotate_ranef.mer
#' @method rotate_ranef lmerMod
#' @S3method rotate_ranef lmerMod
rotate_ranef.lmerMod <- function(.mod, .L, s = NULL, .varimax = FALSE, ...) {
  y <- .mod@resp$y
  X <- lme4::getME(.mod, "X")
  Z <- lme4::getME(.mod, "Z")
  
  n <- nrow(X)
  p <- ncol(X)
  ngrps <- unname( summary(.mod)$ngrps )
  
  vc <- lme4::VarCorr(.mod)
  Di <- bdiag( lme4::VarCorr(.mod) ) / (unname(attr(vc, "sc")))^2
  D  <- kronecker( Diagonal(ngrps), Di )
  
  zdzt <- crossprod( lme4::getME(.mod, "A") )
  V  <- Diagonal( n ) + zdzt
  V.chol <- chol( V )
  Vinv  <- chol2inv( V.chol ) 
  
  XVXinv <- solve( t(X) %*% Vinv %*% X )
  VinvX  <- Vinv %*% X
  M      <- VinvX %*% XVXinv %*% t(VinvX)
  P      <- .Call("cxxmatsub", BB = as.matrix(Vinv), CC = as.matrix(M), 
                  PACKAGE = "HLMdiag")
  
  betahat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
  mr <- y - X %*% betahat
  
  bvec <- D %*% t(Z) %*% Vinv %*% mr
  
  pzdl <- P %*% Z %*% D %*% .L
  A <- crossprod( pzdl )
  B <- t(.L) %*% D %*% t(Z) %*% P %*% Z %*% D %*% .L ## diagnostic se
  W <- try( mcrotate(A, B, s) )
  if( class(W) == "try-error") {W <- NA} else {W <- as.matrix(W)}
  
  if( .varimax == TRUE) {
    W <- try( varimax(W, normalize = FALSE)$loadings )
    if( class(W) == "try-error" ) W <- NA 
  }
  
  return( as.numeric( t(W) %*% as.numeric( t(.L) %*% bvec ) ) )
}


mcrotate <- function(A, B, s) {
  r <- rankMatrix(B)
  if(is.null(s)) s <- r
  
  B.svd <- svd(B)
  Cr.diag <- B.svd$d[1:r]
  Tr <- B.svd$u[, 1:r]
  
  A.star <- Diagonal( x = 1 / sqrt(Cr.diag) ) %*% t(Tr) %*% 
    A %*% 
    Tr %*% Diagonal( x = 1 / sqrt(Cr.diag) )
  
  A.star.svd <- svd( A.star )
  
  index <- seq(r, length.out = s, by = -1)
  index <- sort(index[index >= 0])
  W <- Tr %*% Diagonal( x = 1 / sqrt( Cr.diag ) ) %*% A.star.svd$u[,index]
  
  return(W)
}


#' @export
#' @rdname rotate_ranef.mer
#' @method rotate_ranef lme
#' @S3method rotate_ranef lme
rotate_ranef.lme <- function(.mod, .L, s = NULL, .varimax = FALSE, ...) {
  design.info <- .extract.lmeDesign(.mod)
  
  y <- design.info$y
  X <- design.info$X
  Z <- Matrix( design.info$Z )
  D <- Matrix( design.info$Vr )
  
  V  <- .extractV.lme( .mod )
  V.chol <- chol( V )
  Vinv  <- chol2inv( V.chol ) 
  
  XVXinv <- solve( t(X) %*% Vinv %*% X )
  VinvX  <- Vinv %*% X
  M      <- VinvX %*% XVXinv %*% t(VinvX)
  P      <- .Call("cxxmatsub", BB = as.matrix(Vinv), CC = as.matrix(M), 
                  PACKAGE = "HLMdiag")
  
  betahat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
  mr <- y - X %*% betahat
  
  bvec <- D %*% t(Z) %*% Vinv %*% mr
  
  pzdl <- P %*% Z %*% D %*% .L
  A <- crossprod( pzdl )
  B <- t(.L) %*% D %*% t(Z) %*% P %*% Z %*% D %*% .L ## diagnostic se
  W <- try( mcrotate(A, B, s) )
  if( class(W) == "try-error") {W <- NA} else {W <- as.matrix(W)}
  
  if( .varimax == TRUE) {
    W <- try( varimax(W, normalize = FALSE)$loadings )
    if( class(W) == "try-error" ) W <- NA 
  }
  
  return( as.numeric( t(W) %*% as.numeric( t(.L) %*% bvec ) ) )
}

# I only need this if the browser() stays in the RLRsim code
.extract.lmeDesign <- function (m) 
{
  start.level = 1
  data <- if (any(!complete.cases(m$data))) {
    warning("Removing incomplete cases from supplied data.")
    m$data[complete.cases(m$data), ]
  }
  else m$data
  grps <- nlme::getGroups(m)
  n <- length(grps)
  X <- list()
  grp.dims <- m$dims$ncol
  Zt <- model.matrix(m$modelStruct$reStruct, data)
  cov <- as.matrix(m$modelStruct$reStruct)
  i.col <- 1
  n.levels <- length(m$groups)
  Z <- matrix(0, n, 0)
  if (start.level <= n.levels) {
    for (i in 1:(n.levels - start.level + 1)) {
      if (length(levels(m$groups[[n.levels - i + 1]])) != 
            1) {
        X[[1]] <- model.matrix(~m$groups[[n.levels - 
                                            i + 1]] - 1, contrasts.arg = c("contr.treatment", 
                                                                           "contr.treatment"))
      }
      else X[[1]] <- matrix(1, n, 1)
      X[[2]] <- as.matrix(Zt[, i.col:(i.col + grp.dims[i] - 
                                        1)])
      i.col <- i.col + grp.dims[i]
      Z <- cbind(mgcv::tensor.prod.model.matrix(X), Z)
    }
    Vr <- matrix(0, ncol(Z), ncol(Z))
    start <- 1
    for (i in 1:(n.levels - start.level + 1)) {
      k <- n.levels - i + 1
      for (j in 1:m$dims$ngrps[i]) {
        stop <- start + ncol(cov[[k]]) - 1
        Vr[ncol(Z) + 1 - (stop:start), ncol(Z) + 1 - 
             (stop:start)] <- cov[[k]]
        start <- stop + 1
      }
    }
  }
  X <- if (class(m$call$fixed) == "name" && !is.null(m$data$X)) {
    m$data$X
  }
  else {
    model.matrix(formula(eval(m$call$fixed)), data)
  }
  y <- as.vector(matrix(m$residuals, ncol = NCOL(m$residuals))[, 
                                                               NCOL(m$residuals)] + matrix(m$fitted, ncol = NCOL(m$fitted))[, 
                                                                                                                            NCOL(m$fitted)])
  return(list(Vr = Vr, X = X, Z = Z, sigmasq = m$sigma^2, lambda = unique(diag(Vr)), 
              y = y, k = n.levels))
}