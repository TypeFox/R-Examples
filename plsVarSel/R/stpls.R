#' @title Soft-Threshold PLS (ST-PLS)
#'
#' @description A soft-thresholding step in PLS algorithm (ST-PLS) based on
#' ideas from the nearest shrunken centroid method.
#'
#' @param ... arguments passed on to \code{mvrV}).
#' @param method choice between the default \code{stpls} and alternative \code{model.frame}.
#'
#' @details The ST-PLS approach is more or less identical to the Sparse-PLS presented
#' independently by Lè Cao et al. This implementation is an expansion of code from the
#' pls package.
#'  
#' @return Returns an object of class mvrV, simliar to to mvr object of the pls package.
#'
#' @author Solve Sæbø, Tahir Mehmood, Kristian Hovde Liland.
#'
#' @references S. Sæbø, T. Almøy, J. Aarøe, A.H. Aastveit, ST-PLS: a multi-dimensional 
#' nearest shrunken centroid type classifier via pls, Journal of Chemometrics 20 (2007) 54-62.
#'
#' @seealso \code{\link{VIP}} (SR/sMC/LW/RC), \code{\link{filterPLSR}}, \code{\link{spa_pls}}, 
#' \code{\link{stpls}}, \code{\link{truncation}}, \code{\link{bve_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{ipw_pls}}, \code{\link{ga_pls}}, \code{\link{rep_pls}}.
#'
#' @examples
#' data(yarn, package = "pls")
#' st <- stpls(density~NIR, ncomp=5, shrink=c(0.1,0.2), validation="CV", data=yarn)
#' summary(st)
#' 
#' @export
stpls <- function (..., method = c("stpls", "model.frame")){
  cl <- match.call()
  cl$method <- match.arg(method)
  cl[[1]] <- quote(mvrV)
  res <- eval(cl, parent.frame())
  ## Fix call component
  if (cl$method != "model.frame") res$call[[1]] <- as.name("stpls")
  if (missing(method)) res$call$method <- NULL
  res
}

stpls.fit <- function(X, Y, ncomp, shrink, stripped = FALSE,
                     tol = .Machine$double.eps^0.5, ...){
  ## Initialise
  Y <- as.matrix(Y)    
  if (!stripped) {
    ## Save dimnames
    dnX <- dimnames(X)
    dnY <- dimnames(Y)
  }
  ## Remove dimnames for performance (doesn't seem to matter; in fact,
  ## as far as it has any effect, it hurts a tiny bit in most situations.
  ## dimnames(X) <- dimnames(Y) <- NULL
  nobj  <- dim(X)[1]
  npred <- dim(X)[2]
  nresp <- dim(Y)[2]
  nshrink <- length(shrink)
  
  W <- P <- R <- array(0, dim = c(npred, ncomp, nshrink))
  #    W <- P <- matrix(0, nrow = npred, ncol = ncomp)
  tQ <- array(0, dim = c(ncomp, nresp, nshrink)) # Y loadings; "transposed"
  B  <- array(0, dim = c(npred, nresp, ncomp, nshrink))
  if (!stripped) {
    TT   <- U <- array(0, dim = c(nobj, ncomp, nshrink))
    tsqs <- matrix(0,nrow = ncomp, ncol = nshrink)          # t't
    fitted <- residuals <- array(0, dim = c(nobj, nresp, ncomp, nshrink))
  }
  ## C1
  Xmeans <- colMeans(X)
  globX  <- X - rep(Xmeans, each = nobj)
  Ymeans <- colMeans(Y)
  globY  <- Y - rep(Ymeans, each = nobj)
  ## Must be done here due to the deflation of X
  if (!stripped)
    Xtotvar <- sum(globX * globX)
  
  for(aa in 1:nshrink){
    X <- globX
    Y <- globY
    for(a in 1:ncomp) {
      ## Initial values:
      if (nresp == 1) {               # pls1
        u.a <- Y                    # FIXME: scale?
      } else {                        # pls2
        ## The coloumn of Y with largest sum of squares:
        u.a <- Y[,which.max(colSums(Y * Y))]
        t.a.old <- 0
      }
      repeat {
        ## C2.1
        w.a <- crossprod(X, u.a)
        w.a <- w.a/max(abs(w.a))                        #rescaling before soft-thresholding
        w.a <- as.vector(sign(w.a)*ifelse((abs(w.a)<shrink[aa]),0,abs(w.a)-shrink[aa]))   #Thresholding                        
        w.a <- w.a / sqrt(c(crossprod(w.a)))
        
        ## C2.2
        t.a <- X %*% w.a
        
        ## C2.3
        tsq  <- c(crossprod(t.a))
        t.tt <- t.a / tsq
        
        ## C2.4
        q.a <- crossprod(Y, t.tt)
        
        if (nresp == 1)
          break                   # pls1: no iteration
        
        ## C2.4b-c
        ## Convergence check for pls2:
        if (sum(abs((t.a - t.a.old) / t.a), na.rm = TRUE) < tol)
          break
        else {
          u.a <- Y %*% q.a / c(crossprod(q.a))
          t.a.old <- t.a          # Save for comparison
        }
      }
      
      ## C2.3 contd.
      p.a <- crossprod(X, t.tt)
      
      ## C2.5
      X <- X - t.a %*% t(p.a)
      Y <- Y - t.a %*% t(q.a)
      
      ## Save scores etc:
      W[,a,aa]  <- w.a
      P[,a,aa]  <- p.a
      tQ[a,,aa] <- q.a
      if (!stripped) {
        TT[,a,aa]  <- t.a
        U[,a,aa]   <- u.a
        tsqs[a,aa] <- tsq
        ## (For very tall, slim X and Y, X0 %*% B[,,a] is slightly faster,
        ## due to less overhead.)
        mytQ <- matrix(tQ[,,aa], nrow = ncomp, ncol = nresp)
        fitted[,,a,aa] <- TT[,1:a,aa] %*% mytQ[1:a,,drop=FALSE]
        residuals[,,a,aa] <- Y
      }
    }
    
  }
  
  for(aa in 1:nshrink){
    ## Calculate rotation matrix:
    PW <- crossprod(P[,,aa], W[,,aa])
    ## It is known that P^tW is right bi-diagonal (one response) or upper
    ## triangular (multiple responses), with all diagonal elements equal to 1.
    if (nresp == 1) {
      ## For single-response models, direct calculation of (P^tW)^-1 is
      ## simple, and faster than using backsolve.
      PWinv <- backsolve(PW, diag(ncomp))
      R[,,aa] <- W[,,aa] %*% PWinv
      
      ## Calculate regression coefficients:
      for (a in 1:ncomp) {
        myR <-  matrix(R[,,aa], nrow = npred, ncol = ncomp)
        mytQ <- matrix(tQ[,,aa], nrow = ncomp, ncol = nresp)
        B[,,a,aa] <- myR[,1:a,drop=FALSE] %*% mytQ[1:a,,drop=FALSE]
      }
    } else {
      PWinv <- backsolve(PW, diag(ncomp))
    }
    R[,,aa]   <- W[,,aa] %*% PWinv
    for (a in 1:ncomp) {
      myR  <- matrix( R[,,aa], nrow = npred, ncol = ncomp)
      mytQ <- matrix(tQ[,,aa], nrow = ncomp, ncol = nresp)
      B[,,a,aa] <- myR[,1:a,drop=FALSE] %*% mytQ[1:a,,drop=FALSE]
    }
  }
  
  if (stripped) {
    ## Return as quickly as possible
    list(coefficients = B, Xmeans = Xmeans, Ymeans = Ymeans)
  } else {
    fitted <- fitted + rep(Ymeans, each = nobj) # Add mean
    ## Add dimnames and classes:
    objnames <- dnX[[1]]
    if (is.null(objnames))
      objnames <- dnY[[1]]
    prednames <- dnX[[2]]
    respnames <- dnY[[2]]
    compnames <- paste("Comp", 1:ncomp)
    shrinknames <- paste("Shrink", shrink)
    nCompnames  <- paste(1:ncomp, "comps")
    dimnames(TT) <- dimnames(U) <- list(objnames, compnames, shrinknames)
    dimnames(R)  <- dimnames(W) <- dimnames(P) <-
      list(prednames, compnames, shrinknames)
    dimnames(tQ) <- list(compnames, respnames, shrinknames)
    dimnames(B)  <- list(prednames, respnames, nCompnames, shrinknames)
    dimnames(fitted) <- dimnames(residuals) <-
      list(objnames, respnames, nCompnames, shrinknames)
    class(TT) <- class(U) <- "scores"
    class(P)  <- class(W) <- class(tQ) <- "loadings"
    list(coefficients = B,
         scores = TT, loadings = P,
         loading.weights = W,
         Yscores = U, Yloadings = tQ,
         projection = R,
         Xmeans = Xmeans, Ymeans = Ymeans,
         fitted.values = fitted, residuals = residuals,
         Xvar = colSums(P * P) * tsqs,
         Xtotvar = Xtotvar)
  }
}
