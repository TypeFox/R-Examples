#' @title Trunction PLS
#'
#' @description Distribution based truncation for variable selection in subspace 
#' methods for multivariate regression.
#'
#' @param ... arguments passed on to \code{mvrV}).
#' @param Y.add optional additional response vector/matrix found in the input data.
#' @param weights optional object weighting vector.
#' @param method choice (default = \code{truncation}).
#'
#' @details Loading weights are truncated around their median based on confidence intervals
#' for modelling without replicates (Lenth et al.).
#'  
#' @return Returns an object of class mvrV, simliar to to mvr object of the pls package.
#'
#' @author Kristian Hovde Liland.
#'
#' @references K.H. Liland, M. Høy, H. Martens, S. Sæbø: Distribution based truncation for 
#' variable selection in subspace methods for multivariate regression, Chemometrics and
#' Intelligent Laboratory Systems 122 (2013) 103-111.
#'
#' @seealso \code{\link{VIP}} (SR/sMC/LW/RC), \code{\link{filterPLSR}}, \code{\link{spa_pls}}, 
#' \code{\link{stpls}}, \code{\link{truncation}}, \code{\link{bve_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{ipw_pls}}, \code{\link{ga_pls}}, \code{\link{rep_pls}}.
#'
#' @examples
#' data(yarn, package = "pls")
#' tr <- truncation(density ~ NIR, ncomp=5, data=yarn, validation="CV",
#'  truncation="Lenth", trunc.width=0.95) # Default truncation
#' summary(tr)
#' 
#' @export
truncation <- function(..., Y.add, weights, method = "truncation"){
  cl <- match.call()
  cl$method <- match.arg(method,c("truncation", "model.frame"))
  cl[[1]] <- quote(mvrV)
  res <- eval(cl, parent.frame())
  ## Fix call component
  if (cl$method != "model.frame") res$call[[1]] <- as.name("truncation")
  if (missing(method)) res$call$method <- NULL
  res
}

trunc_fit <- function(X, Y, ncomp, Y.add = NULL, stripped = FALSE,
                      lower = 0.5, upper = 0.5, trunc.pow = FALSE, weights = NULL,
                      truncation = NULL, trunc.width=NULL, trunc.weight=0, reorth=FALSE, symmetric=FALSE, ...)
{
  ## X       - the data matrix
  ## Y       - the primary response matrix
  ## Y.add   - the additional response matrix (optional)
  ## ncomp   - number of components
  ## lower   - lower bounds for power algorithm (default=0.5)
  ## upper   - upper bounds for power algorithm (default=0.5)
  ## weights - prior weighting of observations (optional)
  
  Yprim <- as.matrix(Y)
  Y <- cbind(Yprim, Y.add)
  
  nobj <- dim(X)[1]
  npred <- dim(X)[2]
  nresp <- dim(Yprim)[2]
  
  if(is.null(weights)){
    Xmeans <- colMeans(X)
    X <- X - rep(Xmeans, each = nobj)
  } else {
    Xmeans <- crossprod(weights,X)/sum(weights)
    X <- X - rep(Xmeans, each = nobj)
  }
  
  X.orig <- X
  Ymeans = colMeans(Yprim)
  
  if(!stripped) {
    ## Save dimnames:
    dnX <- dimnames(X)
    dnY <- dimnames(Yprim)
  }
  dimnames(X) <- dimnames(Y) <- dimnames(Yprim) <- NULL
  
  ## Declaration of variables
  W   <- matrix(0,npred,ncomp)    # W-loadings
  TT  <- matrix(0,nobj,ncomp)     # T-scores
  P   <- matrix(0,npred,ncomp)    # P-loadings
  Q   <- matrix(0,nresp,ncomp)    # Q-loadings
  A   <- matrix(0,dim(Y)[2],ncomp)# Column weights for W0 (from CCA)
  cc  <- numeric(ncomp)
  pot <- rep(0.5,ncomp)           # Powers used to construct the w-s in R
  B   <- array(0, c(npred, nresp, ncomp))
  smallNorm <- numeric()
  if(!stripped){
    U <- TT                     # U-scores
    tsqs <- rep.int(1, ncomp)   # t't
    fitted <- array(0, c(nobj, nresp, ncomp))
  }
  
  for(a in 1:ncomp){
    if( length(lower)==1 && lower==0.5 && length(upper)==1 && upper==0.5 ) {
      Rlist <- Rcal(X, Y, Yprim, weights)} 			   # Default CPLS algorithm
    else {
      Rlist <- RcalP(X, Y, Yprim, weights, lower, upper, trunc.pow) # Alternate CPPLS algorithm
      pot[a] <- Rlist$pot }
    cc[a] <- Rlist$cc
    w.a <- Rlist$w
    ifelse(!is.null(Rlist$a), aa <- Rlist$a, aa <- NA)
    A[,a] <- aa
    
    ## Truncation
    if(!is.null(truncation) && !(is.logical(truncation)&&truncation==FALSE)){
      if(is.null(trunc.width))
        stop("Truncation specified with out 'trunc.width'")
      tdf <- (X%*%w.a)^2
      df <- ceiling(sum(tdf)/max(tdf))
      if(truncation == "Lenth"){
        wWeights <- Lenth(w.a,df, trunc.width, trunc.weight, symmetric)
      } else {
        if(truncation == "quantile"){
          wWeights <- quant(w.a,df, trunc.width, trunc.weight, symmetric)
        } else {
          stop(paste('Unknown truncation type: ', truncation, sep=""))
        }
      }
      wWeights[wWeights<pls.options()$w.tol] <- 0
      w.a <- w.a*wWeights
      w.a <- w.a/norm(w.a)
      if(reorth){
        w.a <- w.a - W[,1:(a-1)]%*%crossprod(W[,1:(a-1)],w.a) # Orthogonalisation
      }
    }
    
    ## Make new vectors orthogonal to old ones?
    ## w.a <- w.a - W[,1:(a-1)]%*%crossprod(W[,1:(a-1)], w.a)
    w.a[abs(w.a) < pls.options()$w.tol] <- 0   # Removes insignificant values
    w.a <- w.a / norm(w.a)                     # Normalization
    t.a <- X %*% w.a                           # Score vectors
    tsq <- crossprod(t.a)[1]
    p.a <- crossprod(X,t.a) / tsq
    q.a <- crossprod(Yprim,t.a) / tsq
    X   <- X - tcrossprod(t.a,p.a)             # Deflation
    
    ## Check and compensate for small norms
    mm <- apply(abs(X),2,sum)
    r <- which(mm < pls.options()$X.tol)
    if(length(r)>0){
      for(i in 1:length(r)){
        if(sum(smallNorm==r[i]) == 0){
          ## Add new short small to list
          smallNorm[length(smallNorm)+1] <- r[i]
        }
      }
    }
    X[,smallNorm] <- 0 # Remove collumns having small norms
    
    W[,a]  <- w.a
    TT[,a] <- t.a
    P[,a]  <- p.a
    Q[,a]  <- q.a
    B[,,a] <- W[,1:a, drop=FALSE]%*%tcrossprod(solve(crossprod(P[,1:a, drop=FALSE],W[,1:a, drop=FALSE])),Q[,1:a, drop=FALSE])
    
    if (!stripped) {
      tsqs[a] <- tsq
      ## Extra step to calculate Y scores:
      U[,a] <- Yprim %*% q.a / crossprod(q.a)[1] # Ok for nresp == 1 ??
      ## make u orth to previous X scores:
      if (a > 1) U[,a] <- U[,a] - TT %*% (crossprod(TT, U[,a]) / tsqs)
      fitted[,,a] <- X.orig %*% B[,,a]
    }
  }
  if (stripped) {
    ## Return as quickly as possible
    list(coefficients = B, Xmeans = Xmeans, Ymeans = Ymeans, gammas = pot)
  } else {
    fitted <- fitted + rep(Ymeans, each = nobj) # Add mean
    residuals <- - fitted + c(Yprim)
    
    ## Add dimnames:
    objnames <- dnX[[1]]
    if (is.null(objnames)) objnames <- dnY[[1]]
    prednames <- dnX[[2]]
    respnames <- dnY[[2]]
    compnames <- paste("Comp", 1:ncomp)
    nCompnames <- paste(1:ncomp, "comps")
    dimnames(TT) <- list(objnames, compnames)
    dimnames(W) <- dimnames(P) <-
      list(prednames, compnames)
    dimnames(Q) <- list(respnames, compnames)
    dimnames(B) <- list(prednames, respnames, nCompnames)
    dimnames(fitted) <- dimnames(residuals) <-
      list(objnames, respnames, nCompnames)
    colnames(A) <- compnames
    class(TT) <- "scores"
    class(P) <- class(W) <- class(Q) <- "loadings"
    
    list(coefficients = B,
         scores = TT, loadings = P,
         loading.weights = W,
         Yscores = U, Yloadings = Q,
         projection = W %*% solve(crossprod(P,W)),
         Xmeans = Xmeans, Ymeans = Ymeans,
         fitted.values = fitted, residuals = residuals,
         Xvar = colSums(P * P) * tsqs,
         Xtotvar = sum(X.orig * X.orig),
         gammas = pot,
         canonical.correlations = cc,
         smallNorm = smallNorm,
         A = A, trunc.pow = trunc.pow)
  }
}

#######################
## Rcal function (CPLS)
Rcal <- function(X, Y, Yprim, weights) {
  W0 <- crossprod(X,Y)
  Ar <- cancorr(X%*%W0, Yprim, weights, FALSE) # Computes canonical correlations between columns in XW and Y with rows weighted according to 'weights'
  w  <- W0 %*% Ar$A[,1, drop=FALSE]  # Optimal loadings
  ifelse(exists('Ar'), a <- Ar$A[,1], a <- NA)
  list(w = w, cc = Ar$r^2, a = a)
}


#########################
## RcalP function (CPPLS)
RcalP <- function(X, Y, Yprim, weights, lower, upper, trunc.pow) {
  CS <- CorrXY(X, Y, weights)     # Matrix of corr(Xj,Yg) and vector of std(Xj)
  sng <- sign(CS$C)               # Signs of C {-1,0,1}
  C <- abs(CS$C)                  # Correlation without signs
  mS <- max(CS$S); S <- CS$S / mS # Divide by largest value
  mC <- max(C); C <- C / mC       #  -------- || --------
  
  ## Computation of the best vector of loadings
  lw <- lw_bestpar(X, S, C, sng, Yprim, weights, lower, upper, trunc.pow)
}


################
## lw_bestpar function
lw_bestpar <- function(X, S, C, sng, Yprim, weights, lower, upper, trunc.pow) {
  if(!is.null(weights))
    weights <- sqrt(weights) # Prepare weights for cca
  # Compute for S and each columns of C the distance from the median scaled to [0,1]
  if(trunc.pow){
    medC <- t(abs(t(sng*C)-apply(sng*C,2,median)))
    medC <- t(t(medC)/apply(medC,2,max))
    medS <- abs(S-median(S))
    medS <- medS/max(medS)
  } else {
    medS <- medC <- NULL
  }
  
  #########################
  # Optimization function #
  #########################
  f <- function(p, X, S, C, sng, Yprim, weights, trunc.pow, medS, medC) {
    if(p == 0){         # Variable selection from standard deviation
      S[S < max(S)] <- 0
      W0 <- S
    } else if(p == 1) { # Variable selection from correlation
      C[C < max(C)] <- 0
      W0 <- rowSums(C)
    } else {            # Standard deviation and correlation with powers
      if(trunc.pow){
        ps <- (1-p)/p
        if(ps<1){
          S <- S^ps
        } else {
          S[medS<(1-2*p)] <- 0
        }
        pc <- p/(1-p)
        if(pc<1){
          W0 <- (sng*(C^pc))*S
        } else {
          C[medC<(2*p-1)] <- 0
          W0 <- (sng*C)*S
        }
      } else {
        S <- S^((1-p)/p)
        W0 <- (sng*(C^(p/(1-p))))*S
      }
    }
    Z <- X %*% W0  # Transform X into W0
    -(cancorr(Z, Yprim, weights))^2
  }
  
  #####################################
  # Logic for optimization segment(s) #
  #####################################
  nOpt <- length(lower)
  pot  <- numeric(3*nOpt)
  ca   <- numeric(3*nOpt)
  
  for (i in 1:nOpt){
    ca[1+(i-1)*3]  <- f(lower[i], X, S, C, sng, Yprim, weights, trunc.pow, medS, medC)
    pot[1+(i-1)*3] <- lower[i]
    if (lower[i] != upper[i]) {
      Pc <- optimize(f = f, interval = c(lower[i], upper[i]),
                     tol = 10^-4, maximum = FALSE,
                     X = X, S = S, C = C, sng = sng, Yprim = Yprim,
                     weights = weights, trunc.pow = trunc.pow, medS, medC)
      pot[2+(i-1)*3] <- Pc[[1]]; ca[2+(i-1)*3] <- Pc[[2]]
    }
    ca[3+(i-1)*3]  <- f(upper[i], X, S, C, sng, Yprim, weights, trunc.pow, medS, medC)
    pot[3+(i-1)*3] <- upper[i]
  }
  
  
  ########################################################
  # Computation of final w-vectors based on optimization #
  ########################################################
  cc <- max(-ca)                      # Determine which is more succesful
  cmin <- which.max(-ca)              # Determine which is more succesful
  if (pot[cmin] == 0) {        # Variable selection from standard deviation
    S[S < max(S)] <- 0
    w <- S
  } else if (pot[cmin] == 1) { # Variable selection from correlation
    C[C < max(C)] <- 0
    w <- rowSums(C)
  } else {                     # Standard deviation and correlation with powers
    p <- pot[cmin]                  # Power from optimization
    if(trunc.pow){ # New power algorithm
      ps <- (1-p)/p
      if(ps<1){
        S <- S^ps
      } else {
        S[medS<(1-2*p)] <- 0
      }
      pc <- p/(1-p)
      if(pc<1){
        W0 <- (sng*(C^pc))*S
      } else {
        C[medC<(2*p-1)] <- 0
        W0 <- (sng*C)*S
      }
    } else {
      S <- S^((1-p)/p)
      W0 <- (sng*(C^(p/(1-p))))*S
    }
    
    Z <- X %*% W0                   # Transform X into W
    Ar <- cancorr(Z, Yprim, weights, FALSE) # Computes canonical correlations between columns in XW and Y with rows weighted according to 'weights'
    w <- W0 %*% Ar$A[,1, drop=FALSE]  # Optimal loadings
  }
  pot <- pot[cmin]
  ifelse(exists('Ar'), a <- Ar$A[,1], a <- NA)
  list(w = w, pot = pot, cc = cc, a = a)
}


################
## CorrXY function
CorrXY <- function(X, Y, weights) {
  ##  Computation of correlations between the columns of X and Y
  n  <- dim(X)[1]
  if(is.null(weights)){
    cx <- colMeans(X)
    cy <- colMeans(Y)
    X <- X - rep(cx, each = n)
    Y <- Y - rep(cy, each = n)
  } else {
    cx <- crossprod(weights,X)/sum(weights)
    cy <- crossprod(weights,Y)/sum(weights)
    X  <- X - rep(cx, each = n)
    Y  <- Y - rep(cy, each = n)
    X  <- X * weights
    Y  <- Y * weights
  }
  
  sdX <- sqrt(apply(X^2,2,mean))
  inds <- which(sdX == 0, arr.ind=FALSE)
  sdX[inds] <- 1
  
  ccxy <- crossprod(X, Y) / (n * tcrossprod(sdX, sqrt(apply(Y^2,2,mean))))
  sdX[inds] <- 0
  ccxy[inds,] <- 0
  list(C = ccxy, S = sdX)
}


################
## function norm
norm <- function(vec) {
  sqrt(crossprod(vec)[1])
}


################
## Stripped version of canonical correlation (cancor)
cancorr <- function (x, y, weights, opt = TRUE) {
  nr  <- nrow(x)
  ncx <- ncol(x)
  ncy <- ncol(y)
  if (!is.null(weights)){
    x <- x * weights
    y <- y * weights
  }
  qx <- qr(x, LAPACK = TRUE)
  qy <- qr(y, LAPACK = TRUE)
  qxR <- qr.R(qx)
  # Compute rank like MATLAB does
  dx <- sum(abs(diag(qxR)) > .Machine$double.eps*2^floor(log2(abs(qxR[1])))*max(nr,ncx))
  if (!dx)
    stop("'x' has rank 0")
  qyR <- qr.R(qy)
  # Compute rank like MATLAB does
  dy <- sum(abs(diag(qyR)) > .Machine$double.eps*2^floor(log2(abs(qyR[1])))*max(nr,ncy))
  if (!dy)
    stop("'y' has rank 0")
  dxy <- min(dx,dy)
  if(opt) {
    z <- svd(qr.qty(qx, qr.qy(qy, diag(1, nr, dy)))[1:dx,, drop = FALSE],
             nu = 0, nv = 0)
    ret <- max(min(z$d[1],1),0)
  } else {
    z <- svd(qr.qty(qx, qr.qy(qy, diag(1, nr, dy)))[1:dx,, drop = FALSE],
             nu = dxy, nv = 0)
    A <- backsolve((qx$qr)[1:dx,1:dx, drop = FALSE], z$u)*sqrt(nr-1)
    if((ncx - nrow(A)) > 0) {
      A <- rbind(A, matrix(0, ncx - nrow(A), dxy))
    }
    A[qx$pivot,] <- A
    ret <- list(A=A, r=max(min(z$d[1],1),0))
  }
  ret
}

# Lenth's method for CI without replicates
Lenth <- function(x, df=0, conf=0.95, weight=0, symmetric=FALSE){
  last <- length(x)
  
  m <- Median(x)
  med <- m$med; sx <- m$sx; i <- m$i; n <- m$n
  if(abs(sx[1])>abs(sx[last])){
    max_out <- i[1]
  } else {
    max_out <- i[last]
  }
  wWeights <- rep(1,n)
  
  if(conf == 0){
    alpha <- 0.975
  } else {
    alpha <- conf+(1-conf)/2
  }
  if(symmetric){
    # Contrast without replicate (RV Lenth 1989)
    c0   <- abs(x); s0  <- 1.5*median(c0)
    PSE  <- 1.5*median(c0[c0<2.5*s0])
    offset <- PSE
    if(df == 0){
      offsets <- rep(1,2)*qt(alpha,1000)*PSE
    } else {
      offsets <- rep(1,2)*qt(alpha,ceiling(df/3))*PSE
    }
  } else {
    offsets <- numeric(2)
    for(j in 1:2){
      if(j==1)
        c0 <- abs(x[x>=0])
      else 
        c0 <- abs(x[x<0])
      s0  <- 1.5*median(c0)
      PSE  <- 1.5*median(c0[c0<2.5*s0])
      if(df == 0){
        offsets[j] <- qt(alpha,1000)*PSE
      } else {
        offsets[j] <- qt(alpha,ceiling(df/3))*PSE
      }
    }
    offset <- min(offsets)
  }
  
  InLiers <- (x<(med+max(offset,offsets[2]))) & (x>(med-max(offset,offsets[1])))
  InLiers[max_out] <- FALSE # At least one point outside
  
  if(weight > 0){ # Weights from normal distribution
    xTr <- x(x<(med+offsets[2]) & x>(med-offsets[1])) # Trim asymmetrically around the median
    std <- sd(xTr)
    wWeights <- 1-2*abs(pt((x-med)/std,df)-0.5)
    #         wWeights = pt((x-med)./sd,df);
    wWeights <- wWeights - min(wWeights)
    wWeights <- 1 - wWeights/max(wWeights)
    if(weight != 1){ # Parameterised sharpening of weights towards cut-off
      m1 <- max(wWeights[InLiers])
      m2 <- min(wWeights[!InLiers])
      weight <- 1-weight/2
      w <- weight/(1-weight)
      if(m1 > 0){
        wWeights[InLiers] <- ((wWeights[InLiers]/m1)^w)*m1
      }
      if(m2 < 1){
        wWeights[!InLiers] <- 1-(((1-wWeights[!InLiers])/(1-m2))^w)*(1-m2)
      }
    }
  } else { # Sharp cut-off, weights = 0 or 1
    wWeights[InLiers]  <- 0
    wWeights[!InLiers] <- 1
  }
  InLiers <- which(InLiers)
  #list(InLiers, wWeights)
  wWeights
}

# Extended median function
Median <- function(x){
  i  <- order(x)
  sx <- x[i]
  n  <- length(x)
  if((n%%2) == 0){
    med <- (sx[n/2] + sx[n/2+1]) / 2
  } else {
    med <- sx[(n+1)/2]
  }
  list(med=med,sx=sx,i=i,n=n)
}

# Quantile line function
quant <- function(x, df=0, quant=0.95, weight=0, symmetric=FALSE){
  last <- length(x)
  
  m <- Median(x)
  med <- m$med; sx <- m$sx; i <- m$i; n <- m$n
  if(abs(sx[1])>abs(sx[last])){
    max_out <- i[1]
  } else {
    max_out <- i[last]
  }
  wWeights <- rep(1,n)
  dev    <- max(sx[n]-med, med-sx[1]) # Maximum deviation from the median
  minDev <- max(abs(med - sx[ceiling(n/2) + (-max(ceiling(n/200),3):max(ceiling(n/200),3) )]))
  
  
  ## Minimise the difference between a quartile line
  # and the distribution (asymmetrically)
  minQQ <- function(tr,trF,direct,sx,y,centers,slope,n1){
    if(direct == 0){
      interv <- tr*c(-1, 1) + centers[1]      # Search both directions
    } else {
      if(direct == -1){
        interv <- c(tr, trF)*c(-1, 1) + centers[1] # Search to the left
      } else {
        interv <- c(trF, tr)*c(-1, 1) + centers[1] # Search to the right
      }
    }
    Y <- centers[2] + slope*(interv-centers[1])
    y <- y[sx>interv[1] & sx<interv[2]]
    x <- sx[sx>interv[1] & sx<interv[2]]
    n <- length(x)
    yi <- approx(interv, Y, x)$y
    mean((yi-y)^2)*n1/n
  }
  
  
  # Adapt straight line in quantile plot
  eprob <- ((1:n) - 0.5)/n
  y <- qt(eprob, df)
  q1x <- quantile(x,0.5-quant); q3x <- quantile(x,0.5+quant)
  q1y <- quantile(y,0.5-quant); q3y <- quantile(y,0.5+quant)
  dx <- q3x - q1x;              dy <- q3y - q1y
  slope <- dy/dx
  centerx <- (q1x + q3x)/2;     centery <- (q1y + q3y)/2
  offset <- optimise(minQQ, interval=c(minDev,dev),  trF=0, direct=0, sx=sx, y=y, centers=c(centerx,centery), slope=slope, n1=n)$minimum # Symmetrically
  if(!symmetric){ # Asymmetrically
    offsetL <- optimise(minQQ, interval=c(offset-.Machine$double.eps,dev),  trF=offset, direct=-1, sx=sx, y=y, centers=c(centerx,centery), slope=slope, n1=n)$minimum
    offsetR <- optimise(minQQ, interval=c(offset-.Machine$double.eps,dev),  trF=offsetL, direct=1, sx=sx, y=y, centers=c(centerx,centery), slope=slope, n1=n)$minimum
    offsets <- optim(c(offsetL,offsetR), minQQ, trF=0, direct=0, sx=sx, y=y, centers=c(centerx,centery), slope=slope, n1=n)$par
  } else {
    offsets <- rep(1,2)*offset
  }
  
  
  InLiers <- (x<(med+max(offset,offsets[2]))) & (x>(med-max(offset,offsets[1])))
  InLiers[max_out] <- FALSE # At least one point outside
  
  if(weight > 0){ # Weights from normal distribution
    xTr <- x(x<(med+offsets[2]) & x>(med-offsets[1])) # Trim asymmetrically around the median
    std <- sd(xTr)
    wWeights <- 1-2*abs(pt((x-med)/std,df)-0.5)
    #         wWeights = pt((x-med)./sd,df);
    wWeights <- wWeights - min(wWeights)
    wWeights <- 1 - wWeights/max(wWeights)
    if(weight != 1){ # Parameterised sharpening of weights towards cut-off
      m1 <- max(wWeights[InLiers])
      m2 <- min(wWeights[!InLiers])
      weight <- 1-weight/2
      w <- weight/(1-weight)
      if(m1 > 0){
        wWeights[InLiers] <- ((wWeights[InLiers]/m1)^w)*m1
      }
      if(m2 < 1){
        wWeights[!InLiers] <- 1-(((1-wWeights[!InLiers])/(1-m2))^w)*(1-m2)
      }
    }
  } else { # Sharp cut-off, weights = 0 or 1
    wWeights[InLiers]  <- 0
    wWeights[!InLiers] <- 1
  }
  InLiers <- which(InLiers)
  #list(InLiers, wWeights)
  wWeights
}
