#
#  multinomRob
#
#  Walter R. Mebane, Jr.
#  Cornell University
#  http://macht.arts.cornell.edu/wrm1/
#  wrm1@macht.arts.cornell.edu
#
#  Jasjeet Singh Sekhon 
#  UC Berkeley
#  http://sekhon.polisci.berkeley.edu
#  sekhon@berkeley.edu
#
#  $Id: multinomT.R,v 1.9 2005/09/27 08:04:06 wrm1 Exp $
#
#

##   Yp:  matrix of (overdispersed and contaminated) multinomial proportions
##   Xarray:  array of regressors,
##      dim(Xarray) = c(n observations, n parameters, n categories)
##   xvec:  vector to indicate all the coefficient parameters in the model
##      (parms by ncats):
##      It has a 1 for an estimated parameter and a 0 otherwize.
##      example:
##      > xvec
##           [,1] [,2] [,3] [,4] [,5]
##      [1,]    1    1    1    1    0
##      [2,]    1    1    1    1    0
##      [3,]    1    1    1    1    0
##      [4,]    1    1    1    1    0
##   tvec: parms by ncats matrix (matrix with LQD estimates):
##      example:
##      > tvec
##                          Buchanan        Nader     Gore     Bush Other
##      int               -0.1641034    1.0735560 3.363641 4.151853     0
##      p(r,dg,d,r)96      2.4413780    0.3269827 3.207104 1.676676     0
##      p(cr,g,cd,cr)RV00 15.5333800 1149.4130000 2.039766 1.761392     0
##      pCuban            -6.4083750    0.3546630 1.966287 2.795598     0
##   jacstack:  array of regressors,
##      dim(Xarray) = c(n observations, n UNIQUE parameters, n categories)

multinomT  <- function(Yp, Xarray, xvec, jacstack,
                       start=NA, nobsvec, fixed.df = NA)
  {
    #Y (raw-Y) is assumed to be proportions,
    #the last category to be the contrast
    
    obs    <- dim(Yp)[1];
    cats   <- dim(Yp)[2];
    mcats  <- cats-1;
    nvars  <- dim(Xarray)[2];

    if (any(xvec[,cats] != 0)) {
      stop("(multinomT): invalid specification of Xarray (regressors not allowed for last category");
    }
    smdim <- dim(jacstack);
    stack.index <- matrix(FALSE, smdim[2], smdim[3]);
    for (i in 1:smdim[2]) for (j in 1:smdim[3]) {
      stack.index[i,j] <- !all(jacstack[,i,j] == 0);
    }
    if (sum(stack.index) != sum(xvec != 0)) {
      print("multinomT:  xvec is:"); print(xvec);
      print("multinomT:  stack.index is:"); print(stack.index);
      stop("(multinomT):  jacstack structure check failed");
    }
    kY  <- matrix(nrow=obs,ncol=mcats)
    indx1  <- Yp==0;
    indx1vec  <- apply(indx1,1,sum) > 0;
    if (sum(indx1) > 0)
      {
        cat("multinomT:  Need to remove 0 in multinomT transformation\n");
      }
    Yp[indx1]  <- .5/nobsvec[indx1vec];
    kY  <- log(Yp[,1:mcats]/Yp[,cats])
#    indx  <- is.infinite(kY);

    if (is.na(start[1]))
      {
        mt.obj  <- mt.mle(Xarray=Xarray, xvec=xvec, jacstack=jacstack, y=kY,
                          stack.index=stack.index, fixed.df=fixed.df);
      } else {
        mt.obj  <- mt.mle(Xarray=Xarray, xvec=xvec, jacstack=jacstack, y=kY,
                          stack.index=stack.index, start=start, fixed.df=fixed.df);        
      }

    tvec <- mnl.xvec.mapping(forward=FALSE,xvec,xvec, mt.obj$dp$beta, cats, nvars);
    pred <- mnl.probfunc(Yp, Yp==Yp, Xarray, tvec)
    
    return( list(call=mt.obj$call, logL=mt.obj$logL, deviance=mt.obj$deviance,
           par=mt.obj$dp, se=mt.obj$se, optim=mt.obj$optim, pred=pred))
  }#end of multinomT


#
#this version includes analytical gradients, no bounds on DF, other than df>0
#

mt.mle  <- function (Xarray, xvec, jacstack, y, stack.index,
          start = NA, freq =NA, fixed.df = NA, trace = FALSE,
          method = "BFGS", control = list(maxit = 600,trace=0)) 
{
  nvars  <- dim(Xarray)[2];
  Diag <- function(x) diag(x, nrow = length(x), ncol = length(x))
#  y.name <- deparse(substitute(y))
  y.names <- dimnames(y)[[2]]
  y <- as.matrix(y)
  if (missing(freq) | is.na(freq)) {
#  if (missing(freq)) {    
    freq <- rep(1, nrow(y))
  }
#  x.names <- dimnames(mX)[[2]]

  d <- ncol(y)
  n <- sum(freq)
  m <- sum(xvec == 1) + length(unique(xvec[xvec>1]));

  if (is.na(start[1]))
    {
      cat("mt.mle:  I don't know how to generate starting values in the general case\n");
      stop();
    } #end if
  beta <- start$beta
  Omega <- start$Omega
  if (!is.na(fixed.df)) start$df <- fixed.df;
  df <- start$df

  Oinv <- solve(Omega, tol=1e-100)
  Oinv <- (Oinv + t(Oinv))/2
  upper <- chol(Oinv)
  D <- diag(upper)
  A <- upper/D
  D <- D^2
  if (d > 1)
    {
      param <- c(beta, -0.5 * log(D), A[!lower.tri(A, diag = TRUE)])
    } else {
      param <- c(beta, -0.5 * log(D))
    }
  if (is.na(fixed.df)) 
    param <- c(param, log(df))
  opt <- optim(param, fn = mt.dev, method = method, 
               control = control, hessian = TRUE,
               Xarray=Xarray, xvec=xvec, jacstack=jacstack, y = y,
               stack.index = stack.index, nvars = nvars, freq = freq,
               trace = trace, fixed.df = fixed.df)
  dev <- opt$value
  param <- opt$par
  if (trace) {
    cat("mt.mle:  Message from optimization routine:", opt$message, 
        "\n")
    cat("mt.mle:  deviance:", dev, "\n")
  }
  beta <- param[1:m] ;
  D <- exp(-2 * param[(m + 1):(m + d)])
  if (d > 1) {
    A <- diag(d)
    A[!lower.tri(A, diag = TRUE)] <-
      param[(m + d + 1):(m + d + d * (d - 1)/2)]
    i0 <- m + d + d * (d - 1)/2
  } else {
    i0 <- m + 1
    A <- as.matrix(1)
  }
  if (is.na(fixed.df))
    {
        df <- exp(param[i0 + 1])
      } else {
        df <- fixed.df
      }
  Ainv <- backsolve(A, diag(d))
  Omega <- Ainv %*% Diag(1/D) %*% t(Ainv)
#  omega <- sqrt(diag(Omega))
#  dimnames(beta) <- list(x.names, y.names)
  dimnames(Omega) <- list(y.names, y.names)
  info <- opt$hessian/2
  if (all(is.finite(info))) {
    qr.info <- qr(info)
    info.ok <- (qr.info$rank == length(param))
  } else {
      info.ok <- FALSE
    }
  if (info.ok) {
    se2 <- diag(solve(qr.info))
    if (min(se2) < 0)
      {
        se <- NA
          } else {
            se <- sqrt(se2)
            se.beta <- se[1:m] ;
#            dimnames(se.beta)[2] <- list(y.names)
#            dimnames(se.beta)[1] <- list(x.names)
            se.df <- df * se[i0 + 1]
            se <- list(beta = se.beta, df = se.df, 
                       info = info)
          }
  } else {
    se <- NA
  }
  dp <- list(beta = beta, Omega = Omega, df = df)
  list(call = match.call(), logL = -0.5 * dev, deviance = dev, 
       dp = dp, se = se, optim = opt)
}


mt.dev  <- function (param, Xarray, xvec, jacstack, y, stack.index, nvars, freq,
    fixed.df = NA, trace = FALSE) 
{
    Diag <- function(x) diag(x, nrow = length(x), ncol = length(x))
    d <- ncol(y)
    n <- sum(freq)
    m <- sum(xvec == 1) + length(unique(xvec[xvec>1]));
    beta <- param[1:m];
    D <- exp(-2 * param[(m + 1):(m + d)])
    if (d > 1) {
        A <- diag(d)
        A[!lower.tri(A, diag = TRUE)] <-
          param[(m + d + 1):(m + d + d * (d - 1)/2)]
        i0 <- m + d + d * (d - 1)/2
    } else {
        i0 <- m + 1
        A <- as.matrix(1)
    }
    eta <- rep(0,d);
    if (is.na(fixed.df))
      {
        df <- exp(param[i0 + 1])
      } else {
        df <- fixed.df
      }
    Oinv <- t(A) %*% Diag(D) %*% A
#    cat("beta:\n");
#    print(as.matrix(beta));
    tvec <- mnl.xvec.mapping(forward=FALSE, xvec, xvec, beta, d+1, nvars);
    ylinpred <- y;
    if (dim(tvec)[1] == 1) {
      for (j in 1:d) {
        ylinpred[,j] <- Xarray[,,j] * tvec[,j];
      }
    }
    else {
      for (j in 1:d) {
        ylinpred[,j] <- Xarray[,,j] %*% tvec[,j];
      }
    }
    u <- y - ylinpred ;
#    cat("u:\n")
#    print(u)
#    cat("boo1\n");
    Q <- apply((u %*% Oinv) * u, 1, sum)
    L <- as.vector(u %*% eta)
    logDet <- sum(log(df * pi/D))
    dev <- (n * (2 * lgamma(df/2) + logDet - 2 * lgamma((df + 
        d)/2)) + (df + d) * sum(freq * log(1 + Q/df)) - 2 * sum(freq * 
        log(2 * pt(L * sqrt((df + d)/(Q + df)), df + d))))
    if (trace) 
        cat("mt.dev: ", dev, "\n")
    dev
}


mt.dev.grad  <- function (param, Xarray, xvec, jacstack, y, stack.index, nvars, freq,
                          fixed.df = NA, trace = FALSE) 
{
    Diag <- function(x) diag(x, nrow = length(x), ncol = length(x))
    d <- ncol(y)
    n <- sum(freq)
    m <- sum(xvec == 1) + length(unique(xvec[xvec>1]));
    nvarsunique <- dim(jacstack)[2];
    beta <- param[1:m];    
    D <- exp(-2 * param[(m + 1):(m + d)])
    if (d > 1) {
        A <- diag(d)
        A[!lower.tri(A, diag = TRUE)] <-
          param[(m + d + 1):(m + d + d * (d - 1)/2)]
        i0 <- m + d + d * (d - 1)/2
    }
    else {
        i0 <- m + d
        A <- as.matrix(1)
    }
    eta <- rep(0,d);
    if (is.na(fixed.df)) 
        df <- exp(param[i0 + 1])
    else df <- fixed.df
    tA <- t(A)
    Oinv <- tA %*% Diag(D) %*% A
    tvec <- mnl.xvec.mapping(forward=FALSE, xvec, xvec, beta, d+1, nvars);
    ylinpred <- y;
    if (dim(tvec)[1] == 1) {
      for (j in 1:d) {
        ylinpred[,j] <- Xarray[,,j] * tvec[,j];
      }
    }
    else {
      for (j in 1:d) {
        ylinpred[,j] <- Xarray[,,j] %*% tvec[,j];
      }
    }
    u <- y - ylinpred ;
    Q <- as.vector(apply((u %*% Oinv) * u, 1, sum))
    L <- as.vector(u %*% eta)
    t. <- L * sqrt((df + d)/(Q + df))
    dlogft <- -(df + d)/(2 * df * (1 + Q/df))
#    dt.dL <- sqrt((df + d)/(Q + df))
    dt.dQ <- (-0.5) * L * sqrt(df + d)/(Q + df)^1.5
    T. <- pt(t., df + d)
    dlogT. <- dt(t., df + d)/T.
    u.freq <- u * freq
#    foo1 <- foo2 <- foo3 <- matrix(0, nvarsunique, d+1);
    fooA <- matrix(0, nvarsunique, d);
    for (j in 1:d) {
      fooA[,j] <- t(jacstack[,,j]) %*% (u.freq * (dlogft + dlogT. * dt.dQ));
#      foo1[,j] <- t(jacstack[,,j]) %*% (u.freq * dlogft);
#      foo3[,j] <- t(jacstack[,,j]) %*% (dlogT. * dt.dQ * u.freq);
#      foo2[,j] <- t(jacstack[,,j]) %*% (dlogT. * dt.dL * freq);
    }
#    foo1 <- t(mX) %*% (u.freq * dlogft);
#    foo2 <- t(mX) %*% (dlogT. * dt.dL * freq);
#    foo3 <- t(mX) %*% (dlogT. * dt.dQ * u.freq);
#    Dbeta <- (-2 * (foo1 + foo3) %*% Oinv) - outer(as.vector(foo2), eta)
    Dbeta <- -2 * fooA %*% Oinv
    if (d > 1) {
        M <- 2 * (Diag(D) %*% A %*% t(u * dlogft) %*% u.freq + 
            Diag(D) %*% A %*% t(u * dlogT. * dt.dQ) %*% u.freq)
        DA <- M[!lower.tri(M, diag = TRUE)]
    }
    else DA <- NULL
    M <- (A %*% t(u * dlogft) %*% u.freq %*% tA + A %*% t(u * 
        dlogT. * dt.dQ) %*% u.freq %*% tA)
    if (d > 1) 
        DD <- diag(M) + 0.5 * n/D
    else DD <- as.vector(M + 0.5 * n/D)
    grad <- (-2) * c(Dbeta[stack.index[,-(d+1)]], DD * (-2 * D), DA)
    if (is.na(fixed.df)) {
        dlogft.ddf <- 0.5 * (digamma((df + d)/2) - digamma(df/2) - 
            d/df + (df + d) * Q/((1 + Q/df) * df^2) - log(1 + 
            Q/df))
        eps <- 1e-04
        T.eps <- pt(L * sqrt((df + eps + d)/(Q + df + eps)), 
            df + eps + d)
        dlogT.ddf <- (log(T.eps) - log(T.))/eps
        Ddf <- sum((dlogft.ddf + dlogT.ddf) * freq)
        grad <- c(grad, -2 * Ddf * df)
    }
    if (trace) 
        cat("mt.dev.grad: norm is ", sqrt(sum(grad^2)), "\n")
    return(grad)
}#end of mt.dev.grad
