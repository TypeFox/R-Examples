hierNet <- function(x, y, lam, delta=1e-8, strong=FALSE, diagonal=TRUE, aa=NULL, zz=NULL, center=TRUE, stand.main=TRUE, stand.int=FALSE, 
                    rho=nrow(x), niter=100, sym.eps=1e-3,
                    step=1, maxiter=2000, backtrack=0.2, tol=1e-5,
                    trace=0) {
  # Main Hiernet function for fitting at a single parameter lambda.
  # Note: L1 penalty terms have parameter lam.l1 = lambda * (1-delta)
  #       and L2 penalty has parameter lam.l2 = lambda * delta.
  #
  # stand.main and stand.int refer to scaling
  stopifnot(nrow(x) == length(y), lam >= 0, delta >= 0, delta <= 1)
  stopifnot(!is.null(step) && !is.null(maxiter))
  if (strong) stopifnot(!is.null(niter))
  stopifnot(class(y) == "numeric")
  stopifnot(class(lam) == "numeric")
  stopifnot(class(delta) == "numeric")
  stopifnot(class(step) == "numeric", step > 0, maxiter > 0)
  stopifnot(is.finite(x), is.finite(y), is.finite(lam), is.finite(delta))
  this.call <- match.call()

  if (!center) cat("WARNING: center=FALSE should almost never be used.  This option is available for special uses only.", fill=TRUE)
  # center and (maybe) scale variables
  x <- scale(x, center=center, scale=stand.main)
  mx <- attr(x, "scaled:center")
  sx <- attr(x, "scaled:scale") # may be NULL
  if (center) {
    my <- mean(y)
    y <- y - my
  } else my <- NULL
  
  if (is.null(zz)) {
    if (trace > 0) cat("Computing zz...", fill=TRUE)
    zz <- compute.interactions.c(x, diagonal=diagonal)
  }
  if (is.matrix(zz)) {
    zz <- scale(zz, center=center, scale=stand.int)
    mzz <- attr(zz, "scaled:center")
    szz <- attr(zz, "scaled:scale") # may be NULL
    zz <- as.numeric(zz)
  } else {
    mzz <- szz <- NULL
    #cat("Provided zz is not a matrix, so it's assumed to be already centered.", fill=TRUE)
  }

  xnum <- as.numeric(x)
  p <- ncol(x)
  lam.l1 <- lam * (1 - delta)
  lam.l2 <- lam * delta
  if (strong) {
    # strong hierarchy -- use ADMM4
    if (is.null(rho)) rho <- as.numeric(nrow(x))
    stopifnot(is.numeric(rho), is.finite(rho))
    aa <- admm4(x, xnum, y, lam.l1=lam.l1, lam.l2=lam.l2, diagonal=diagonal, zz=zz,
                 rho=rho, niter=niter, aa=aa, sym.eps=sym.eps, # ADMM params
                 stepsize=step, backtrack=backtrack, maxiter=maxiter, tol=tol, # GG params
                 trace=trace)
    # lack of symmetry in theta means that sometimes strong hierarchy will be (very slightly violated)
    ii <- aa$bp + aa$bn == 0
    # note aa$th[ii, ] = 0 since weak hierarchy holds for sure
    if (sum(ii) > 0 & sum(ii) < p) {
      thr <- max(abs(aa$th[!ii, ii]))
      if (thr > 0) {
        cat("  thr = ",thr, fill=TRUE)
        if (thr > 1e-3)
          warning("Had to change ADMM's 'th' by more than 0.001 to make strong hier hold! Increase niter (and/or rho). ")
        aa$th[abs(aa$th) <= thr] <- 0
      }
    }
  } else {
    # weak hierarchy -- a single call to generalized gradient descent
    if (is.null(aa)) {
      aa <- list(th=matrix(0, p, p), bp=rep(0, p), bn=rep(0, p))
    } else {
      stopifnot(dim(aa$th) == c(p,p), length(aa$bp) == p, length(aa$bn) == p)
    }
    # this could be improved by not actually creating V...
    V <- matrix(0, p, p)
    rho <- 0
    aa <- ggdescent.c(x=x, xnum=xnum, zz=zz, y=y, lam.l1=lam.l1, lam.l2=lam.l2, diagonal=diagonal,
       	  	      rho=rho, V=V,
                      stepsize=step, backtrack=backtrack, maxiter=maxiter, tol=tol,
                      aa=aa, trace=trace)
  }

  aa$lam <- lam
  aa$delta <- delta
  aa$type <- "gaussian"
  aa$diagonal <- diagonal
  aa$strong <- strong
  aa$obj <- Objective(aa=aa, x=x, y=y, lam.l1=lam.l1, lam.l2=lam.l2, xnum=xnum, zz=zz, strong=strong)
  aa$step <- step
  aa$maxiter <- maxiter
  aa$backtrack <- backtrack
  aa$tol <- tol
  if (strong) {
    # ADMM parameters:
    aa$rho <- rho
    aa$niter <- niter
    aa$sym.eps <- sym.eps
  }
  aa$mx <- mx
  aa$sx <- sx
  aa$my <- my
  aa$mzz <- mzz
  aa$szz <- szz
  aa$call <- this.call
  class(aa) <- "hierNet"
  return(aa)
}

print.hierNet <- function(x, ...) {
  cat("Call:\n")
  dput(x$call)
  th=(x$th+t(x$th))/2
  o2=colSums(th^2)!=0
  b=x$bp-x$bn
  o=b!=0
  b=b[o]
  if (any(o2)) {
    # model has interactions
    th=th[o,o2,drop=FALSE]
    tight <- rowSums(abs(th)) >= x$bp[o] + x$bn[o] - 1e-9
    tt <- rep("", length(tight))
    tt[tight] <- "*"
    mat=cbind(b,th)
    mat=round(mat,4)
    mat <- cbind(mat, tt)
    cat("\n")
    cat("Non-zero coefficients:",fill=T)
    cat("  (Rows are predictors with nonzero main effects)",fill=T)
    cat("  (1st column is main effect)", fill=T)
    cat("  (Next columns are nonzero interactions of row predictor)", fill=T)
    cat("  (Last column indicates whether hierarchy constraint is tight.)",fill=T)
    cat("\n")
    dimnames(mat)=list(as.character(which(o)),c("Main effect",as.character(which(o2)),"Tight?"))
    print(mat, quote = FALSE)
  } else {
    mat <- matrix(round(b,4), length(b), 1)
    cat("\n")
    cat("Non-zero coefficients:",fill=T)
    cat("  (No interactions in this model)",fill=T)
    cat("\n")
    dimnames(mat)=list(as.character(which(o)),"Main effect")
    print(mat, quote = FALSE)    
  }
  invisible()
}

print.hierNet.path <- function(x, ...) {
  cat("Call:\n")
  dput(x$call)
  b=x$bp-x$bn
  mat=cbind(round(x$lam,2),round(x$obj,2),colSums(b!=0),apply(x$th!=0,3,function(a) sum(diag(a)) + sum((a+t(a)!=0)[upper.tri(a)])))

  dimnames(mat)=list(NULL,c("Lambda", "Objective", "Number of main effects","Number of interactions"))
  cat("\n")
  print(mat, quote = FALSE)
  invisible()
}

print.hierNet.cv <- function(x, ...) {
  cat("Call:\n")
  dput(x$call)
mat=cbind(round(x$lamlist,2),x$nonzero,round(x$cv.err,2),round(x$cv.se,2))

  dimnames(mat)=list(NULL,c("Lambda", "Number of nonzero","Mean CV error", "SE"))
 cat("\n")
  print(mat, quote = FALSE)
 cat("\n")
cat(c("lamhat=",round(x$lamhat,2),"lamhat.1se=",round(x$lamhat.1se,2)),fill=T)

  invisible()
}


hierNet.path <- function(x, y, lamlist=NULL, delta=1e-8, minlam=NULL, maxlam=NULL, nlam=20, flmin=.01,
                         diagonal=TRUE, strong=FALSE, aa=NULL, zz=NULL,
                         stand.main=TRUE, stand.int=FALSE,
                         rho=nrow(x), niter=100, sym.eps=1e-3,# ADMM params
                         step=1, maxiter=2000, backtrack=0.2, tol=1e-5, # GG descent params
                         trace=0) {
  # Main Hiernet function for fitting at a sequence of lambda values.
  # Note: L1 penalty terms have parameter lam.l1 = lambda * (1-delta)
  #       and L2 penalty has parameter lam.l2 = lambda * delta.
  #
  # Always centers both x and zz (unless zz is provided in as.numeric form)
  # stand.main and stand.int refer to whether main effects and interactions should have norm sqrt(n-1)
  
  # center and (maybe) scale variables
  this.call <- match.call()
  x <- scale(x, center=TRUE, scale=stand.main)
  mx <- attr(x, "scaled:center")
  sx <- attr(x, "scaled:scale") # may be NULL
  my <- mean(y)
  y <- y - my
  
  if (is.null(maxlam)) {
    if (!is.null(minlam)) stop("Cannot have maxlam=NULL if minlam is non-null.")
  #  maxlam <- max(abs(t(x) %*% y)/colSums(x^2))
    maxlam <- max(abs(t(x) %*% y))
  #  temp <- t(scale(t(x), center=FALSE, scale=1/y))
  #  temp2 <- apply(temp, 2, twonorm)
  #  maxlam <- max(max(temp2), maxlam)
    minlam <- maxlam * flmin
  }
  if (is.null(minlam)) minlam <- maxlam * flmin
  if (is.null(lamlist))
    lamlist <- exp(seq(log(maxlam),log(minlam),length=nlam))
  nlam <- length(lamlist)
  
  if (is.null(zz))
    zz <- compute.interactions.c(x, diagonal=diagonal)
  else
    stopifnot(is.matrix(zz))
  
  # center and (maybe) scale zz
  zz <- scale(zz, center=TRUE, scale=stand.int)
  mzz <- attr(zz, "scaled:center")
  szz <- attr(zz, "scaled:scale") # may be NULL
  
  zz <- as.numeric(zz)
  p <- ncol(x)
  cp2 <- choose(p, 2)
  bp <- bn <- matrix(NA, nrow=p, ncol=nlam)
  th <- array(NA, c(p, p, nlam))
  obj <- rep(NA, nlam)
  aa <- NULL
  for (i in seq(nlam)) {
    cat(c("i,lam=", i, round(lamlist[i],2)), fill=TRUE)
    aa <- hierNet(x, y, lam=lamlist[i], delta=delta, strong=strong, diagonal=diagonal, aa=aa, zz=zz,
                  stand.main=FALSE, stand.int=FALSE, # have already standardized
                  rho=rho, niter=niter, sym.eps=sym.eps,
                  step=step, maxiter=maxiter, backtrack=backtrack, tol=tol, trace=trace)
    bp[, i] <- aa$bp
    bn[, i] <- aa$bn
    th[, , i] <- aa$th
    obj[i] <- aa$obj
  }
  dimnames(bp) <- dimnames(bn) <- list(as.character(1:p), NULL)
  dimnames(th) <- list(as.character(1:p), as.character(1:p), NULL)

  out <- list(bp=bp, bn=bn, th=th, obj=obj, lamlist=lamlist, delta=delta, mx=mx, sx=sx, mzz=mzz, szz=szz, my=my,
              type="gaussian", diagonal=diagonal, strong=strong,
              step=step, maxiter=maxiter, backtrack=backtrack, tol=tol,    
              call=this.call)
  if (strong) {
    # ADMM parameters:
    out$rho <- rho
    out$niter <- niter
    out$sym.eps <- sym.eps
  }
  class(out) <- "hierNet.path"
  out
}

predict.hierNet <- function(object, newx, newzz=NULL, ...) {
  n <- nrow(newx)
  if (is.null(object$sx))
    newx <- scale(newx, center=object$mx, scale=FALSE)
  else
    newx <- scale(newx, center=object$mx, scale=object$sx)    
  if (is.null(newzz))
    newzz <- compute.interactions.c(newx, diagonal=object$diagonal)
  if (is.null(object$szz))
    newzz <- scale(newzz, center=object$mzz, scale=FALSE)
  else
    newzz <- scale(newzz, center=object$mzz, scale=object$szz)
  newzz <- as.numeric(newzz)
  newx <- as.numeric(newx)
  stopifnot(is.finite(newzz), is.finite(newx))
  if (class(object$bp) == "numeric")
    yhatt <- Compute.yhat.c(newx, newzz, object) + object$my
  else {
    nlam <- ncol(object$bp)
    yhat <- matrix(NA, n, nlam)
    # this could be made more efficient
    for (i in seq(nlam)) {
      bb <- list(bp=object$bp[, i], bn=object$bn[, i], th=object$th[, , i], diagonal=object$diagonal)
      yhat[, i] <- Compute.yhat.c(newx, newzz, bb)
    }
    yhatt <- yhat + object$my
  }
  if (object$type == "logistic") {
    # predict from hierNet.logistic object object
    b0 <- object$b0
    if(is.matrix(yhatt))
      b0 <- matrix(b0, nrow=nrow(yhatt), ncol=ncol(yhatt), byrow=T)
    yhatt <- b0 + yhatt
    pr <- 1 / (1 + exp(-yhatt))
    return(list(prob=pr, yhat=1*(pr>.5)))
  }
  return(yhatt)
}

predict.hierNet.path <- function(object, newx, newzz=NULL, ...){
 predict.hierNet(object, newx, newzz, ...)
}

admm4 <- function(x, xnum, y, lam.l1, lam.l2, diagonal, zz=NULL, rho, niter, aa=NULL, sym.eps=1e-3, trace=1, ...) {
  # Performs ADMM4.
  # Note: xnum is the matrix x as a numeric.  Both are passed to avoid having to call as.numeric too
  # many times.
  p <- ncol(x)
  if (is.null(zz)) {
    if (trace > 0) cat("Computing zz...", fill=TRUE)
    zz <- as.numeric(compute.interactions.c(x, diagonal=diagonal))
  }
  else if (class(zz) == "matrix") zz <- as.numeric(zz)
  if (is.null(aa)) {
    aa <- list(u=matrix(0, p, p),
               th=matrix(0, p, p),
               bp=rep(0, p),
               bn=rep(0, p),
               tt=matrix(0, p, p),
               diagonal=diagonal)
  } else {
    stopifnot(diagonal == aa$diagonal)
  }
  if (is.null(aa$tt) || is.null(aa$u)) {
    aa$tt <- 0.5 * (aa$th + t(aa$th))
    aa$u <- matrix(0, p, p)
  }
  obj <- Objective(aa=aa, x=x, y=y, lam.l1=lam.l1, lam.l2=lam.l2, xnum=xnum, zz=zz, strong=TRUE, sym.eps=sym.eps)
  ll <- NULL
  for (i in seq(niter)) {
    if (trace > 0) cat(i, " ")
    ll <- c(ll, ADMM4.Lagrangian(aa, xnum, zz, y, lam.l1=lam.l1, lam.l2=lam.l2, diagonal=diagonal, rho))
    V <- aa$u - rho * aa$tt
    gg <- ggdescent.c(x, xnum, zz, y, lam.l1=lam.l1, lam.l2=lam.l2, diagonal=diagonal, rho, V, trace=trace-1, aa=aa, ...)
    aa$th <- gg$th
    aa$bp <- gg$bp
    aa$bn <- gg$bn
    aa$tt <- (aa$th + t(aa$th)) / 2 + (aa$u + t(aa$u)) / (2 * rho)
    aa$u <- aa$u + rho * (aa$th - aa$tt)
    obj <- c(obj, Objective(aa=aa, x=x, y=y, lam.l1=lam.l1, lam.l2=lam.l2, xnum=xnum, zz=zz, strong=TRUE, sym.eps=sym.eps))
    if (trace > 0) cat(obj[i+1], fill=TRUE)
  }
  if (max(abs(aa$th-t(aa$th))) > sym.eps)
    cat("Attention: th not symmetric within the desired sym.eps.  Run ADMM for more iterations. And try increasing rho.")
  aa$obj <- obj
  aa$lagr <- ll
  aa
}

Objective <- function(aa, x, y, lam.l1, lam.l2, xnum=NULL, zz=NULL, strong=TRUE, sym.eps=1e-3) {
  # evaluates the NewYal objective at aa.
  if (strong) {
    if (max(aa$th-t(aa$th)) > sym.eps) {
      cat("Theta is not symmetric.", fill=TRUE)
      return(Inf)
    }
  }
  if (any(rowSums(abs(aa$th)) > aa$bp + aa$bn + 1e-5)) {
    cat("hierarchy violated.", fill=TRUE)
    return(Inf)
  }
  if (any(aa$bp < -1e-5)||any(aa$bn < -1e-5)) {
    cat("Non-negative of bp or bn violated.", fill=TRUE)
    return(Inf)
  }
  if (aa$diagonal == FALSE)
    if (any(abs(diag(aa$th)) > 1e-8)) {
      cat("Zero diagonal violated.", fill=TRUE)
      return(Inf)
    }
  if (is.null(zz)) {
    zz <- as.numeric(compute.interactions.c(x, diagonal=aa$diagonal))
  }
  if (is.null(xnum)) xnum <- as.numeric(x)
  r <- y - Compute.yhat.c(xnum, zz, aa)
  pen <- lam.l1 * sum(aa$bp + aa$bn) + lam.l1 * sum(abs(aa$th))/2 + lam.l1 * sum(abs(diag(aa$th)))/2
  pen <- pen + lam.l2 * (sum(aa$bp^2) + sum(aa$bn^2) + sum(aa$th^2))
  sum(r^2)/2 + pen
}

Objective.logistic <- function(aa, x, y, lam.l1, lam.l2, xnum=NULL, zz=NULL, strong=TRUE, sym.eps=1e-3) {
  # evaluates the logistic hiernet objective at aa.
  stopifnot(y %in% c(0,1))
  stopifnot("diagonal" %in% names(aa))
  if (aa$diagonal == FALSE)
    if (any(abs(diag(aa$th)) > 1e-8)) {
      cat("Diagonal of Theta is nonzero.", fill=TRUE)
      return(Inf)
    }
  if (strong) {
    if (max(aa$th-t(aa$th)) > sym.eps) {
      cat("Theta is not symmetric.", fill=TRUE)
      return(Inf)
    }
  }
  if (any(rowSums(abs(aa$th)) > aa$bp + aa$bn + 1e-5)) {
    cat("hierarchy violated.", fill=TRUE)
    return(Inf)
  }
  if (any(aa$bp < -1e5)||any(aa$bn < -1e5)) {
    cat("Non-negative of bp or bn violated.", fill=TRUE)
    return(Inf)
  }
  if (is.null(zz)) {
    zz <- as.numeric(scale(compute.interactions.c(x, diagonal=aa$diagonal), center=TRUE, scale=FALSE))
  }
  if (is.matrix(zz)) zz <- as.numeric(zz)
  if (is.null(xnum)) xnum <- as.numeric(x)
  phat <- Compute.phat.c(xnum, zz, aa)
  loss <- -sum(y*log(phat)) - sum((1-y)*log(1-phat))
  pen <- lam.l1 * sum(aa$bp + aa$bn) + lam.l1 * sum(abs(aa$th))/2 + lam.l1 * sum(abs(diag(aa$th)))/2
  pen <- pen + lam.l2 * (sum(aa$bp^2) + sum(aa$bn^2) + sum(aa$th^2))
  loss + pen
}

compute.interactions.c <- function(x, diagonal=TRUE) {
  # Returns (uncentered) n by cp2 matrix of interactions.
  # The columns of zz are in standard order (11), 12,13,14,...,(22),23,...
  # z's (jk)th column is x_j * x_k
  n <- nrow(x)
  p <- ncol(x)
  cp2 <- p * (p - 1) / 2
  if (diagonal) {
    cp2 <- cp2 + p
    out <- .C("ComputeInteractionsWithDiagWithIndices",
              as.double(x),
              as.integer(n),
              as.integer(p),
              z=rep(0, n * cp2),
              i1=as.integer(rep(0, cp2)),
              i2=as.integer(rep(0, cp2)), dupl = FALSE, PACKAGE="hierNet")
  }
  else {
    out <- .C("ComputeInteractionsWithIndices",
              as.double(x),
              as.integer(n),
              as.integer(p),
              z=rep(0, n * cp2),
              i1=as.integer(rep(0, cp2)),
              i2=as.integer(rep(0, cp2)), dupl = FALSE, PACKAGE="hierNet")
  }
  z <- matrix(out$z, n, cp2)
  rownames(z) <- rownames(x)
  if (is.null(colnames(x))) {
    colnames(z) <- paste(out$i1, out$i2, sep=":")
  }
  else {
    colnames(z) <- paste(colnames(x)[out$i1], colnames(x)[out$i2], sep=":")
  }
  z
}

compute.full.interactions.c <- function(x) {
  # Returns (uncentered) n by p^2 matrix of interactions.
  # The columns of zz are in standard order 11,12,13,14,...,23,...
  # z's (jk)th column is x_j * x_k
  n <- nrow(x)
  p <- ncol(x)
  out <- .C("ComputeFullInteractions",
            as.double(x),
            as.integer(n),
            as.integer(p),
            z=rep(0, n * p^2),
            dupl = FALSE, PACKAGE="hierNet")
  matrix(out$z, n, p^2)
}


Compute.yhat.c <- function(xnum, zz, aa) {
  # aa: list containing bp, bn, th, diagonal
  # note: zz is the n by cp2 matrix, whereas z is the n by p^2 one.
  p <- length(aa$bp)
  n <- length(xnum) / p
  stopifnot(n==round(n))
  stopifnot("diagonal" %in% names(aa))
  if (aa$diagonal) stopifnot(length(zz) == n * (choose(p,2) + p))
  else stopifnot(length(zz) == n * choose(p,2))

  out <- .C("compute_yhat_zz_R",
            xnum,
            as.integer(n),
            as.integer(p),
            zz,
            as.integer(aa$diagonal),
            as.double(aa$th),
            aa$bp,
            aa$bn,
            yhat=rep(0, n),
            DUP=FALSE, PACKAGE="hierNet")
  out$yhat
}

Compute.phat.c <- function(xnum, zz, aa) {
  # aa: list containing b0, bp, bn, th
  # note: zz is the n by cp2 matrix, whereas z is the n by p^2 one.
  stopifnot(c("b0","bp","bn","th","diagonal") %in% names(aa))
  p <- length(aa$bp)
  n <- length(xnum) / p
  if (is.matrix(xnum)) xnum <- as.numeric(xnum)
  stopifnot(n == round(n))
  if (aa$diagonal) stopifnot(length(zz) == n * (choose(p,2) + p))
  else stopifnot(length(zz) == n * choose(p,2))
  #void compute_phat_zz_R(double *x, int *n, int *p, double *zz, int *diagonal,
  #		       double *b0, double *th, double *bp, double *bn, double *phat) {
  out <- .C("compute_phat_zz_R",
            xnum,
            as.integer(n),
            as.integer(p),
            zz,
            as.integer(aa$diagonal),
            as.double(aa$b0),
            as.double(aa$th),
            aa$bp,
            aa$bn,
            phat=rep(0, n),
            DUP=FALSE, PACKAGE="hierNet")
  out$phat
}

ggdescent.c <- function(x, xnum, zz, y, lam.l1, lam.l2, diagonal, rho, V, stepsize, backtrack=0.2, maxiter=100,
                             tol=1e-5, aa=NULL, trace=1) {
  # See ADMM4 pdf for the problem this solves.
  # 
  # x, xnum, zz, y: data (note: zz is a length n*cp2 vector, not a matrix) xnum is x as a vector
  # lam.l1: l1-penalty parameter
  # lam.l2: l2-penalty parameter
  # rho: admm parameter
  # V: see ADMM4 pdf
  # stepsize: step size to start backtracking with
  # backtrack: factor by which step is reduced on each backtrack.
  # maxiter: number of generalized gradient steps to take.
  # tol: stop gg descent if change in objective is below tol.
  # aa: initial estimate of (th, bp, bn)
  # trace: how verbose to be
  #
  # void ggdescent_R(double *x, int *n, int *p, double *zz, int *diagonal, double *y, 
  #		     double *lamL1, double*lamL2, double *rho, double *V, int *maxiter, 
  #		     double *curth, double *curbp, double *curbn,
  #		     double *t, int *stepwindow, double *backtrack, double *tol, int *trace,
  #		     double *th, double *bp, double *bn) {
  n <- length(y)
  p <- ncol(x)
  stepwindow <- 10
  if (is.null(aa)) aa <- list(th=matrix(0,p,p), bp=rep(0,p), bn=rep(0,p))
  out <- .C("ggdescent_R",
            xnum,
            as.integer(n),
            as.integer(p),
            zz,
            as.integer(diagonal),
            y,
            as.double(lam.l1),
	    as.double(lam.l2),
            as.double(rho),
            as.double(V),
            as.integer(maxiter),
            as.double(aa$th),
            aa$bp,
            aa$bn,
            stepsize,
            as.integer(stepwindow),
            backtrack,
            tol,
            as.integer(trace),
            th=rep(0, p*p),
            bp=rep(0, p),
            bn=rep(0, p),
            DUP=FALSE, PACKAGE="hierNet")
  list(bp=out$bp, bn=out$bn, th=matrix(out$th, p, p))
}


hierNet.logistic <- function(x, y, lam, delta=1e-8, diagonal=TRUE, strong=FALSE, aa=NULL, zz=NULL, center=TRUE,
                             stand.main=TRUE, stand.int=FALSE,
                             rho=nrow(x), niter=100, sym.eps=1e-3,# ADMM params
                             step=1, maxiter=2000, backtrack=0.2, tol=1e-5, # GG descent params
                             trace=1) {
  # Solves the logistic regression hiernet. Returns (b0, bp, bn, th)
  this.call <- match.call()
  n <- nrow(x)
  p <- ncol(x)
  stopifnot(y %in% c(0,1))
  stopifnot(length(y) == n, lam >= 0, delta >= 0, delta <= 1)
  stopifnot(!is.null(step) && !is.null(maxiter))
  stopifnot(class(lam) == "numeric")
  stopifnot(class(delta) == "numeric")
  stopifnot(class(step) == "numeric", step > 0, maxiter > 0)
  stopifnot(is.finite(x), is.finite(y), is.finite(lam), is.finite(delta))
  lam.l1 <- lam * (1 - delta)
  lam.l2 <- lam * delta
  if (!center) 
    cat("WARNING: center=FALSE should almost never be used.  This option is available for special uses only.", fill = TRUE)
  x <- scale(x, center = center, scale = stand.main)
  mx <- attr(x, "scaled:center")
  sx <- attr(x, "scaled:scale")
  if (is.null(aa)) aa <- list(b0=0, bp=rep(0, p), bn=rep(0, p), th=matrix(0, p, p), diagonal=diagonal)
  if (is.null(zz)) {
    if (trace > 0) cat("Computing zz...", fill=TRUE)
    zz <- compute.interactions.c(x, diagonal=diagonal)
  }
  if (is.matrix(zz)) {
    zz <- scale(zz, center=center, scale=stand.int)
    mzz <- attr(zz, "scaled:center")
    szz <- attr(zz, "scaled:scale")
    zz <- as.numeric(zz)
  } else {
    mzz <- szz <- NULL
    #cat("Provided zz is not a matrix, so it's assumed to be already centered.", fill = TRUE)
  }
  xnum <- as.numeric(x)
  if (strong) {
    # strong hierarchy -- use ADMM4 (logistic regression version)
    stopifnot(is.numeric(rho), is.finite(rho))
    out <- admm4.logistic(x, xnum, y, lam.l1, lam.l2, diagonal=diagonal, zz=zz,
                          rho=rho, niter=niter, aa=aa, sym.eps=sym.eps, # ADMM params
                          stepsize=step, backtrack=backtrack, maxiter=maxiter, tol=tol, # GG params
                          trace=trace)
    ii <- out$bp + out$bn == 0
    # note out$th[ii, ] = 0 since weak hierarchy holds for sure
    sumii <- sum(ii)
    if (sumii > 0 && sumii < p) {
      thr <- max(abs(out$th[!ii, ii]))
      if (thr > 0) {
        cat("  thr = ",thr, fill=TRUE)
        if (thr > 1e-3)
          warning("Had to change ADMM's 'th' by more than 0.001 to make strong hier hold! Increase niter (and/or rho). ")
        aa$th[abs(aa$th) <= thr] <- 0
      }
    }
  } else {
    out <- ggdescent.logistic(xnum=xnum, zz=zz, y=y, lam.l1=lam.l1, lam.l2=lam.l2, diagonal=diagonal, rho=0, V=matrix(0,p,p),
                                  stepsize=step, backtrack=backtrack, maxiter=maxiter,
                                  tol=tol, aa=aa, trace=trace)
  }
  out$call <- this.call
  out$lam <- lam
  out$delta <- delta
  out$type <- "logistic"
  out$diagonal <- diagonal
  out$strong <- strong
  if (strong) {
    # ADMM parameters:
    out$rho <- rho
    out$niter <- niter
    out$sym.eps <- sym.eps
  }
  out$step <- step
  out$maxiter <- maxiter
  out$backtrack <- backtrack
  out$tol <- tol
  out$obj <- critf.logistic(x, y, lam.l1, lam.l2, out$b0, out$bp, out$bn, out$th)
  out$mx <- mx
  out$my <- 0
  out$sx <- sx
  out$mzz <- mzz
  class(out) <- "hierNet"
  return(out)
}


admm4.logistic <- function(x, xnum, y, lam.l1, lam.l2, diagonal, zz=NULL, rho=10, niter, aa=NULL, sym.eps=1e-3, trace=1, ...) {
  # Performs ADMM4 for logistic loss.
  # Note: xnum is the matrix x as a numeric.  Both are passed to avoid having to call as.numeric too
  # many times.
  p <- ncol(x)
  if (is.null(zz)) {
    if (trace > 0) cat("Computing zz...", fill=TRUE)
    zz <- as.numeric(compute.interactions.c(x, diagonal=diagonal))
  }
  else if (class(zz) == "matrix") zz <- as.numeric(zz)
  if (is.null(aa)) {
    aa <- list(u=matrix(0, p, p),
               th=matrix(0, p, p),
               bp=rep(0, p),
               bn=rep(0, p),
               tt=matrix(0, p, p),
               diagonal=diagonal)
  }
  if (is.null(aa$tt) || is.null(aa$u)) {
    aa$tt <- 0.5 * (aa$th + t(aa$th))
    aa$u <- matrix(0, p, p)
  }
  obj <- Objective.logistic(aa=aa, x=x, y=y, lam.l1=lam.l1, lam.l2=lam.l2, xnum=xnum, zz=zz, strong=TRUE, sym.eps=sym.eps)
  for (i in seq(niter)) {
    if (trace > 0) cat(i, " ")
    V <- aa$u - rho * aa$tt
    gg <- ggdescent.logistic(xnum, zz, y, lam.l1=lam.l1, lam.l2=lam.l2, diagonal=diagonal, rho, V, trace=trace-1, aa=aa, ...)
    aa$th <- gg$th
    aa$bp <- gg$bp
    aa$bn <- gg$bn
    aa$tt <- (aa$th + t(aa$th)) / 2 + (aa$u + t(aa$u)) / (2 * rho)
    aa$u <- aa$u + rho * (aa$th - aa$tt)
    obj <- c(obj, Objective.logistic(aa=aa, x=x, y=y, lam.l1=lam.l1, lam.l2=lam.l2, xnum=xnum, zz=zz, strong=TRUE, sym.eps=sym.eps))
    if (trace > 0) cat(obj[i+1], fill=TRUE)
  }
  if (max(abs(aa$th-t(aa$th))) > sym.eps)
    cat("Attention: th not symmetric within the desired sym.eps.  Run ADMM for more iterations. And try increasing rho.")
  aa$obj <- obj
  aa
}



ggdescent.logistic <- function(xnum, zz, y, lam.l1, lam.l2, diagonal, rho, V, stepsize, backtrack=0.2, maxiter=100,
                             tol=1e-5, aa=NULL, trace=1) {
  # See ADMM4 pdf and logistic.pdf for the problem this solves.
  # 
  # xnum, zz, y: data (note: zz is a length n*cp2 vector, not a matrix) xnum is x as a (n*p)-vector
  # lam.l1: l1-penalty parameter
  # lam.l2: l2-penalty parameter
  # rho: admm parameter
  # V: see ADMM4 pdf
  # stepsize: step size to start backtracking with
  # backtrack: factor by which step is reduced on each backtrack.
  # maxiter: number of generalized gradient steps to take.
  # tol: stop gg descent if change in objective is below tol.
  # aa: initial estimate of (b0, th, bp, bn)
  # trace: how verbose to be
  #
  #void ggdescent_logistic_R(double *x, int *n, int *p, double *zz, int * diagonal, double *y, 
  #			     double *lamL1, double *lamL2, double *rho, double *V, int *maxiter, 
  #			     double *curb0, double *curth, double *curbp, double *curbn,
  #			     double *t, int *stepwindow, double *backtrack, double *tol, int *trace,
  #			     double *b0, double *th, double *bp, double *bn) {
  n <- length(y)
  p <- length(xnum) / n
  stopifnot(p == round(p))
  if (diagonal) stopifnot(length(zz) == n * (choose(p,2)+p))
  else stopifnot(length(zz) == n * choose(p,2))
  stepwindow <- 10
  if (is.null(aa)) aa <- list(b0=0, th=matrix(0,p,p), bp=rep(0,p), bn=rep(0,p))
  out <- .C("ggdescent_logistic_R",
            xnum,
            as.integer(n),
            as.integer(p),
            zz,
            as.integer(diagonal),
            as.double(y), # convert from integer to double
            as.double(lam.l1),
            as.double(lam.l2),
            as.double(rho),
            as.double(V),
            as.integer(maxiter),
            as.double(aa$b0),
            as.double(aa$th),
            aa$bp,
            aa$bn,
            as.double(stepsize),
            as.integer(stepwindow),
            as.double(backtrack),
            as.double(tol),
            as.integer(trace),
            b0=as.double(0),
            th=rep(0, p*p),
            bp=rep(0, p),
            bn=rep(0, p),
            DUP=FALSE, PACKAGE="hierNet")
  list(b0=out$b0, bp=out$bp, bn=out$bn, th=matrix(out$th, p, p))
}


ADMM4.Lagrangian <- function(aa, xnum, zz, y, lam.l1, lam.l2, diagonal, rho) {
  # aa: list with (th, bp, bn, tt, u)
  # zz is a vector not a matrix
  if (aa$diagonal == FALSE)
    if (any(abs(diag(aa$th)) > 1e-8)) {
      cat("Diagonal of Theta is nonzero.", fill=TRUE)
      return(Inf)
    }
  if (max(aa$tt-t(aa$tt)) > 1e-8) {
    cat("Theta is not symmetric.", fill=TRUE)
    return(Inf)
  }
  if (any(rowSums(abs(aa$th)) > aa$bp + aa$bn + 1e-5)) {
    cat("hierarchy violated.", fill=TRUE)
    return(Inf)
  }
  if (any(aa$bp < -1e-5)||any(aa$bn < -1e-5)) {
    cat("Non-negative of bp or bn violated.", fill=TRUE)
    return(Inf)
  }
  if (diagonal == FALSE)
    if (any(abs(diag(aa$th)) > 1e-5)) {
      cat("Zero diagonal violated.", fill=TRUE)
      return(Inf)
    }

  V <- aa$u - rho * aa$tt

  r <- y - Compute.yhat.c(xnum, zz, aa)
  admm <- sum(aa$u*(aa$th-aa$tt)) + (rho/2) * sum((aa$th-aa$tt)^2)
  #admm <- sum(V*aa$th) + (rho/2) * sum(aa$th^2) + (rho/2)*sum(aa$tt^2) - sum(aa$u*aa$tt)
  pen <- lam.l1 * (sum(aa$bp + aa$bn) + sum(abs(aa$th))/2)
  pen <- pen + lam.l2 * (sum(aa$bp^2) + sum(aa$bn^2) + sum(aa$th^2))
  sum(r^2)/2 + pen + admm
}


predict.hierNet.logistic <- function(object, newx, newzz=NULL, ...) {
  predict.hierNet(object, newx, newzz, ...)
}

critf.logistic <- function(x, y, lam.l1, lam.l2, b0, bp, bn, th) {
  yhat <- b0 + x %*% (bp - bn) + 0.5 * diag(x %*% th %*% t(x))
  p <- 1 / (1 + exp(-yhat))
  val <- -sum(y * log(p) + (1 - y) * log(1 - p)) 
  val <- val + lam.l1 * sum(bp + bn) + lam.l1 * sum(abs(th))/2 + lam.l1 * sum(abs(diag(th)))/2
  val <- val + lam.l2 * (sum(bp^2) + sum(bn^2) + sum(th^2))
  return(val)
}

twonorm <- function(x) {sqrt(sum(x * x))}

hierNet.logistic.path <- function (x, y, lamlist=NULL, delta=1e-8, minlam=NULL, maxlam=NULL, flmin=.01, nlam=20, 
                                   diagonal=TRUE, strong=FALSE, aa=NULL, 
                                   zz=NULL, stand.main=TRUE, stand.int=FALSE,
                                   rho=nrow(x), niter=100, sym.eps=1e-3,# ADMM params
                                   step=1, maxiter=2000, backtrack=0.2, tol=1e-5, # GG params
                                   trace=0) {
  this.call=match.call()
  stopifnot(y %in% c(0, 1))
  x <- scale(x, center=TRUE, scale=stand.main)
  mx <- attr(x, "scaled:center")
  sx <- attr(x, "scaled:scale")
  
  if (is.null(maxlam)) {
    if (!is.null(minlam)) stop("Cannot have maxlam=NULL if minlam is non-null.")
    maxlam <- max(abs(t(x) %*% y))
    minlam <- maxlam * flmin
  }
  if (is.null(minlam)) minlam <- maxlam * flmin
  if (is.null(lamlist))
    lamlist <- exp(seq(log(maxlam), log(minlam), length=nlam))
  nlam <- length(lamlist)
  
  if (is.null(zz)) 
    zz <- compute.interactions.c(x, diagonal=diagonal)
  else stopifnot(is.matrix(zz))
  zz <- scale(zz, center=TRUE, scale=stand.int)
  mzz <- attr(zz, "scaled:center")
  szz <- attr(zz, "scaled:scale")
  zz <- as.numeric(zz)
  p <- ncol(x)
  cp2 <- choose(p, 2)
  b0 <- rep(NA, nlam)
  bp <- bn <- matrix(NA, nrow=p, ncol=nlam)
  th <- array(NA, c(p, p, nlam))
  obj <- rep(NA, nlam)
  aa <- NULL
  for (i in seq(nlam)) {
    cat(c("i,lam=", i, round(lamlist[i],2)), fill=TRUE)
    aa <- hierNet.logistic(x, y, lam=lamlist[i], delta=delta, diagonal=diagonal, strong=strong, 
                           aa=aa, zz=zz, stand.main=FALSE, stand.int=FALSE,
                           rho=rho, niter=niter, sym.eps=sym.eps,
                           step=step, maxiter=maxiter, backtrack=backtrack, tol=tol, 
                           trace=trace)
    b0[i] <- aa$b0
    bp[, i] <- aa$bp
    bn[, i] <- aa$bn
    th[, , i] <- aa$th
    obj[i] <- aa$obj
  }
  dimnames(bp) <- dimnames(bn) <- list(as.character(1:p), NULL)
  dimnames(th) <- list(as.character(1:p), as.character(1:p), NULL)
  out <- list(b0=b0, bp=bp, bn=bn, th=th, obj=obj, lamlist=lamlist, delta=delta,
              mx=mx, my=0, sx=sx, mzz=mzz, szz=szz,
              type="logistic", diagonal=diagonal, strong=strong,
              step=step, maxiter=maxiter, backtrack=backtrack, tol=tol,
              call=this.call)   
  if (strong) {
    # ADMM parameters:
    out$rho <- aa$rho
    out$niter <- niter
    out$sym.eps <- sym.eps
  }
  class(out) <- "hierNet.path"
  out
}

balanced.folds <- function(y, nfolds=min(min(table(y)), 10)) {
  totals <- table(y)
  fmax <- max(totals)
  nfolds <- min(nfolds, fmax)
  # makes no sense to have more folds than the max class size
  folds <- as.list(seq(nfolds))
  yids <- split(seq(y), y)
  # nice we to get the ids in a list, split by class
  ###Make a big matrix, with enough rows to get in all the folds per class
  bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
  for(i in seq(totals)) {
    bigmat[seq(totals[i]), i] <- sample(yids[[i]])
  }
  smallmat <- matrix(bigmat, nrow = nfolds)       # reshape the matrix
  ### Now do a clever sort to mix up the NAs
  smallmat <- permute.rows(t(smallmat))   ### Now a clever unlisting
  # the "clever" unlist doesn't work when there are no N
  #       apply(smallmat, 2, function(x)
  #        x[!is.na(x)])
  res <-vector("list", nfolds)
  for(j in 1:nfolds) {
    jj <- !is.na(smallmat[, j])
    res[[j]] <- smallmat[jj, j]
  }
  return(res)
}

permute.rows <-function(x) {
  dd <- dim(x)
  n <- dd[1]
  p <- dd[2]
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}

hierNet.cv <- function(fit, x, y, nfolds=10, folds=NULL, trace=0) {
  this.call <- match.call()
  stopifnot(class(fit) == "hierNet.path")
  if(fit$type=="gaussian"){errfun=function(y,yhat){(y-yhat)^2}} 
  if(fit$type=="logistic"){errfun=function(y,yhat){1*(y!=yhat)}} 
  n <- length(y)
  if(is.null(folds)) {
    folds <- split(sample(1:n), rep(1:nfolds, length = n))
  }
  else {
    stopifnot(class(folds)=="list")
    nfolds <- length(folds)
  } 
  lamlist=fit$lamlist

  # get whether fit was standardized based on fit$sx and fit$szz...
  if (is.null(fit$mx)) stop("hierNet object was not centered.  hierNet.cv has not been written for this (unusual) case.")
  stand.main <- !is.null(fit$sx)
  stand.int <- !is.null(fit$szz)
  
  n.lamlist <- length(lamlist)        ### Set up the data structures
  size <- double(n.lamlist)
  err2=matrix(NA,nrow=nfolds,ncol=length(lamlist))
  for(ii in 1:nfolds) {
    cat("Fold", ii, ":")
    if(fit$type=="gaussian"){
      a <- hierNet.path(x[-folds[[ii]],],y=y[-folds[[ii]]], 
                        lamlist=lamlist, delta=fit$delta, diagonal=fit$diagonal, strong=fit$strong, trace=trace,
                        stand.main=stand.main, stand.int=stand.int,
                        rho=fit$rho, niter=fit$niter, sym.eps=fit$sym.eps, # ADMM parameters (which will be NULL if strong=F)
                        step=fit$step, maxiter=fit$maxiter, backtrack=fit$backtrack, tol=fit$tol) # GG descent params
      
      yhatt=predict.hierNet(a,newx=x[folds[[ii]],])
    }
    if(fit$type=="logistic"){
      a <- hierNet.logistic.path(x[-folds[[ii]],],y=y[-folds[[ii]]], 
                                 lamlist=lamlist, delta=fit$delta, diagonal=fit$diagonal, strong=fit$strong,
                                 trace=trace, stand.main=stand.main, stand.int=stand.int,
                                 rho=fit$rho, niter=fit$niter, sym.eps=fit$sym.eps, # ADMM parameters (which will be NULL if strong=F)
                                 step=fit$step, maxiter=fit$maxiter, backtrack=fit$backtrack, tol=fit$tol) # GG descent params                                 
      yhatt=predict.hierNet.logistic(a,newx=x[folds[[ii]],])$yhat
    }
    
    temp=matrix(y[folds[[ii]]],nrow=length(folds[[ii]]),ncol=n.lamlist)
    err2[ii,]=colMeans(errfun(yhatt,temp))
    cat("\n")
  }
  errm=colMeans(err2)
  errse=sqrt(apply(err2,2,var)/nfolds)
  o=which.min(errm)
  lamhat=lamlist[o]
  oo=errm<= errm[o]+errse[o]
  lamhat.1se=lamlist[oo & lamlist>=lamhat][1]
  
  nonzero=colSums(fit$bp-fit$bn!=0) + apply(fit$th!=0, 3, function(a) sum(diag(a)) + sum((a+t(a)!=0)[upper.tri(a)]))
  obj <- list(lamlist=lamlist, cv.err=errm,cv.se=errse,lamhat=lamhat, lamhat.1se=lamhat.1se,
             nonzero=nonzero, folds=folds,
             call = this.call)
  class(obj) <- "hierNet.cv"
  obj
}

plot.hierNet.cv <- function(x, ...) {
  par(mar = c(5, 5, 5, 1))
  yrang=range(c(x$cv.err-x$cv.se,x$cv.err+x$cv.se))
  plot(log(x$lamlist), x$cv.err, xlab="log(lambda)",
       ylab = "Cross-validation Error", type="n",ylim=yrang)
  axis(3, at = log(x$lamlist), labels = paste(x$nonzero), srt = 90, adj = 0)
  mtext("Number of features", 3, 4, cex = 1.2)
  axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8))
  error.bars(log(x$lamlist), x$cv.err - x$cv.se, x$cv.err + x$cv.se, width = 0.01, col = "darkgrey")
  points(log(x$lamlist), x$cv.err, col=2, pch=19)
  abline(v=log(x$lamhat), lty=3)
  abline(v=log(x$lamhat.1se), lty=3)
  invisible()
}

error.bars <-function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}

hierNet.varimp <- function(fit,x,y, ...) {
  # NOTE: uses 0.5 cutoff for logistic case
  lam=fit$lam
  if(fit$type=="gaussian"){errfun=function(y,yhat){(y-yhat)^2}}
  if(fit$type=="logistic"){
    errfun=function(y,yhat){
      term1=y*log(yhat);term1[yhat==0]=0
      term2=(1-y)*log(1-yhat);term2[yhat==1]=0
      val=-sum(term1+term2)
      return(val)
    }}
  yhat=predict(fit,x)
  rss=sum(errfun(y,yhat))
  varsum=fit$bp-fit$bn+rowSums(abs(fit$th))
  oo=which(abs(varsum)>1e-6)
  imp=rss2=rep(NA,ncol(x))
  for(j in oo){
    cat(j)
    fit0=fit;fit0$bp=fit$bp[-j];fit0$bn=fit$bn[-j];fit0$th=fit$th[-j,-j]
    if(fit$type=="gaussian"){ fit2=hierNet(x[,-j],y,lam,delta=fit$delta,diagonal=fit$diagonal,aa=fit0)}
    if(fit$type=="logistic"){ fit2=hierNet.logistic(x[,-j],y,lam,delta=fit$delta,diagonal=fit$diagonal,aa=fit0)}
    yhat2=predict(fit2,x[,-j])
    rss2[j]=sum(errfun(y,yhat2))
    imp[j]=(rss2[j]-rss)/rss2[j]
  }
  imp[-oo]=0
  res=cbind(1:ncol(x),round(imp,3))
  ooo=order(-imp)
  dimnames(res)=list(NULL,c("Predictor","Importance"))
  cat("",fill=T)
  return(res[ooo,])
}
