##' spinyReg
##'
##' Computethe path of solution of a spinyReg fit.
##'
##' @param X matrix of features. Do NOT include intercept.
##'
##' @param Y matrix of responses.
##'
##' @param alpha numeric scalar; prior value for the alpha parameter
##' (see the model's details). Default is 0.1.
##'
##' @param gamma numeric scalar; prior value for the gamma parameter
##' (see the model's details). Default is 1.
##'
##' @param z numeric vector; prior support of active variable. Default
##' is \code{rep(1,p)}, meaning all variable activated
##'
##' @param intercept logical; indicates if a vector of intercepts
##' should be included in the model. Default is \code{TRUE}.
##'
##' @param normalize logical; indicates if predictor variables should
##' be normalized to have unit L2 norm before fitting.  Default is
##' \code{TRUE}.
##'
##' @param verbose integer; activate verbose mode from '0' (nothing)
##' to '2' (detailed output).
##'
##' should be included in the model. Default is \code{TRUE}.
##'
##' @param recovery logical; indicates if the full path of models
##' should be inspected for model selection. Default is \code{TRUE}.
##'
##' @param maxit  integer; the maximal number of iteration
##' (i.e. number of alternated optimization between each parameter)
##' in the Expectation/Maximization algorithm.
##'
##' @param eps a threshold for convergence. Default is \code{1e-10}.
##'
##' @return an object with class \code{spinyreg}, see the
##' documentation page \code{\linkS4class{spinyreg}} for details.
##'
##' @seealso See also \code{\linkS4class{spinyreg}}.
##' @name spinyReg
##' @rdname spinyReg
##' @keywords models, regression
##'
##' @examples \dontrun{
##' data <- read.table(file="http://statweb.stanford.edu/~tibs/ElemStatLearn/datasets/prostate.data")
##' x <- data[, 1:8]
##' y <- data[, 9]
##' out <- spinyreg(x,y,verbose=2)
##' }
##' @importFrom stats optim lm.fit
##' @export
spinyreg <-function(X,
                    Y,
                    alpha     = 1e-1,
                    gamma     = 1,
                    z         = rep(1,ncol(X)),
                    intercept = TRUE,
                    normalize = TRUE,
                    verbose   = 1,
                    recovery  = TRUE,
                    maxit     = 1000,
                    eps       = 1e-10) {

  X <- as.matrix(X)

  n <- nrow(X); p <- ncol(X)
  x <- X
  y <- Y
  if (is.null(colnames(X))) {colnames(X) <- 1:p}
  names <- colnames(X)

  ## ======================================================
  ## INTERCEPT TREATMENT
  if (intercept) {
    names <- c("intercept",names)
    x.bar <- colMeans(x)
    x     <- scale(x,x.bar,FALSE)
    y.bar <- mean(y)
    y <- y-y.bar
  }

  ## ======================================================
  ## NORMALIZATION

  ## normalizing the data
  if (normalize) {
    normx <- sqrt(drop(colSums(x^2)))
    x    <- scale(x, FALSE, normx)
  } else {
    normx <- rep(1,p)
  }

  xtx <- crossprod(x)
  yty <- crossprod(y)
  xty <- crossprod(x,y)

  loglik <- c()

  cond <- FALSE
  i <- 0
  while(!cond) {
    i <- i+1

    ## E-step
    if(verbose){cat('\n E step (parameter estimation)')}
    ## Call to relevance
    out <- E.step(z, x, xtx, yty, xty, alpha, gamma)
    alpha <- out$alpha
    gamma <- out$gamma
    loglik <- c(loglik,out$loglik)
    Mn <- out$Mn
    Sigma <- out$Sigma
    if(verbose){cat('\t alpha, gamma =', alpha, gamma,"\n")}

    ## M step
    if(verbose){cat('\n M step: support recovery')}
    z <- M.step(z, gamma, Mn, Sigma, xty, xtx)
    if(verbose){cat('\t card(z) =', sum(abs(z) > .Machine$double.eps), "\n")}

    if(i>1) {
      cond <- (i >= maxit) | ((loglik[i]-loglik[i-1])^2/loglik[i]^2 < eps)
    }

    if (is.na(cond)) {
      return(list(q=NA,w=rep(NA,p),coef=rep(NA,p),w=NA,delta=NA,zrelax=NA,alpha=NA,gamma=NA,loglik=NA,margin=NA))
    }
  }

  # Binarisation selon le chemin de sparsite + astuce de reorganisation
  id <- order(z,decreasing=TRUE)
  zstart <- as.numeric(z == 1)
  q0 = sum(zstart)
  if (q0 == 0) {zstart<-rep(0, p); zstart[id[1]]=1; q0 = 1}
  margin <- rep(E.step( zstart, x, xtx, yty, xty, alpha, gamma)$loglik,q0)

  if(verbose){cat("\n Model selection and sparsity path\n")}
  for(q in (q0+1):round(3*p/4)) {
    zstart<-rep(0, p); zstart[id[1:q]]=1
    margin<-c(margin, E.step(zstart, x, xtx, yty, xty, alpha, gamma)$loglik)
    if (recovery) {
      if ((q <= round(p/2)) & (margin[q-1] > margin[q])){
        Etmp = c()
        for(j in 1:round(p/10)){
          ztmp = rep(0, p); ztmp[id[1:(q-1)]]=1
          ztmp[id[(q+j)]]=1
          Etmp = c(Etmp,E.step(ztmp, x, xtx, yty, xty, alpha, gamma)$loglik)
        }
        if (max(Etmp) > margin[q-1]) {
          pos = which.max(Etmp)
          idtmp = id[q]; id[q] = id[q+pos]; id[q+pos] = idtmp
          margin[q] = max(Etmp)
        }
      }
    }
  }
  margin[!is.finite(margin)] = -Inf
  q = which.max(margin)
  zopt <- as.numeric(z == 1)
  zopt[id[1:q]]=1

  ## ======================================================
  ## RESCALING, INTERCEPT TREATMENT
  if (verbose >1) {
    cat("\n Normalize back to the original scale")
  }

  ## coef resestimated with OLS
  coef <- rep(0,p)
  coef[which(zopt != 0)] <- lm.fit(x[, which(zopt != 0)],y)$coefficients

  ## re-scaling
  if (normalize) {
    coef     <- coef/normx
  }
  ## intercept
  if (intercept) {
    coef <- c(y.bar - crossprod(coef,x.bar), coef)
    fitted     <- c(cbind(1,X) %*% coef)
  } else {
    fitted     <- c(X %*% coef)
  }
  names(coef) <- names

  residuals <- as.numeric(Y - fitted)
  r.squared <- 1-sum(residuals^2)/sum(Y^2)

  if (verbose >1) {
    cat("\n DONE!\n")
  }

  ## Then we are done
  return(
      new("spinyreg",
          coefficients = coef,
          alpha        = alpha      ,
          gamma        = gamma      ,
          residuals    = residuals  ,
          fitted       = fitted     ,
          r.squared    = r.squared  ,
          normx        = normx      ,
          intercept    = intercept,
          monitoring   = list(loglik   = c(loglik),
          margin   = c(margin),
          iterates = i)
          )
      )
}

####################################################################################################
E.step <- function(z, x, xtx, yty, xty, alpha, gamma) {

  n <- nrow(x); p <- ncol(x)
  Z <- as.matrix(diag(z))
  cond <- FALSE
  loglik <- c()

  ## Ensure an O(p^2 min(n,p)) complexity
  if (p>n) { # woodburry
    XZ <- sweep(x,2,z,"*",check.margin=FALSE)
    chol.out <- try(chol( diag(n)/gamma + tcrossprod(XZ)/alpha ), TRUE)
    if (!is.character(chol.out)) {
      Sn <- diag(p)/alpha - t(XZ) %*% chol2inv(chol.out) %*% (XZ) / (alpha^2)
    } else {
      Sninv <- gamma * Z %*% xtx %*% Z + alpha * diag(p)
      Sn <- chol2inv(chol(Sninv))
    }
  } else {
    Sninv <- gamma * Z %*% xtx %*% Z + alpha * diag(p)
    Sn <- chol2inv(chol(Sninv))
  }

  Mn <- gamma * Sn %*% (z*xty)
  Sigma <- Sn + tcrossprod(Mn)
  dSn <- det(Sn)
  if(dSn == 0) dSn <- .Machine$double.eps
  if(dSn == Inf) dSn <- .Machine$double.xmax

  Mn.z.xty <- crossprod(Mn,z*xty)
  loglik <- .5*(gamma*(Mn.z.xty-yty) + log(dSn) + p*log(alpha) + n*log(gamma/(2*pi)))

  alpha <- p/sum(diag(Sigma))
  gamma <- as.numeric(n/(yty + crossprod(crossprod(xtx*Sigma, z),z) - 2*Mn.z.xty))

  return(list(alpha=alpha, gamma=gamma, loglik=loglik, Mn=Mn, Sigma=Sigma))
}

################################################################################################
M.step <- function(z, gamma, Mn, Sigma, xty, xtx) {

  fn <- function(x) {return(gamma*(.5*crossprod(crossprod(xtx*Sigma, x),x) - crossprod(x, Mn*xty)))}
  gr <- function(x) {return((gamma*xtx*Sigma) %*% x - gamma*Mn*xty)}

  res <- optim(z, fn=fn, gr=gr, lower=0, upper=1, method="L-BFGS-B")
  return(res$par)
}
