#' @include rng.R
NULL
#' Invert the matrix lam*t(lam)+diag(u) via Woodbury identity
#' @param lam p by k matrix
#' @param u p dimensional vector
#' @return The matrix inverse
#' @export
woodbury = function(lam, u) {
  k = ncol(lam)
  if(is.null(k)) { 
    k=1 
    lam=matrix(lam, nrow=length(u)) 
  }
  if(length(u)==1) {
    out = 1/(lam%*%t(lam)+1/u)
  } else {
    Ui = diag(1/u)
    M = diag(1, k)+t(lam) %*% Ui %*% lam
    out = Ui - Ui %*% lam %*% solve(M, t(lam)) %*% Ui
  }
  return(out)
}

#' Compute samples of the correlation matrix
#' @param model \code{bfa} model object
#' @return A p x p x (no. of mcmc samples) array containing samples of the correlation matrix
#' @export
corr_samp = function(model) {
  pl=model$post.loadings
  ns = dim(pl)[3]
  out = array(NA, dim=c(model$P, model$P, ns))
  lam2 = apply(model$post.loadings, c(1,3), function(x) sum(x^2))
  U = 1/(t(model$post.sigma2)+lam2)
  for (i in 1:ns) {
    pl[, , i]  = pl[, , i]*sqrt(U[,i]) #/sqrt(1 + rowSums(pl[, , i]^2))
    out[, , i] = pl[, , i]%*%t(pl[, , i])+diag(U[,i])
  }
  return(out)
}

#' Compute samples of the covariance matrix
#' @param model \code{bfa} model object
#' @return A p x p x (no. of mcmc samples) array containing samples of the covariance matrix
#' @export
cov_samp = function(model) {
  pl=model$post.loadings
  ns = dim(pl)[3]
  out = array(NA, dim=c(model$P, model$P, ns))
  for (i in 1:ns) {
    out[, , i] = pl[, , i]%*%t(pl[, , i])+diag(model$post.sigma2[i,])
  }
  return(out)
}

#' Posterior predictive and univariate conditional posterior predictive distributions
#' 
#' Posterior predictive and univariate conditional posterior predictive distributions, currently
#' implemented only for Gaussian and copula models. If resp.var is not NA, returns an estimate of the conditional
#' cdf at every observed data point for each MCMC iterate. If resp.var is NA, returns draws from the
#' joint posterior predictive.
#' @param object \code{bfa} model object
#' @param resp.var Either a character vector (length 1) with name of the response variable for 
#' conditional, or NA for draws from the joint posterior predictive. 
#' @param cond.vars Conditioning variables; either a list like list(X1=val1, X2=val2) with
#' X1, X2 variables in the original data frame, or a P length vector with either the conditioning
#' value or NA (for marginalized variables). Ignored if resp.var is NA
#' @param numeric.as.factor Treat numeric variables as ordinal when conditioning
#' @param ... Ignored
#' @return A matrix where each row is either a sample of the conditional posterior predictive 
#' cdf at each datapoint, or a single sample from the joint posterior predictive.
#' @method predict bfa
#' @export
predict.bfa = function(object, resp.var=NA, cond.vars=NA, numeric.as.factor=TRUE, ...) {
  mtype = attr(object, "type")
  if(mtype!="copula") stop("Prediction only available for copula models (see also ?reg_samp)")
  
  n.mcmc=dim(object$post.loadings)[3]
  
  if(!is.na(resp.var)) {
    y.idx = which(resp.var %in% colnames(object$original.data))
    
    if(is.list(cond.vars)) {
      x.var = names(cond.vars)
      x.idx = which(colnames(object$original.data) %in% x.var)
    } else {
      x.idx = which(!is.na(cond.vars))
    }
    
    lo = hi = rep(NA, length(x.var))
    ji=1
    for(j in x.idx) {
      if(is.factor(object$original.data[,j]) || numeric.as.factor) {
        cv = as.numeric(factor(cond.vars[[ji]], 
                        levels=sort(unique(object$original.data[,j]))))
        f = ecdf(object$original.data[,j])
        lo[ji] = qnorm(f(cv-1))
        hi[ji] = qnorm(f(cv))
      } else {
        f = ecdf(object$original.data[,j])
        u = f(cond.vars[[ji]])
        if(u==1) u = (object$N-1)/object$N
        lo[ji] = hi[ji] = qnorm(u)
      }
      ji=ji+1
    }
    y.uq = sort(unique(object$original.data[,y.idx]))
    y.samp   = matrix(NA, nrow=n.mcmc, ncol=length(y.uq))
    y.cutpts = qnorm(ecdf(object$original.data[,y.idx])(y.uq))

    for(i in 1:n.mcmc) {
      A = matrix(object$post.loadings[,,i], nrow=object$P)
      u = 1/(1+rowSums(A^2))
      A = A*sqrt(u)
      k = ncol(A)
      eta = rnorm(k)
      xs = matrix(rnorm(length(x.idx)), ncol=1)
      Ax = matrix(A[x.idx,], nrow=length(x.idx))
      for(h in 1:length(x.idx)) {
        j = x.idx[h]
        xs[h] = rtnorm(lo[h], hi[h], crossprod(A[j,], eta), sqrt(u[j]))
      }
      miv = woodbury(Ax, u[x.idx])
      rc = tcrossprod(A[y.idx[1],], Ax)%*%miv
      mu.y = rc%*%xs
      var.y = 1-A[y.idx[1],] %*% t(Ax) %*% miv %*% Ax %*% matrix(A[y.idx[1],],nrow=k)
      y.samp[i,] = pnorm(y.cutpts, mu.y, sqrt(var.y))
    }
  } else {
    y.samp = matrix(NA, nrow=n.mcmc, ncol=object$P)
    for(i in 1:n.mcmc) {
      A = matrix(object$post.loadings[,,i], nrow=object$P)
      ru = sqrt(1+rowSums(A^2))
      A = A/ru
      k = ncol(A)
      eta = rnorm(k)
      z = A%*%eta+rnorm(object$P, 0, 1/ru)
      y.samp[i,] = z
    }
    y.samp = data.frame(pnorm(y.samp))
    colnames(y.samp) = object$varlabel
    for(j in 1:object$P) {
      y.samp[,j] = quantile(object$ranks[j,]+1, y.samp[,j], type=1)
      y.samp[,j] = object$original.data[match(y.samp[,j], object$ranks[j,]+1), j]
    }
  }
  return(y.samp)
}
