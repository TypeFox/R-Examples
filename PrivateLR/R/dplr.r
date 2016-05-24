################################################################################
#
# File:         linclass.r
# RCS:          $Header: $
# Description:  differentially private logistic regression
# Author:       Staal Vinterbo
# Created:      Thu May 16 11:41:22 2013
# Modified:     Mon Oct 27 18:48:00 2014 (Staal Vinterbo) staal@mats
# Language:     ESS[S]
# Package:      N/A
# Status:       Experimental
#
# (c) Copyright 2013-2014, Staal Vinterbo, all rights reserved.
#
# dplr.r is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# dplr.r is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with linclass.r; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
################################################################################
#
# Revisions:
#
# Fri Oct 24 18:31:11 2014 (Staal Vinterbo) staal@mats
#  Fixed error in noise variance, and implemented full algorithm such that 
#  too small regularizer value does no longer fall back on output perturbation.
################################################################################

# L2 norm
enorm = function(x) sqrt(sum(x^2))

# sample X ~ Erlang(d, l)
rerlang = function(d, lam) {
  if (lam == 0) return (0)
  sum(rexp(d, lam))
}

# sample r in R^d proportionally to exp(alpha * enorm(r))
sampleR = function(d, alpha) {
  dir = rnorm(d)
  dir = dir/enorm(dir)
  rerlang(d, alpha) * dir
}


# average logistic loss
loss.logistic = function(w, x, y) {
  lc = x %*% w
  mean(log1p(exp(-y * lc)))
}

# d(average logistic loss)/dw
dw.loss.logistic = function(w, x, y) {
  lc = x %*% w
  den = as.vector((exp(y*lc) + 1))
  dd = -(y * x) / den
  apply(dd, 2, mean)
}

## # hinge loss
## loss.hinge = function(w, x, y) {
##   lc = x %*% w  
##   mean(pmax(0, 1 - y * lc))
## }

# privacy penalty incurred by regularizer lambda
penalty = function(lambda, n, cc = 0.25) 2 * log(1 + cc/(n * lambda))

# suggest regularizer for given alpha if we want penalty = alpha - alpha.use
lambda.eff = function(alpha, alpha.use, n, cc=0.25)
  cc/n*(1/(exp(0.5 * (alpha - alpha.use)) - 1))

#regularized logistic regression estimation function
# y is {-1, 1} outcome
# x is the design matrix
# alpha is the differential privacy (often called epsilon)
# op if TRUE chooses objective perturbation
# cc is the regularization privacy loss constant
#   for average logistic loss, cc = 0.25
logreg.w = function(y, x, lambda, alpha, op = TRUE, cc = 0.25, verbose = 0, method = 'BFGS') {

  if (lambda <= 0)
    stop('Cannot use non-positive regularizer.')

  status = 'ok'

  n = nrow(x)
  d = ncol(x) # x is the model matrix, it already contains the intercept
  
  if (length(unique(y)) < 2) {
    warning('dplr: unique outcome!', immediate. = TRUE)
    par = c(0.5, rep(0, d - 1))
    names(par) = colnames(x)
    res = list(converge=1,
      par=par, 
      pred = function(x) 0.5, status='unique.outcome')
    class(res) = 'dplr'
    res$eps = alpha
    res$eps.eff = alpha
    res$penalty = penalty
    res$lambda = lambda
    res$n = n
    res$d = d
    res$coefficients = res$par
    res$call = match.call()
    res$CIndex = 0.5    
    return(res)
  }

  # function to do output perturbation
  do.output = function(alpha.use, lambda) {
      f = function(w) (0.5 * lambda * sum(w^2)) + loss.logistic(w, x, y)
      df = function(w) lambda * w + dw.loss.logistic(w, x, y)
      res = optim(rep(0,d), f, df, method=method)
      if (alpha.use > 0) # alpha.use = 0 => non-private 
          res$par = res$par + 2/(alpha.use*lambda*n) * sampleR(d, 1)
      res
  }

  # function to do objective perturbation
  do.perturb = function(alpha.use, lambda) {
      R = sampleR(d, 1) * 2/(alpha.use * n)
      f = function(w)
          (0.5 * lambda * sum(w^2)) + loss.logistic(w, x, y) + (R %*% w)
      df = function(w) lambda * w + dw.loss.logistic(w, x, y) + R
      optim(rep(0,d), f, df, method=method)
  }

  if(op) {
      penalty = 2 * log(1 + cc/(n * lambda))
      alpha.use = alpha - penalty
      if (verbose) {
          cat('dplr: privacy penalty due to regularizer: ', penalty, '\n')
          cat('dplr: alpha - penalty: ', alpha.use, '\n')
      }
      if (alpha.use > 0) {
          res = do.perturb(alpha.use, lambda)
      } else {
          warning('dplr: regularizer too small, using alternative.')
          lambda.adjust = cc/(n * (exp(alpha/4) - 1) )
          res = do.perturb(alpha/2, lambda.adjust)
          alpha.use = alpha/2
          lambda = lambda.adjust
          status = 'adjusted lambda'
      }
  }

  if (!op) { # output perturbation
      alpha.use = alpha
      res = do.output(alpha.use, lambda)
  }

  if (res$convergence > 0) { # positive means no convergence
    status = paste(status, 'nonconvergence')
  }

  names(res$par) = colnames(x)
  class(res) = 'dplr'

  pred.p = 1/(1 + exp(-(x %*% res$par)))
  
  CI = function(y, p, n1 = sum(y == 1))
    c(CI=(mean(rank(p)[y == 1]) - (n1+1)/2)/(length(y)-n1))

  res$CIndex = CI(y, pred.p)
  res$eps = alpha
  res$eps.eff = alpha.use
  res$penalty = penalty
  res$lambda = lambda
  res$n = n
  res$d = d
  res$coefficients = res$par
  res$status = status
  res$call = match.call()
  
  return(res)
}


######### random projections

# compute projection matrix P
randproj = function(p, q) {
  rval = c(-sqrt(3), 0, sqrt(3))
  rprob = c(1/6, 2/3, 1/6)
  rr = sample(rval, p * q, replace=TRUE, prob = rprob)
  m = matrix(rr,
    ncol = q)
}

# do the projection with P on data X
rproject = function(X, P) {
  q = ncol(P)
  1/sqrt(q) * X %*% P
}

# suggest a dimension to project onto
# JL-lemma: preserve distances with factor 1+/-gamma if d = O(gamma^-2 log(n))
rp.suggest = function(n, gamma=0.5, tuning=0.5)
  round(tuning * gamma^(-2) * log(n))


######### dplr interface


# scale covariates to lie within the unit L2 ball
scaled = function(fml, data) {
  mf = model.frame(fml, data)
  mm = model.matrix(terms(mf), mf)
  X = mm[,-1]
  X = scale(X)
  sf = max(apply(X, 1, enorm))
  X = X/sf
  Y = data.frame(X, model.response(mf))
  vnames = attributes(mf)$names
  colnames(Y)[ncol(Y)] = vnames[attributes(terms(mf))$response]
  list(scale = 1/sf, data=Y)
}
    
# recode responses into {fail, succ}
recode.y = function(y, fail = -1, succ = 1) {
  failure = NULL
  if (is.factor(y)){
    levs = levels(y)
    failure = y == levs[1]
  }
  else if (is.logical(y)) {
    failure = !y
  }
  else failure = (y == sort(unique(y))[1])

  y = rep(0, length(y))
  y[failure] = fail
  y[!failure] = succ
  
  y
}

# fitting function used by dplr.formula
dplr.fit = logreg.w

dplr = function(object, ...) UseMethod("dplr")

# main interface that all others are using
dplr.formula = function(object, data, lambda=NA, eps=1, verbose=0, rp.dim = 0, threshold='fixed', do.scale=FALSE,...) {

  scalef = 1

  if(is.na(lambda)) {
    if (eps > 0)
      lambda = lambda.eff(eps, eps - eps/10, nrow(data))
    else
      lambda = 0.001
    if (verbose) cat('dplr: set lambda to: ', lambda, '\n')
  }
  
  mf = model.frame(object, data)
  mm = model.matrix(terms(mf), mf)

  if (do.scale) {
    X = mm[,-1]
    X = scale(X)
    scalef = 1/max(apply(X, 1, enorm))
    mm[,-1] = X*scalef
  }
  
  ndim = ncol(mm[,-1])
  did.rp = FALSE
  if(rp.dim < 0)
    rp.dim = rp.suggest(nrow(data))
  if (rp.dim > 0 && rp.dim < ndim) {
    if(verbose) cat('dplr: randomly projecting onto ', rp.dim, ' dimensions.\n')
    randp = randproj(ndim, rp.dim)
    X = rproject(mm[,-1], randp)
    mm = cbind(mm[, 1], X)
    colnames(mm) = c("(Intercept)",
                paste('dim', 1:ncol(X), sep='.'))
    did.rp = TRUE
  } else {
    if (rp.dim > 0){
      warning('dplr: rp.dim larger than data dimension, not projecting.')
      if (verbose)
        cat('dplr: rp.dim larger than data dimension, not projecting.\n')
    }
  }

  resp = model.response(mf)

  if (is.factor(resp))
    resp.vals = as.factor(levels(resp))
  else if (is.logical(resp))
    resp.vals = c(FALSE, TRUE)
  else resp.vals = sort(unique(resp))
  
  resp = recode.y(resp)
  
  res = dplr.fit(resp, mm, lambda=lambda, alpha=eps, verbose=verbose,...)

  res$call = match.call()
  res$scaled = do.scale
  res$scalef = scalef
  res$terms = terms(mf)
  res$formula = formula(mf)
  res$pred = function(newdata, type = 'probabilities') # nicety 
    predict(res, newdata, type=type)
  res$did.rp = did.rp
  if (did.rp) {
    res$rp.p = randp
    res$rp.dim = rp.dim
  }

  if (threshold == 'youden' || threshold == 'topleft') {
      if (verbose)
          cat('dplr: learning non-private threshold with method', threshold, '\n')
      p = predict(res, data)
      res$p.tr = p.tr(p, resp)
  } else {
      res$p.tr = 0.5
  }
  res$resps = resp.vals
  
  return(res)
}

predict.dplr = function(object, data, type = 'probabilities',  ...) {
  terms = delete.response(terms(object))
  data = data.frame(data)
  mm = model.matrix(terms, data)

  if (object$scaled) { # rescale input
    X = scale(mm[,-1])
    scalef = 1/max(apply(X, 1, enorm))
    mm[,-1] = X*scalef    
  }
  
  if (object$did.rp) {
    X = rproject(mm[,-1], object$rp.p)
    mm = cbind(mm[, 1], X)
  }

  p = as.vector(1/(1 + exp(-(mm %*% object$par))))
  if (type[1] == 'probabilities') 
    return (p)
  i = (p > object$p.tr) + 1
  return (object$resps[i])
  
}


summary.dplr = function(object,...) {
    l = list(formula=list(
                 formula = paste(format(object$formula), collapse='\n')),
        common=list(epsilon = paste(' Privacy level:', object$eps),
        lambda = paste(' Regularization parameter:', object$lambda),
        convergence=paste(' Convergence of optimization :', object$convergence == 0),
        c.index=paste(' C-index (AUC):', object$CIndex),
        status=paste(' Status:', object$status)            
        ))

    if (object$did.rp)
        l[[2]] = c(l[[2]], 
             rp.dim = paste(' Number of dimensions in random projection:',
                 object$rp.dim))
                             
    l = c(l, list(list(coefficients=object$coefficients)))
    class(l) = 'summary.dplr'
    l
}

print.summary.dplr = function(x,...){
    for (l in x) {
        for (ll in l){
            if (class(ll) == 'character')
                cat(ll, '\n')
            else
                print(ll)
        }
        cat('\n')
    }
    invisible(x)
}

print.dplr = function(x,...) {
  print(summary(x))
  invisible(x)
}
# other interfaces:

# y is a vector of outcomes,
# x is a matrix/data frame of parameters.
# assumes that no column in x is called dprl.class.
dplr.numeric = function(object, x, ...) {
  df = data.frame(cbind(dplr.class=object, x))
  return(dplr(dplr.class ~ ., data=df, ...))
}

# true/false outcomes
dplr.logical = function(object, x, ...) dplr.numeric(object, x, ...)

# factor outcomes
dplr.factor = function(object, x, ...) dplr.numeric(object, x, ...)


# data frame with target
dplr.data.frame = function(object, target=ncol(object), ...){
  fmla = as.formula(paste(names(object[target]), '~ .'))
  dplr(fmla, data=object,...)
}

# matrix with target
dplr.matrix = function(object, target=ncol(object), ...)
  dplr(data.frame(object), target=target, ...)
  
  
############### ROC analysis helpers

# produce sensitivity and specificities
sesp = function(p, y) {

  r = rank(p, ties.method='first')
  r1 = sort(r[y == 1])

  p1 = sort(p)[r1]

  n0 = sum(y != 1)
  n1 = sum(y == 1)

  yiv = 0
  yi = 0

  cl = 1
  cli = 0

  r = matrix(ncol=3)
  
  for (i in 1:length(r1)) {

    pbelow = i - 1
    pabove = n1 - i

    nbelow = r1[i] - pbelow
    nabove = n0 - nbelow

    #fn = zbelow
    tp = pabove + 1 # count current
    #fp = oabove
    tn = nbelow

    se = tp/n1
    sp = tn/n0

    r = rbind(r, c(p1[i], se, sp))
  }
  r
}

# compute best threshold
p.tr = function(p, y, method='youden'){
  x = sesp(p, y)
  if (method == 'youden')
    i = which.max(apply(x[,-1], 1, sum))
  else # point closest to (0,1)
    i = which.min(apply(x[,-1], 1, function(z)
      (1 - z[1])^2 + (1 - z[2])^2))
  x[i, 1]
}





  
    
