rmvnorm = mvtnorm::rmvnorm

dot = function(x,y) sum(x*y)

## should be faster/more numerically stable than 'solve' for inverting covariance matrices
symmsolve = function(V) chol2inv(chol(V))

rwishart = function(dof, V) {
  ## NB the final parameter usually included in the definition of the Wishart distribution is implicit in the size of V
  ## To sample from the inverse-Wishart distribution, use solve(rwishart(dof, solve(V)))
  A = rmvnorm(dof, sigma=V)
  t(A) %*% A
}

rdirichlet = function(alpha) {
  y = rgamma(length(alpha), alpha)
  y / sum(y)
}

Gibbs.regression = function(x.in, y.in, M, Nsamples, Ngauss=1, dirichlet=FALSE, M.inv=NULL, intercept=TRUE, trace='bs', fix='', start=list(), B.prior.mean=NULL, B.prior.cov=NULL, Sigma.prior.scale=NULL, Sigma.prior.dof=-1, dp.prior.alpha=NULL, dp.prior.beta=NULL, mention.every=NA, save.every=NA, save.to=NA) {
  ## n data points, p covariates and m responses
  ## x.in should be n*p, and y.in n*m
  ## but also allow x.in and/or y.in to be vectors (i.e. p=1 and/or m=1)
  ## To fit only intercepts, pass something of length 0, e.g. double(), as x.in
  ## M should be a length n list of (p+m)^2 covariances.
  ## Alternatively, M.inv can be set to their inverses.
  ## If both M and M.inv have length 0 (NULL), do not update X or Y (with a warning).
  ## Ngauss=0 is a special case where we take the limit of a uniform prior on x.
  ## (This is not advised.) If dirichlet=TRUE, we use a Dirichlet process to
  ## effectively learn Ngauss from the data.
  ## 
  ##
### Sort through the options and set up everything
  ##
  ##
  update.X = !grepl('x', tolower(fix))
  update.Y = !grepl('y', tolower(fix))
  update.B = !grepl('b', tolower(fix))
  update.Sigma = !grepl('s', tolower(fix))
  update.G = !grepl('g', tolower(fix))
  update.pi = !grepl('p', tolower(fix))
  update.mu = !grepl('m', tolower(fix))
  update.Tau = !grepl('t', tolower(fix))
  update.mu0 = !grepl('z', tolower(fix))
  update.U = !grepl('u', tolower(fix))
  update.W = !grepl('w', tolower(fix))
  update.alpha = !grepl('a', tolower(fix))
  if (is.null(dim(y.in))) yy = matrix(y.in, length(y.in), 1) else yy = y.in
  if (is.null(dim(x.in))) {
    if (length(x.in) == 0) {
      if (!intercept) warning('Gibbs.regression: intercept forced to TRUE because no covariate measurements were passed.')
      intercept = TRUE
      xx = matrix(nrow=nrow(yy), ncol=0)
    } else {
      xx = matrix(x.in, length(x.in), 1)
    }
  } else {
    xx = x.in
  }
  n = nrow(xx)
  p = ncol(xx)
  m = ncol(yy)
  if (nrow(yy) != nrow(xx)) stop('Gibbs.regression: number of covariate and response measurements do not match.')
  if (length(M) > 0) {
    if (nrow(yy) != length(M)) stop('Gibbs.regression: number of measurement covariance entries does not match number of measurements.')
    for (i in 1:n) if (!is.matrix(M[[i]]) | nrow(M[[i]]) != p+m | ncol(M[[i]]) != p+m) stop(paste('Gibbs.regression: measurement covariance', i, 'is not a matrix or has the wrong size/shape.'))
    M.inv = list()
    for (i in 1:n) M.inv[[i]] = symmsolve(M[[i]])
  } else if (length(M.inv) > 0) {
    if (nrow(yy) != length(M.inv)) stop('Gibbs.regression: number of inverse measurement covariance entries does not match number of measurements.')
    for (i in 1:n) if (!is.matrix(M.inv[[i]]) | nrow(M.inv[[i]]) != p+m | ncol(M.inv[[i]]) != p+m) stop(paste('Gibbs.regression: inverse measurement covariance', i, 'is not a matrix or has the wrong size/shape.'))
  } else {
    ## neither M nor M.inv were passed
    if (update.X | update.Y) warning('Gibbs.regression: fixing X and Y because no measurement covariances information was passed.')
    update.X = FALSE
    update.Y = FALSE
  }
  Ngauss = as.integer(Ngauss)
  if (dirichlet) {
    Ngauss = 1
    ## defaults from R.M. Dorazio / Journal of Statistical Planning and Inference 139 (2009) 3384 -- 3390, Table 1
    if (length(dp.prior.alpha) == 0) dp.prior.alpha = approx(seq(5,50,5), c(0.541, 0.525, 0.512, 0.501, 0.490, 0.486, 0.480, 0.475, 0.470, 0.467), n, rule=2)$y else if (length(dp.prior.alpha) != 1 | dp.prior.alpha <= 0) stop('Gibbs.regression: invalid value of dp.prior.alpha argument')
    if (length(dp.prior.beta) == 0) dp.prior.beta = approx(seq(5,50,5), c(0.096, 0.046, 0.029, 0.021, 0.015, 0.013, 0.010, 0.009, 0.008, 0.007), n, rule=2)$y else if (length(dp.prior.beta) != 1 | dp.prior.beta <= 0) stop('Gibbs.regression: invalid value of dp.prior.beta argument')
    dp.concentration = rgamma(1, dp.prior.alpha, rate=dp.prior.beta)
    if ('alpha' %in% names(start)) {
      dp.concentration = start$alpha
      if (length(dp.concentration) != 1 | dp.concentration <= 0) stop('Gibbs.regression: invalid value of start$alpha')
    }
  }
  if (Ngauss < 0) stop('Gibbs.regression: invalid value of Ngauss argument')
  pin = 0 # extra covariate corresponding to an intercept term
  if (intercept) pin = 1
  if (p > 0) prange = pin+1:p else prange = 0
  if ('X' %in% names(start)) {
    X = start$X
    if (nrow(X) != n | ncol(X) != pin+p) stop('Gibbs.regression: dimensions of start$X are wrong.')
  } else {
    if (intercept) {
      X = cbind(rep(1, n), xx) # true x values - start at measurements; n*(p+pin)
    } else {
      X = xx # true x values - start at measurements; n*(p+pin)
    }
  }
  if ('Y' %in% names(start)) {
    Y = start$Y
    if (nrow(Y) != n | ncol(Y) != m) stop('Gibbs.regression: dimensions of start$Y are wrong.')
  } else {
    Y = yy # true y values - start at measurements; n*m
  }
  ## need a starting guess for either coefficients or intrinsic scatter
  if ('B' %in% names(start)) {
    B = start$B
    if (nrow(B) != pin+p | ncol(B) != m) stop('Gibbs.regression: dimensions of start$B are wrong.')
  } else {
    if (intercept) {
      B = rbind(colMeans(Y), matrix(0, p, m)) # craptastic guess for coefficients; (p+pin)*m
    } else {
      B = matrix(0, p, m) # craptastic guess for coefficients; (p+pin)*m
    }
  }
  if ('Sigma' %in% names(start)) {
    Sigma = start$Sigma
    if (nrow(Sigma) != m | ncol(Sigma) != m) stop('Gibbs.regression: dimensions of start$Sigma are wrong.')
    Sigma.inv = symmsolve(Sigma)
  } else {
    if (!update.Sigma) stop('Gibbs.regression: a starting value must be passed if Sigma is fixed.')
  }
  ##
  res = list()
  return.X = grepl('x', tolower(trace))
  return.Y = grepl('y', tolower(trace))
  return.B = grepl('b', tolower(trace))
  return.Sigma = grepl('s', tolower(trace))
  return.G = grepl('g', tolower(trace))
  return.pi = grepl('p', tolower(trace))
  return.mu = grepl('m', tolower(trace))
  return.Tau = grepl('t', tolower(trace))
  return.mu0 = grepl('z', tolower(trace))
  return.U = grepl('u', tolower(trace))
  return.W = grepl('w', tolower(trace))
  return.alpha = grepl('a', tolower(trace))
  ##
  nogauss = TRUE
  if (p > 0) {
    nogauss = FALSE
    if (Ngauss == 0) {
      nogauss = TRUE
      Ngauss = 1
    }
    if (dirichlet) {
      uniqueX = which(!duplicated(X))
      nclusters = length(uniqueX)
      nG = integer(nclusters)
      G = integer(n)
      for (j in 1:nclusters) {
        gg = which(apply(X, 1, function(x) all(x == X[uniqueX[j],])))
        G[gg] = j
        nG[j] = length(gg)
      }
    } else {
      update.alpha = FALSE
      return.alpha = FALSE
      if ('G' %in% names(start)) {
        G = start$G
        if (length(G) != n) stop('Gibbs.regression: dimensions of start$G are wrong.')
      } else {
        G = sample(1:Ngauss, n, replace=TRUE) # index of which gaussian component each data point belongs to
        ## note that this is different from (more compact than) the representation of G defined in Kelly 2007
      }
      nG = array(dim=Ngauss)
      for (k in 1:Ngauss) nG[k] = length(which(G == k))
    }
    ppii = nG / sum(nG)
    if ('pi' %in% names(start)) {
      ppii = start$pi
      if (length(ppii) != Ngauss) stop('Gibbs.regression: dimensions of start$pi are wrong.')
    }
    if ('mu0' %in% names(start)) {
      mu0 = start$mu0
      if (length(mu0) != p) stop('Gibbs.regression: dimensions of start$mu0 are wrong.')
    } else {
      mu0 = colMeans(xx) # length p
    }
    if ('U' %in% names(start)) {
      U = start$U
      if (nrow(U) != p | ncol(U) != p) stop('Gibbs.regression: dimensions of start$U are wrong.')
    } else {
      if (p > 1) {
        U = diag(apply(xx, 2, 'var')) # p^2
      } else {
        U = var(xx)
      }
    }
    U.inv = symmsolve(U)
    if ('mu' %in% names(start)) {
      mu = start$mu
      if (nrow(mu) != Ngauss | ncol(mu) != p) stop('Gibbs.regression: dimensions of start$mu are wrong.')      
    } else {
      mu = rmvnorm(Ngauss, mu0, U) # Ngauss*p
    }
    Tau = list() # length Ngauss of p^2 matrices
    Tau.inv = list()
    if (nogauss) {
      Tau.inv[[1]] = matrix(0, nrow=p, ncol=p)
      update.mu = FALSE
      update.Tau = FALSE
      return.mu = FALSE
      return.Tau = FALSE
    } else {
      if ('Tau' %in% names(start)) {
        if (!all(dim(start$Tau) == c(p,p,Ngauss))) stop('Gibbs.regression: dimensions of start$Tau are wrong.')
        for (k in 1:Ngauss) Tau[[k]] = start$Tau[,,k]
        Tau.inv[[k]] = symmsolve(Tau[[k]])
      } else {
        if (!update.Tau) stop('Gibbs.regression: a starting value must be passed if Tau is fixed and Ngauss > 0.')
      }
    }
    if (Ngauss > 1) {
      if ('W' %in% names(start)) {
        W = start$W
      if (nrow(W) != p | ncol(W) != p) stop('Gibbs.regression: dimensions of start$W are wrong.')              
      } else {
        W = cov(xx) # p^2
      }
      nu.Tau = p
    } else {
      update.G = FALSE
      update.pi = FALSE
      update.mu0 = FALSE
      update.U = FALSE
      update.W = FALSE
      if (!dirichlet) return.G = FALSE
      return.pi = FALSE
      return.mu0 = FALSE
      return.U = FALSE
      return.W = FALSE
      W = matrix(0, nrow=p, ncol=p)
      nu.Tau = 0
      U.inv = matrix(0, nrow=p, ncol=p)
    }
  } else {
    update.X = FALSE
    update.G = FALSE
    update.pi = FALSE
    update.mu = FALSE
    update.Tau = FALSE
    update.mu0 = FALSE
    update.U = FALSE
    update.W = FALSE
    update.alpha = FALSE
    return.X = FALSE
    return.G = FALSE
    return.pi = FALSE
    return.mu = FALSE
    return.Tau = FALSE
    return.mu0 = FALSE
    return.U = FALSE
    return.W = FALSE
    return.alpha = FALSE
  }
  ##
  Bprior = FALSE
  if (length(B.prior.cov) != 0) {
    Bprior = TRUE
    if (!is.matrix(B.prior.cov) | nrow(B.prior.cov) != (pin+p)*m | ncol(B.prior.cov) != (pin+p)*m) stop(paste('Gibbs.regression: B.prior.cov argument must be a', (pin+p)*m, 'x', (pin+p)*m, 'matrix'))
    if (length(B.prior.mean) == 0) {
      B.prior.mean = rep(0, (pin+p)*m)
    } else if (length(B.prior.mean) != (pin+p)*m) stop(paste('Gibbs.regression: B.prior.mean argument must have length', (pin+p)*m))
    B.prior.icov = symmsolve(B.prior.cov)
  }
  if (length(Sigma.prior.scale) == 0) {
    Sigma.prior.scale = matrix(0, m, m)
  } else {
    if (!is.matrix(Sigma.prior.scale) | nrow(Sigma.prior.scale) != m | ncol(Sigma.prior.scale) != m) stop(paste('Gibbs.regression: Sigma.prior.scale argument must be a', m, 'x', m, 'matrix'))
  }
  ##
  if (return.X) res$X = array(dim=c(n, p+pin, Nsamples))
  if (return.Y) res$Y = array(dim=c(n, m, Nsamples))
  if (return.B) res$B = array(dim=c(p+pin, m, Nsamples))
  if (return.Sigma) res$Sigma = array(dim=c(m, m, Nsamples))
  if (return.G) res$G = array(dim=c(n, Nsamples))
  if (return.pi) res$pi = array(dim=c(Ngauss, Nsamples))
  if (return.mu) res$mu = array(dim=c(Ngauss, p, Nsamples))
  if (return.Tau) res$Tau = array(dim=c(p, p, Ngauss, Nsamples))
  if (return.mu0) res$mu0 = array(dim=c(p, Nsamples))
  if (return.U) res$U = array(dim=c(p, p, Nsamples))
  if (return.W) res$W = array(dim=c(p, p, Nsamples))
  if (return.alpha) res$alpha = array(dim=c(Nsamples))
  ##
  ##
### Run the Gibbs sampler Nsamples times
  ##
  ##
  for (Isample in 1:Nsamples) {
    ##
    ## update the covariances of the gaussian components 
    if (!nogauss & update.Tau) {
      for (k in 1:Ngauss) {
        if (dirichlet) { # NB Ngauss=1 autmatically in this case
          ii = uniqueX
        } else {
          ii = which(G == k)
        }
        S = W
        for (i in ii) {
          z = X[i, prange] - mu[k,]
          S = S + z %*% t(z)
        }
        Tau.inv[[k]] = rwishart(length(ii) + nu.Tau, symmsolve(S))
        Tau[[k]] = symmsolve(Tau.inv[[k]])
      }
    }
    ##
    ## update the component proportions
    if (!nogauss & update.pi) {
      ppii = rdirichlet(1 + nG)
    }
    ##
    ## update the intrinsic scatter
    if (update.Sigma) {
      E = Y - X %*% B # n*m
      Sigma.inv = rwishart(n+Sigma.prior.dof, symmsolve(t(E)%*%E + Sigma.prior.scale)) # m^2
      Sigma = symmsolve(Sigma.inv)
    }
    ##
    ## update the coefficients
    if (update.B) {
      Y.tilde = as.matrix(as.vector(Y), nrow=n*m) # nm*1
      B.tilde.cal = matrix(0, (p+pin)*m, 1) # (p+pin)m*1
      S.tilde.cal = matrix(0, (p+pin)*m, (p+pin)*m) # [(p+pin)m]^2
      ## old code; replacement below should be more stable, in principle.
      ## XtXinv = solve(t(X) %*% X) # (p+pin)^2
      ## XtXinvXt = XtXinv %*% t(X)
      ## for (j in 1:m) B.tilde.cal[1:(p+pin)+(j-1)*(p+pin),1] = XtXinvXt %*% Y.tilde[1:n+(j-1)*n,1]
      ## .. and now the real thing:
      qx = qr(X) # QR decomposition of X
      for (j in 1:m) B.tilde.cal[1:(p+pin)+(j-1)*(p+pin),1] = qr.solve(qx, Y.tilde[1:n+(j-1)*n,1]) # each chunk (p+pin)*1 is the solution (B) to X*B=Y
      XtXinv = chol2inv(qx$qr) # fancy replacement for solve(t(X) %*% X); (p+pin)^2
      for (i in 1:m) for (j in 1:m) S.tilde.cal[1:(p+pin)+(i-1)*(p+pin), 1:(p+pin)+(j-1)*(p+pin)] = XtXinv * Sigma[i,j]
      if (Bprior) {
        S.tilde.cal = symmsolve(S.tilde.cal)
        B.tilde.cal = S.tilde.cal %*% B.tilde.cal + B.prior.icov %*% B.prior.mean
        S.tilde.cal = symmsolve(S.tilde.cal + B.prior.icov)
        B.tilde.cal = S.tilde.cal %*% B.tilde.cal
      }
      B.tilde = as.vector(rmvnorm(1, B.tilde.cal, S.tilde.cal))
      B = matrix(B.tilde, p+pin, m)
    }
    ##
    ## update the covariates
    if (update.X) {
      B.Sinv = B[prange,] %*% Sigma.inv # p*m
      if (dirichlet) {
        ## algorithm 2 of Neal 2000 (http://www.jstor.org/stable/1390653)
        B.Sinv.B = B.Sinv %*% t(B[prange,]) # p*p
        Tinv.mu = Tau.inv[[1]] %*% mu[1,]
        if (pin == 0) {
          alpha = rep(0, m)
        } else {
          alpha = B[pin,]
        }
        for (i in 1:n) {
          Fcov.inv = M.inv[[i]][1:p,1:p] + B.Sinv.B
          Fcov = symmsolve(Fcov.inv)
          zi = c(xx[i,], yy[i,]-Y[i,])
          Fpart = (M.inv[[i]] %*% zi)[1:p] + B.Sinv %*% (Y[i,] - alpha)
          Fmean = Fcov %*% Fpart
          inds = which(nG > 0)
          if (nG[G[i]] == 1) inds = inds[-which(inds == G[i])]
          q = rep(0, nclusters)
          for (j in inds) q[j] = dmvnorm(X[uniqueX[j],prange], Fmean, Fcov)
          n.minus.1 = nG
          n.minus.1[G[i]] = pmax(0, n.minus.1[G[i]] - 1)
          q = q * n.minus.1
          qsum = sum(q)
          r = dp.concentration * dmvnorm(mu[1,], Fmean, Fcov+Tau[[1]])
          b = 1/(qsum + r)
          if (runif(1) <= b*qsum) {
            newg = sample(1:nclusters, 1, prob=q)
            X[i,] = X[uniqueX[newg],]
            nG[G[i]] = nG[G[i]] - 1
            nG[newg] = nG[newg] + 1
            G[i] = newg
          } else {
            Fcov = symmsolve( Fcov.inv + Tau.inv[[1]] )
            Fmean = Fcov %*% (Fpart + Tinv.mu)
            X[i,prange] = rmvnorm(1, Fmean, Fcov)
            nclusters = nclusters + 1
            nG[G[i]] = nG[G[i]] - 1
            G[i] = nclusters
            nG[nclusters] = 1
            uniqueX = c(uniqueX, i)
          }
        }
        uniqueX = which(!duplicated(X))
        nclusters = length(uniqueX)
        nG = integer(nclusters)
        if (return.G) for (j in 1:nclusters) {
          gg = which(apply(X, 1, function(x) all(x == X[uniqueX[j],])))
          G[gg] = j
          nG[j] = length(gg)
        }
        for (j in 1:nclusters) {
          Fcov.inv = Tau.inv[[1]]
          Fpart = Tinv.mu
          gg = which(G == j)
          for (i in gg) {
            Fcov.inv = Fcov.inv + M.inv[[i]][1:p,1:p] + B.Sinv.B
            zi = c(xx[i,], yy[i,]-Y[i,])
            Fpart = Fpart + (M.inv[[i]] %*% zi)[1:p] + B.Sinv %*% (Y[i,] - alpha)
          }
          Fcov = symmsolve(Fcov.inv)
          Fmean = Fcov %*% Fpart
          xnew = rmvnorm(1, Fmean, Fcov)
          for (i in gg) X[i,prange] = xnew
        }
      } else {
        B.Sinv.B.j = rep(NA, p)
        for (j in 1:p) B.Sinv.B.j[j] = dot(B.Sinv[j,], B[pin+j,])
        xi.hat = rep(NA, p)
        s2 = rep(NA, p)
        for (i in 1:n) {
          ## Stripped-down pin=1, p=1, m=1, k=0 version for sanity checking
          ## j = 1
          ## rho = M[[i]][1,2]/sqrt(M[[i]][1,1]*M[[i]][2,2])
          ## s22 = 1/( 1/(M[[i]][1,1]*(1-rho^2)) + B[2,1]^2/Sigma[1,1]  )
          ## thing = xx[i,1] + M[[i]][1,2]/M[[i]][2,2] * (Y[i,1] - yy[i,1])
          ## xi.hat2 = ( thing/(M[[i]][1,1]*(1-rho^2)) + B[2,1]*(Y[i,1]-B[1,1])/Sigma[1,1] ) * s22
          ## X[i,2] = rnorm(1, xi.hat2, sqrt(s22))
          ## .. and now the real thing:
          zi = c(xx[i,], yy[i,]) - c(X[i,prange], Y[i,]) # p+m
          mui = mu[G[i],] - X[i,prange] # p
          for (j in 1:p) {
            s2[j] = 1/(M.inv[[i]][j,j] + Tau.inv[[G[i]]][j,j] + B.Sinv.B.j[j])
            zi.star = zi # p+m
            zi.star[j] = xx[i,j]
            mui.star = mui # p
            mui.star[j] = mu[G[i],j]
            pred = matrix(X[i,-(pin+j)], 1, pin+p-1) %*% matrix(B[-(pin+j),], pin+p-1, m) # 1*m; "matrix" is needed to avoid "non-conformable arguments" error when p=1, m>1
            xi.hat[j] = s2[j] * (dot(M.inv[[i]][j,], zi.star) + dot(Tau.inv[[G[i]]][j,], mui.star) + dot(B.Sinv[j,], Y[i,]-pred))
          } # next covariate
          X[i,prange] = rnorm(p, xi.hat, sqrt(s2))
        } # next data point
      }
    }
    ##
    ## update the responses
    if (update.Y) {
      pred = X %*% B # n*m
      q = pred - Y # n*m
      eta.hat = rep(NA, m)
      for (i in 1:n) {
        ## Stripped-down pin=1, p=1, m=1 version for sanity checking
        ## rho = M[[i]][1,2]/sqrt(M[[i]][1,1]*M[[i]][2,2])
        ## s2 = 1/( 1/(M[[i]][2,2]*(1-rho^2)) + 1/Sigma[1,1] )
        ## thing = yy[i,1] + M[[i]][1,2]*(X[i,2] - xx[i,1]) / M[[i]][1,1]
        ## eta.hat = s2 * ( thing/(M[[i]][2,2]*(1-rho^2)) + (B[1,1]+B[2,1]*X[i,2])/Sigma[1,1] )
        ## Y[i,1] = rnorm(1, eta.hat, sqrt(s2))
        ## .. and now the real thing:
        ## if p=0, this should properly skip the xx and X below
        zi = c(xx[i,], yy[i,]) - c(X[i,prange], Y[i,]) # NB skipping of intercept column in X, if any; p+m
        s2 = 1/(diag(M.inv[[i]])[p+1:m] + diag(Sigma.inv)) # m
        ## note: I can't see a way of collapsing this further
        for (j in 1:m) {
          zi.star = zi # p+m
          zi.star[p+j] = yy[i,j]
          qi.star = q[i,] # m
          qi.star[j] = pred[i,j]
          eta.hat[j] = s2[j] * (dot(M.inv[[i]][p+j,], zi.star) + dot(Sigma.inv[j,], qi.star)) # scalar (index over m)
        } # next response j
        Y[i,] = rnorm(m, eta.hat, sqrt(s2))
      } # next data point
    }
    ##
    ## update the Dirichlet process concentration parameter
    if (update.alpha) {
      ln.eta = log(rbeta(1, dp.concentration+1, n))
      pi.eta = (dp.prior.alpha + nclusters - 1) / (dp.prior.alpha + nclusters - 1 + n * (dp.prior.beta - ln.eta))
      if (runif(1) <= pi.eta) {
        alpha = dp.prior.alpha + nclusters
      } else {
        alpha = dp.prior.alpha + nclusters - 1
      }
      beta = dp.prior.beta - ln.eta
      dp.concentration = rgamma(1, alpha, rate=beta)
    }
    ## update the mixture labels
    if (update.G) {
      q = array(dim=Ngauss)
      for (i in 1:n) {
        for (k in 1:Ngauss) q[k] = ppii[k] * dmvnorm(X[i,prange], mu[k,], Tau[[k]])
        ## q = q / sum(q)   -  not necessary; see help for rmultinom
        G[i] = which(rmultinom(1, 1, q) == 1)
      }
      for (k in 1:Ngauss) nG[k] = length(which(G == k))
    }
    ##
    ## update the mixture component means
    if (update.mu) {
      for (k in 1:Ngauss) {
        if (dirichlet) { # NB Ngauss=1 in this case
          nk = nclusters
          gg = uniqueX
        } else {
          nk = nG[k]
          gg = which(G==k)
        }
        S.mu = symmsolve(U.inv + nk*Tau.inv[[k]])
        if (nk == 0) {
          mu.hat = S.mu %*% (U.inv %*% mu0)
        } else if (p == 1) {
          mu.hat = S.mu %*% (U.inv %*% mu0 + nk * Tau.inv[[k]] * mean(X[gg, pin+1]))
        } else {
          mu.hat = S.mu %*% (U.inv %*% mu0 + nk * Tau.inv[[k]] %*% colMeans(matrix(X[gg, prange], nrow=length(gg))))
        }
        mu[k,] = rmvnorm(1, mu.hat, S.mu)
      }
    }
    ##
    ## update mean of mixture component means
    if (update.mu0) mu0 = as.vector(rmvnorm(1, colMeans(mu), U/Ngauss))
    ##
    ## update covariance of mixture component means
    if (update.U) {
      S = W
      for (k in 1:Ngauss) {
        z = mu[k,] - mu0
        S = S + z %*% t(z)
      }
      U.inv = rwishart(Ngauss + p, symmsolve(S))
      U = symmsolve(U.inv)
    }
    ##
    ## update hyperprior of mixture component covariances
    if (update.W) {
      S = U.inv
      for (k in 1:Ngauss) S = S + Tau.inv[[k]]
      W = rwishart((Ngauss+2)*p+1, symmsolve(S))
    }
    ## save this sample
    if (return.X) res$X[,,Isample] = X
    if (return.Y) res$Y[,,Isample] = Y
    if (return.B) res$B[,,Isample] = B
    if (return.Sigma) res$Sigma[,,Isample] = Sigma
    if (return.G) res$G[,Isample] = G
    if (return.pi) res$pi[,Isample] = ppii
    if (return.mu) res$mu[,,Isample] = mu
    if (return.Tau) for (k in 1:Ngauss) res$Tau[,,k,Isample] = Tau[[k]]
    if (return.mu0) res$mu0[,Isample] = mu0
    if (return.U) res$U[,,Isample] = U
    if (return.W) res$W[,,Isample] = W
    if (return.alpha) res$alpha[Isample] = dp.concentration
    ##
    if (!is.na(mention.every) & Isample%%mention.every==0) message(paste('Done with Gibbs iteration', Isample))
    if (!is.na(save.every) & !is.na(save.to) & Isample%%save.every==0) save(res, file=save.to)
  } # end of Gibbs loop
  res
}


## casts the output of Gibbs.regression as a data frame
## todo: could be an object-oriented specialization of as.data.frame instead
Gibbs.post2dataframe = function(p) {
  if (length(p) == 0) {
    f = NULL
  } else {
    n = dim(p[[1]])
    n = n[length(n)]
    f = data.frame(row.names=1:n)
    for (n in names(p)) {
      d = dim(p[[n]])
      if (length(d) == 1) {
        f[,n] = p[[n]]
      } else {
        if (length(d) == 2) {
          for (i in 1:d[1]) f[,paste(n, i, sep='')] = p[[n]][i,]
        } else {
          if (length(d) == 3) {
            for (i in 1:d[1]) {
              for (j in 1:d[2]) {
                if (j < i & (n=='Sigma' | n=='U' | n=='W')) next # symmetric matrices
                f[,paste(n, i, j, sep='')] = p[[n]][i,j,]
              }
            }
          } else {
            ## remaining case is Tau (4D, 1st 2 dimensions are a symmetric matrix)
            for (i in 1:d[1])
              for (j in i:d[2])
                for (k in 1:d[3])
                  f[,paste(n, i, j, k, sep='')] = p[[n]][i,j,k,]
          } # 3D
        } # 2D
      } # 1D
    } # next list item
  } # p NULL
  f[which(!is.na(f[,1])),]
}
