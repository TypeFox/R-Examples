
simulateData <- function(d=2, k=4, n=100, m=10, mu=NULL, cvar=NULL, s.pi = rep(1/k, k), b.min=0.02, mean = "D", between = "D", within = "D", cov = "D", n.labels = k) {
  if (m>n) 
     stop("M (number of knowns) MUST be larger than n (number of all observations)")
  if (m<=0) 
     stop("M MUST be positive")
     
  Ytrue = c(sample.int(n.labels, m, T, prob=s.pi[1:n.labels]), sample.int(k, n-m, T, prob=s.pi))
  # means
  if (is.null(mu))
    mu = matrix(rnorm(k*d)*7, k, d)
  if (mean == "E")  
    mu[] = rep(colMeans(mu), each=k)
  # variances
  if (is.null(cvar)) {
    cvar = array(0, c(k, d, d))
    for (i in 1:k) {
        tmp = matrix(rnorm(2*d^2)*rexp(d),d,2*d)
        cvar[i,,] <- cov(t(tmp))
#       cvar[i,,] <- diag(d)*rexp(d)
    }
  }
  # are variances equal?
  if (between=="E") {
      # averaging among clusters
     ncvar = matrix(0, d, d)
     for (i in 1:k) 
        ncvar = ncvar + cvar[i, , ] 
     for (i in 1:k) 
        cvar[i, , ] = ncvar
  }
  if (within=="E" && d>1) {
      # averaging among variables
     for (i in 1:k) {
        ndiag = sum(diag(cvar[i, , ]))
        sdiag = ndiag/d
        noutd = min(sdiag, (sum(cvar[i, , ])-ndiag)/(d*(d-1)))
        cvar[i, , ] = noutd
        diag(cvar[i, , ]) = sdiag
     }
  }
  # are covariance equal to 0?
  if (cov=="0") {
   for (i in 1:k) 
      cvar[i, , ] = diag(diag(cvar[i, , ]), nrow=d)
  }
  
  X = matrix(0, n, d)
  for (i in unique(Ytrue)) {
    nobs = sum(Ytrue==i)
    X[Ytrue==i, ] = rmvnorm(nobs, mu[i,], as.matrix(cvar[i,,]))
  }
  B = matrix(b.min, m, min(k, n.labels))
  for (i in 1:m) 
     B[i,Ytrue[i]] = 1-b.min*(min(k, n.labels)-1)
  
  colnames(X) = paste("Col",1:ncol(X))
  model.params = list(pi = s.pi, mu = mu, cvar = cvar, m=m, n=n, k=k, d=d)
  res = list(X = X[-(1:m),,drop=F], knowns=X[1:m,,drop=F], B=B, model.params=model.params, Ytrue=Ytrue)
  class(res) = "simulated"
  res
}

