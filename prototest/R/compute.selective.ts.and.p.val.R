compute.selective.ts.and.p.val <-
function(y, x, type, A, b, mu, sigma){
  if (type == "ELR"){
    ## compute the test statistic
    U = svd(x)$u
    ts = rcpp_compute_lr_stat (U=U, y=as.matrix(y), mu=ifelse(is.null(mu), Inf, mu), sigma=ifelse(is.null(sigma), Inf, sigma), exact=TRUE, verbose=FALSE, tol=10^-8, maxit=10000)
    
    ## find the limits
    limits = limits.exact.lr(U=U, y=y, A=A, b=b, mu=mu, sigma=sigma, M=ncol(U))
    
    ## compute the p-value
    p.Q = pchisq (limits[2], df=1, log.p=TRUE)
    p.ts = pchisq (ts, df=1, log.p=TRUE)
    p.q = pchisq (limits[1], df=1, log.p=TRUE)
    
    p.val = exp(log(1 - exp(p.ts-p.Q)) - log (1 - exp(p.q-p.Q)))
    #p.val = pchisq (ts, df=1, lower.tail=FALSE)
    return (list(ts=ts, p.val=p.val))
  }
  if (type == "ALR"){
    ## compute the test statistic
    U = svd(x)$u
    M = ncol(U)
    ts = rcpp_compute_lr_stat (U=U, y=as.matrix(y), mu=ifelse(is.null(mu), Inf, mu), sigma=ifelse(is.null(sigma), Inf, sigma), exact=FALSE, verbose=FALSE, tol=10^-8, maxit=10000)
    
    ## find the limits
    limits = limits.approx.lr(U=U, y=y, A=A, b=b, mu=mu, sigma=sigma)
    
    ## compute p-value
    p.large.top = pchisq (min(M + sigma*sqrt(2*M*ts), limits[2]), df=M, log.p=TRUE)
    p.small.top = pchisq (max(M - sigma*sqrt(2*M*ts), limits[1]), df=M, log.p=TRUE)
    p.large.bottom = pchisq (limits[2], df=M, log.p=TRUE)
    p.small.bottom = pchisq (limits[1], df=M, log.p=TRUE)
    
    p.val = 1 - exp(p.large.top + log (1-exp(p.small.top-p.large.top)) - p.large.bottom - log (1-exp(p.small.bottom-p.large.bottom)))
    return (list(ts=ts, p.val=p.val))
  }
  if(type == "F"){
    QRS.obj = compute.QRS.vectors (y=y, X.E=x, A=A, b=b, mu=mu) # test stat and per line intervals
    F.trunc.interval = find.overall.truncation.interval(QRS.obj$Q, QRS.obj$R, QRS.obj$S, verbose=FALSE)/QRS.obj$c # overall interval
    F.test.p.val = compute.trunc.F.test.p.value (QRS.obj$F.stat, F.trunc.interval, df1=QRS.obj$df1, df2=QRS.obj$df2) # p-value
    return (list(ts=QRS.obj$F.stat, p.val=F.test.p.val))
  }
  if (type == "MS"){
    # test statistic
    ts = as.numeric(compute.mc.t.statistic(x=x, y=y, mu=mu, sigma=sigma))
    
    # bounds
    if (is.null(mu)){
      x.tilde = x - mean(x)
      z = x.tilde/sqrt(sum(x.tilde^2))
    }else{
      z = x/sqrt(sum(x^2))
    }
    limits = sqrt(sum(x^2))*find.limits(y=y, u=z, A=A, b=b)
    
    # p-value
    p.val = compute.mc.t.p.val (ts, limits)
    return (list(ts=ts, p.val=p.val))
  }
}
