compute.test.statistic <-
function(x, y, U, type, mu, sigma, verbose, tol){
  if (type == "ELR"){
    U = svd(x)$u
    return(list(ts=rcpp_compute_lr_stat (U=U, y=y, mu=ifelse(is.null(mu), Inf, mu), sigma=ifelse(is.null(sigma), Inf, sigma), exact=TRUE, verbose=verbose, tol=tol, maxit=10000)))
  }
  if (type == "ALR"){
    U = svd(x)$u
    return(list(ts=rcpp_compute_lr_stat (U=U, y=y, mu=ifelse(is.null(mu), Inf, mu), sigma=ifelse(is.null(sigma), Inf, sigma), exact=FALSE, verbose=verbose, tol=tol, maxit=10000)))
  }
  if (type == "F"){
    ts = compute.F.statistic (x=x, y=y, mu=mu)
    #print (ts)
    return(ts)
  }
  if (type == "MS"){
    return(compute.mc.t.statistic (x=x, y=y, mu=mu, sigma=sigma))
  }
}
