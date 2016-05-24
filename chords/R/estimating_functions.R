estimate.b.k.2 <- function(k, A.k, B.k, n.k, n.k.count, k.ind){
  ## Initialize:
  result <- list(N.k=NA, converge=FALSE)
  N.k <- NA
  
  target <- function(N.k){
    const1 <- N.k-(B.k/A.k)
    pre.const2 <- (N.k - n.k[which(k.ind)-1])
    
    # Deal with impossible N.k:
    if(any(pre.const2<0)) return(-Inf)
    
    const2 <- sum(1/pre.const2)
    const1*const2 - n.k.count
  }
  
  roots <- NULL
  #   xs <- 0:1000; ys <- sapply(xs, target);   range(ys[is.finite(ys)])
  #   plot(y=ys,x=xs, type='h');abline(0,0)
  #   target(n.k.count)
  #   target(n.k.count*100)
  .interval <- n.k.count*c(1, 10*max(1,1/A.k))
  try(roots <- uniroot(f = target, interval =.interval), silent = TRUE)
  
  if(length(roots)>1) {
    result$N.k <- roots$root
    result$converge <- 0
  } else {
    result$N.k <- n.k.count
    result$converge <- sign(target(n.k.count))
  }
  
  result$N.k <- as.integer(ceiling(result$N.k))
  return(result)
}
## Testing:
# save(n.k, k.ind,A.k, B.k, file='temp/testing_setup.RData')
# load(file='temp/testing_setup.RData')
# estimate.b.k.2(k = 3, A.k = A.k, B.k = B.k, n.k =n.k, k.ind = k.ind )






# estimate beta_k from sampled degrees and snowball matrix:
estimate.b.k<- function (rds.object, 
                         const=1, 
                         impute.Nks=TRUE) {
  ### Sketch:
  # Generate estimable parameters vector.
  # Optimized parameter-wise.
  
  
  ### Verifications:
  if(length(rds.object$estimates)>0) {
    message('Overwriting existing estimates in rds.object.')  
  }
  
  ### Initialize:
  arrival.times <- formatArrivalTimes(rds.object$rds.sample$interviewDt)
  arrival.intervals <- diff(arrival.times)
  
  arrival.degree<- rds.object$rds.sample$NS1
  max.observed.degree<- max(arrival.degree)
  ## TODO: should the degree of the first (seed) be removed?
  degree.counts<- table(arrival.degree)
  
  # Sequences per degree
  I.t <- rds.object$I.t
  degree.in <- rds.object$degree.in
  degree.out <- rds.object$degree.out
  
  likelihood <- NA
  
  ### Estimate:
  Nk.estimates<- rep(9999L, max.observed.degree) 
  names(Nk.estimates)<- seq_len(max.observed.degree)
  log.bk.estiamtes<- rep(NA, max.observed.degree) 
  A.ks<- rep(NA, max.observed.degree) 
  B.ks<- rep(NA, max.observed.degree) 
  n.k.counts<- rep(NA, max.observed.degree)
  convergence<- rep(NA, max.observed.degree)
  
  uniques<- as.integer(names(degree.counts))
  Nk.estimates[-uniques]<- 0
  for(k in uniques){
    # k <- uniques[[1]]
    k.ind <- arrival.degree==k
    k.ind[1] <- FALSE # dealing with sample kickoff
    
    n.k <- cumsum((arrival.degree==k))
    #     n.k <- cumsum((degree.in==k) - (degree.out==k))
    n.k.count <- degree.counts[paste(k)]
    
    #     head(cbind(arrival.degree, n.k, I.t, arrival.intervals),20)
    #     head(cbind(arrival.degree[-1], n.k[-1], I.t[-1], arrival.intervals))
    #     head(cbind(arrival.degree[-1], head(n.k,-1), head(I.t,-1), arrival.intervals))
    
    A.k <- sum( head(I.t,-1) * arrival.intervals, na.rm=TRUE)
    B.k <- sum( head(I.t,-1) * arrival.intervals * head(n.k,-1), na.rm=TRUE)    
    
    .temp <- estimate.b.k.2(k=k, A.k=A.k*const, B.k=B.k*const, n.k=n.k, 
                            n.k.count= n.k.count, k.ind=k.ind)    
    
    Nk.estimates[k]<-.temp$N.k
    log.bk.estiamtes[k] <- log(n.k.count) - log(.temp$N.k *  A.k - B.k)
    A.ks[k] <- A.k
    B.ks[k] <- B.k
    n.k.counts[k] <- n.k.count
    convergence[k] <- .temp$converge
  }
  
  if(impute.Nks) {
    Nk.estimates <- imputeEstimates(Nk.estimates, n.k.counts, convergence)
  }

  likelihood.val <- likelihood(log.bk = log.bk.estiamtes, 
                           Nk.estimates = Nk.estimates, 
                           I.t = I.t, 
                           n.k.counts = n.k.counts, 
                           degree.in = degree.in , 
                           degree.out = degree.out , 
                           arrival.intervals = arrival.intervals, 
                           arrival.degree = arrival.degree)
  
  result<- list(
    call=sys.call(),
    Nk.estimates=Nk.estimates, 
    log.bk.estimates=log.bk.estiamtes,
    A.ks=A.ks,
    B.ks=B.ks,
    n.k.counts=n.k.counts,
    arrival.intervals=arrival.intervals,
    arrival.degree=arrival.degree, 
    convergence=convergence,
    likelihood=likelihood.val)
  
  return(result)    					
}






