gtcorr.hierarchical <- function(n, m=1, p, sigma=0, se=1, sp=1, arrangement=c("nested", "random"), model=c('beta-binomial', 'Madsen', 'Morel-Neerchal'), ...) {
  if (any(c(p, sigma, se, sp)<0) || any(c(p, sigma, se, sp) >1)) {
    stop("p, sigma, se, and sp should be between 0 and 1")
  }
  if (n[length(n)] != 1) {
    n <- c(n,1)
  }
  if (length(n)==1) {
    stop("The largest pool size should be greater than 1")
  }
  if (!all(n==as.integer(n))) {
    stop("All elements of n should be integers")
  }
  if (!all(0==n[-length(n)]%%n[-1])) {
    stop("n[s-1] should be divisible by n[s] for all s in 2:length(n)")
  }
  if (any(m>n[1])) {
    stop("m should not be greater than n[1]")
  }
  if (sum(c(length(m),
            length(p),
            length(sigma),
            length(se),
            length(sp))>1)>1) {
    stop("Only one of m, p, sigma, se, and sp can have a length > 1")
  }
  h <- length(n)
  model=match.arg(model)
  arrangement=match.arg(arrangement)
   
  w <- n[1]/n
  
  param.grid <- expand.grid(p=p, sigma=sigma, se=se, sp=sp, m=m)
  eff <- mapply(gtcorr.h.internal,
                p=param.grid$p,
                sigma=param.grid$sigma,
                se=param.grid$se,
                sp=param.grid$sp,
                m=param.grid$m,
                MoreArgs=list(n=n,
                              h=h,
                              w=w,
                              model,
                              arrangement,
                              ...)
                )
  
  result=list(n=n, param.grid=param.grid, arrangement=arrangement, model=model, efficiency=eff)
  result
}


gtcorr.h.internal <- function(p, sigma, se, sp, n, m, h, w, model, arrangement, ...) {

  #Calculate pr(V^T_s1=0) and pr(V^T_s1=1)
  prob.neg <- get(paste("prob.neg.h.", arrangement, sep=""))
  pr.v.t.0 <- prob.neg(p, sigma, n, m, h, w, model, ...)
  #Probabilities that a pool has at least one positive unit for each stage
  pr.v.t.1 <- 1 - pr.v.t.0

  #Calculate pr(V^T_(s-1)=u|V^T_s=0) for s=2, ..., h-1
  #(pr(true state of parent pool given true state of child))
  pr.par.t.giv.ch.t.0 <- list()
  #pr.par.t.giv.ch.t[[s]][u]=pr(V^T_(s-1)1=u|V^T_s1=0)
  #first element [[1]] is NULL
  if (h > 2) {
    for (s in 2:(h-1)) {
      pr.par.t.giv.ch.t.0[[s]] <- rep(NA, 2)
      pr.par.t.giv.ch.t.0[[s]][1] <- pr.v.t.0[s-1]/pr.v.t.0[s]
      pr.par.t.giv.ch.t.0[[s]][2] <- (pr.v.t.1[s-1]-pr.v.t.1[s])/pr.v.t.0[s]
    }
  }

  #Calculate pr(V_(s-1)=1|V^T_s=v) and pr(V_(s-1)=1|V^T_(s-1)=u
  pr.par.1.giv.ch.t <- list()
  #pr.par.1.giv.ch.t[[s]][v+1]=pr(V_(s-1)=1|V^T_s=v)

  if (h > 2) {
    pr.par.1.giv.ch.t <- list()
    pr.par.1.giv.ch.t[[2]] <- rep(NA, 2)
    pr.par.1.giv.ch.t[[2]][1] <- se    *pr.par.t.giv.ch.t.0[[2]][2]+
                                 (1-sp)*pr.par.t.giv.ch.t.0[[2]][1]
    pr.par.1.giv.ch.t[[2]][2] <- se
  }
  
  if (h > 3) {
    for (s in 3:(h-1)) {
      pr.par.1.giv.ch.t[[s]] <- rep(NA, 2)
      pr.par.1.giv.ch.t[[s]][1] <- se    *pr.par.1.giv.ch.t[[s-1]][2]*pr.par.t.giv.ch.t.0[[s]][2]+
                                   (1-sp)*pr.par.1.giv.ch.t[[s-1]][1]*pr.par.t.giv.ch.t.0[[s]][1]
      pr.par.1.giv.ch.t[[s]][2] <- se*pr.par.1.giv.ch.t[[s-1]][2]
    }
  }

  #Calculate pr(V_s=1)
  pr.v.1 <- rep(NA, h-1)
  pr.v.1[1] <- se*pr.v.t.1[1]+(1-sp)*pr.v.t.0[1]
  if (h > 2) {
    for (s in 2:(h-1)) {
      pr.v.1[s] <- se    *pr.par.1.giv.ch.t[[s]][2]*pr.v.t.1[s] +
                   (1-sp)*pr.par.1.giv.ch.t[[s]][1]*pr.v.t.0[s]
    }
  }

  et <- c(1, w[2:(h)]*pr.v.1[1:(h-1)])

  sum(et)/n[1]
}

prob.neg.h.nested <- function(p, sigma, n, m, h, w, model, ...) {
  #returns a vector of length h-1 where the sth element is the probability that a pool at
  #stage s has no positive units
  h.prime <- which.max(n<=m)
  
  if (c(!all(n[1:(h.prime-1)]%%m==0, m%%n[h.prime:h]==0))) {
    stop(paste("n[s] should be divisible by m for s < h' and n[s] should be divisible by m for s>=h'. h'=",
               h.prime, " in this case", sep=""))
  }

  #Vector of probabilities that a pool has no positive units for each stage
  #pr.v.t.0 <- rep(NA, h)  
  pr.v.t.0 <- rep(NA, h-1)
  
  
  #Probability that a cluster has no positive units for cluster sizes up to m
  q <- sapply(1:m, function(m) prob.neg.clust(p, sigma, m, model))
  
  for (s in 1:(h.prime-1)) {
    pr.v.t.0[s] <- q[m]^(n[s]/m)
  }
#  for (s in h.prime:h) {
#    pr.v.t.0[s] <- q[n[s]]
#  }  
  if (h.prime <= h-1) {
    for (s in h.prime:(h-1)) {
      pr.v.t.0[s] <- q[n[s]]
    }
  }
  
  pr.v.t.0
}

prob.neg.h.random <- function(p, sigma, n, m, h, w, model, ...) {
  #returns a vector of length h where the sth element is the probability that a pool at
  #stage s has no positive units
  pr.v.t.0 <- rep(NA, h-1)

  for (s in 1:(h-1)) {
    pr.v.t.0[s] <- prob.neg.pool.rand(p=p, sigma=sigma, m=m, size=n[s], k=n[1]/m, model=model, ...)
  }
  pr.v.t.0
}

gtcorr.hierarchical.user <- function(n,clusters,p,sigma=0,se=1,sp=1, model=c("beta-binomial", "Morel-Neerchal", "Madsen")){
  if (any(c(p, sigma, se, sp)<0) || any(c(p, sigma, se, sp) >1)) {
    stop("p, sigma, se, and sp should be between 0 and 1")
  }
  if (n[length(n)] != 1) {
    n <- c(n,1)
  }
  if (length(n)==1) {
    stop("The largest pool size should be greater than 1")
  }
  if (!all(n==as.integer(n))) {
    stop("All elements of n should be integers")
  }
  if (!all(0==n[-length(n)]%%n[-1])) {
    stop("n[s-1] should be divisible by n[s] for all s in 2:length(n)")
  }

  bad <- FALSE
  l <- max(clusters)
  if (l != as.integer(l)) bad <- TRUE
  if (l < 1) bad <- TRUE
  for (i in 1:l) {
    if (!(i %in% clusters)) bad <- TRUE
  }
  if (bad) stop("clusters should be a vector of integers from 1 to the largest cluster number")

  if (length(clusters) != n[1]) stop("clusters should have length n[1]")

  if (length(p)==1) p <- rep(p, l)
  if (length(p)!=l) stop("p should have length 1 or the largest cluster number")
  if (length(sigma)==1) sigma <- rep(sigma, l)
  if (length(sigma)!=l) stop("sigma should have length 1 or the largest cluster number")

  model <- match.arg(model)
  
  h <- length(n)	# total number of stages (scalar)
  n1 <- length(clusters)# total number of units (scalar)
  w <- n1/n		# total number of possible pools in stage s (vector)
  et <- 1		# expected number of tests at stage 1
  for (ss in 2:h){	# cycle through stages. see Section 3.1 of revision.
    summ <- 0
    for (ii in 1:w[ss-1]) summ <- summ + gtcorr.h.prv(ii,ss-1,n,clusters,p,se,sp,sigma, model)
    etss <- w[ss]/w[ss-1]*summ
    et <- et + etss	# add number of tests expected at stage ss
  }
  eff <- et/n1	# final efficiency 
}

gtcorr.h.prv <- function(ii,ss,n,M,prev,se,sp,sigma, model){
		# prob pool ii in stage ss tests positive. cf equation (4) and (5) of revision
  ans.prv <- 0
  if (ss==1){
  # cf equation (4) revision
    prv1iteq1 <- gtcorr.h.prvsiteq1(ii,1,n,M,prev,sigma, model)
    ans.prv <- se*prv1iteq1 + (1-sp)*(1-prv1iteq1)
  }
  if (ss>1){		
    # cf equation (5) of revision
    pt <- gtcorr.h.prvsiteq1(ii,ss,n,M,prev,sigma, model)
    ps <- c(1-pt,pt)
    for (vv in 0:1) {
      ans.prv <- ans.prv + se^vv * (1-sp)^(1-vv) * gtcorr.h.cp1(vv,ii,ss,n,M,prev,se,sp,sigma, model) * ps[vv+1]
    }
  }
  ans.prv
}	

gtcorr.h.prvsiteq1 <- function(ii,ss,n,M,prev,sigma, model){	# equation (3) revision
  # probability pool ii in stage ss is truly positive
  # prev vector (of length l) of cluster specific prevalences
  # msik vector (of length l) with number of units from cluster k in pool ii in stage ss
  l <- length(prev)
  msik <- matrix(NA,1,l)
  lower <- n[ss]*(ii-1)+1	# first unit in pool ii in stage ss
  upper <- lower+n[ss]-1
  index <- lower:upper # units in pool ii in stage ss
  for (kk in 1:l) msik[kk] <- sum(M[index]==kk)
  myprod <- 1
  for (kk in 1:l){
    myprod <- myprod*prob.neg.clust(prev[kk], sigma[kk], msik[kk], model)
  }
  ans <- 1-myprod
  ans
}

gtcorr.h.cp1 <- function(vv,ii,ss,n,M,prev,se,sp,sigma, model) {
  # conditional probability in equation (5) of revision
  # uses recursion
  if (ss==2){		# eq (5.1) of revision
    if (vv==1) ans <- se
    if (vv==0){
      ans <- 0
      for (uu in 0:1) ans <- ans + se^uu * (1-sp)^(1-uu) * gtcorr.h.cpt(uu,ii,ss,n,M,prev,sigma, model)
        # here cpt is conditional prob on right side of equation (5.1) of revision
    }  
  }
  else{			# eq (5.2) of revision
    iim1 <- ceiling(ii/(n[ss-1]/n[ss]))	# pool in stage ss-1 containing units from pool ii in stage ss
    if (vv==1) ans <- se*gtcorr.h.cp1(1,iim1,ss-1,n,M,prev,se,sp,sigma, model)	# recursion
    if (vv==0){
      ans <- 0
      for (uu in 0:1) ans <- ans + se^uu * (1-sp)^(1-uu) * gtcorr.h.cpt(uu,ii,ss,n,M,prev,sigma, model)* 
        gtcorr.h.cp1(uu,iim1,ss-1,n,M,prev,sigma, model)	# recursion
    }
  }
  ans
} # end cp1	

gtcorr.h.cpt <- function(uu,ii,ss,n,M,prev,sigma,model){
		# conditional probability on right side of equation (5.1) and (5.2) of revision only involving truth
  iim1 <- ceiling(ii/(n[ss-1]/n[ss]))	# pool in stage ss-1 containing units from pool ii in stage ss
  num <- 1-gtcorr.h.prvsiteq1(iim1,ss-1,n,M,prev,sigma, model)
  den <- 1-gtcorr.h.prvsiteq1(ii,ss,n,M,prev,sigma, model)
  if (uu==0) ans1 <- num/den
  else ans1 <- 1-num/den
  ans1
}


###End Hierarchical functions###

###Matrix functions###

gtcorr.matrix <- function(r, c, m=1, p, sigma=0, se=1, sp=1, r.prime, c.prime, arrangement=c("rectangular", "diagonal", "random"), model=c('beta-binomial', 'Madsen', 'Morel-Neerchal'), ...) {
  if (any(c(p, sigma, se, sp)<0) || any(c(p, sigma, se, sp) >1)) {
    stop("p, sigma, se, and sp should be between 0 and 1")
  }
  if (any(r*c %% m !=0)) {
    stop("r*c should be divisible by m")
  }
  arrangement <- match.arg(arrangement)
  model <- match.arg(model)
  if (arrangement == "diagonal" &&
      (r != c || r != m)) {
    stop("For a diagonal arranement, r, c, and m should be equal.")
  }
  if (sum(c(length(p),
            length(sigma),
            length(se),
            length(sp))>1)>1) {
    stop("Only one of p, sigma, se, and sp can have a length > 1")
  }
  if (length(m)>1) {
    stop("m should have length 1")
  }

  
  param.grid <- expand.grid(p=p, sigma=sigma, se=se, sp=sp)
  eff <- mapply(gtcorr.m.internal,
                p=param.grid$p,
                sigma=param.grid$sigma,
                se=param.grid$se,
                sp=param.grid$sp,
                MoreArgs=list(r=r,
                              c=c,
                              m=m,
                              r.prime=r.prime,
                              c.prime=c.prime,
                              model,
                              arrangement,
                              ...)
                )
  
  result=list(r=r, c=c, m=m, r.prime=r.prime, c.prime=c.prime, param.grid=param.grid, arrangement=arrangement, model=model, efficiency=eff)
  result
}

gtcorr.m.internal <- function(p, sigma, se, sp, m, r, c, r.prime, c.prime, model, arrangement, ...) {

  #Get pr(R_i^T=0), pr(C_j^T=0) and pr(R_i^T=C_j^T=0)
  prob.neg <- get(paste("prob.neg.m.", arrangement, sep=""))
  probs <- prob.neg(p=p, sigma=sigma, r=r, c=c, m=m, r.prime=r.prime, c.prime=c.prime, model=model, ...)
  pr.rt.0 <- probs$pr.rt.0
  pr.ct.0 <- probs$pr.ct.0
  pr.rt.ct.0 <- probs$pr.rt.ct.0
  

  #The u+1th, v+1th element of pr.rt.ct is pr(R_i^T=u, C_j^T=v)
  pr.rt.ct <- matrix(NA, 2, 2)
  pr.rt.ct[1, 1] <- pr.rt.ct.0  
  pr.rt.ct[2, 1] <- pr.ct.0 - pr.rt.ct.0
  pr.rt.ct[1, 2] <- pr.rt.0 - pr.rt.ct.0
  pr.rt.ct[2, 2] <- 1 - (pr.rt.0 + pr.ct.0 - pr.rt.ct.0)
  
  pr.r.c.1 <- 0
  for (u in 0:1) {
    for (v in 0:1) {
      pr.r.c.1 <- pr.r.c.1 + se ^ (u + v) * (1 - sp) ^ (2 -u - v) * pr.rt.ct[u+1, v+1]
    }
  }
  et <- r+c+r*c*pr.r.c.1
  et/(r*c)
}

prob.neg.m.rectangular <- function(p, sigma, r, c, m, r.prime, c.prime, model, ...) {
  pr.rt.0 <- prob.neg.clust(p, sigma, c.prime, model) ^ (c / c.prime)
  pr.ct.0 <- prob.neg.clust(p, sigma, r.prime, model) ^ (r / r.prime)
  pr.rt.ct.0 <- prob.neg.clust(p, sigma, c.prime + r.prime - 1, model) *
                prob.neg.clust(p, sigma, c.prime, model) ^ (c / c.prime - 1) *
                prob.neg.clust(p, sigma, r.prime, model) ^ (r / r.prime - 1)
  list(pr.rt.0=pr.rt.0, pr.ct.0=pr.ct.0, pr.rt.ct.0=pr.rt.ct.0)
}

prob.neg.m.diagonal <- function(p, sigma, r, model, ...) {
  #For a diagonal arrangment, r=c=m so can ignore other parameters
  pr.rt.0 <- (1 - p) ^ (r)
  pr.ct.0 <- pr.rt.0
  pr.rt.ct.0 <- (1 - p) * prob.neg.clust(p, sigma, 2, model) ^ (r - 1)
  list(pr.rt.0=pr.rt.0, pr.ct.0=pr.ct.0, pr.rt.ct.0=pr.rt.ct.0)
}

prob.neg.m.random <- function(p, sigma, r, c, m, model, ...) {
  pr.rt.0 <- prob.neg.pool.rand(p=p, sigma=sigma, m=m, size=c, k=r*c/m, model=model, ...)
  if (r == c) { #If the matrix is square, pr(R_i^T=0)=pr(C_j^T=0)
    pr.ct.0 <- pr.rt.0
  } else {
    pr.ct.0 <- prob.neg.pool.rand(p=p, sigma=sigma, m=m, size=r, k=r*c/m, model=model, ...)
  }
  pr.rt.ct.0 <- prob.neg.pool.rand(p=p, sigma=sigma, m=m, size=r+c-1, k=r*c/m, model=model, ...)
  list(pr.rt.0=pr.rt.0, pr.ct.0=pr.ct.0, pr.rt.ct.0=pr.rt.ct.0)  
}

gtcorr.matrix.user <- function(clusters, p, sigma=0, se=1, sp=1, model=c('beta-binomial', 'Madsen', 'Morel-Neerchal')) {
  bad <- FALSE
  l <- max(clusters)
  if (l != as.integer(l)) bad <- TRUE
  if (l < 1) bad <- TRUE
  for (i in 1:l) {
    if (!(i %in% clusters)) bad <- TRUE
  }
  if (bad) stop("clusters should be a matrix of integers from 1 to the largest cluster number")

  if (length(p)==1) p <- rep(p, l)
  if (length(p)!=l) stop("p should have length 1 or the largest cluster number")
  if (length(sigma)==1) sigma <- rep(sigma, l)
  if (length(sigma)!=l) stop("sigma should have length 1 or the largest cluster number")
  if (any(c(p, sigma, se, sp)<0) || any(c(p, sigma, se, sp) >1)) {
    stop("p, sigma, se, and sp should be between 0 and 1")
  }

  model <- match.arg(model)
    
  r <- nrow(clusters)
  c <- ncol(clusters)
  et <- r+c
  for (ii in 1:r) {
    for (jj in 1:c) {
      ricj <- 0
      for (uu in 0:1) {
        for (vv in 0:1) {
          ricj <- ricj + se^(uu+vv) * (1-sp)^(2-uu-vv) * gtcorr.m.joint(ii,jj,uu,vv,clusters,p,sigma, model)
        }
      }
      et <- et + ricj
    }
  }
  eff <- et/(r*c)
  eff
}

gtcorr.m.joint <- function(ii,jj,uu,vv,clusters,p,sigma, model){
  pr0 <- pc0 <- prc0 <- 1
  l <- length(p)
  for (kk in 1:l){
    mik <- sum(clusters[ii,]==kk)	# number of units in row ii from cluster kk
    mjk <- sum(clusters[,jj]==kk)
    mijk <- sum(clusters[ii,]==kk) + sum(clusters[-ii,jj]==kk)
    pr0 <- pr0*prob.neg.clust(p[kk], sigma[kk], mik, model)
    pc0 <- pc0*prob.neg.clust(p[kk], sigma[kk], mjk, model)
    prc0 <- prc0*prob.neg.clust(p[kk], sigma[kk], mijk, model)
  }
  if (uu==1 & vv==1) ans <- 1-(pr0+pc0-prc0)
  if (uu==1 & vv==0) ans <- pc0-prc0
  if (uu==0 & vv==1) ans <- pr0-prc0
  if (uu==0 & vv==0) ans <- prc0
  ans
}

###End Matrix functions###

###General functions###

prob.neg.clust <- function(p, sigma, m, model) {
  #Returns the probability that a cluster of size m has no positive units
  q <- 1-p

  if (m==0) return(1)

  if (sigma==1) return(q)
    
  if (model=='beta-binomial') {
    gamma <- sigma/(1-sigma)
    f <- 1
    for (j in 0:(m-1)) {
      f <- f * #previous product
           #(1 + (q + gamma * j - 1) * (0 <= j & j <= m-1 ) ) * #second term in numerator if j between Xdot and m - 1
           #(1 + gamma * j ) ** -1 ; #denominator
           (1 + (q + gamma * j - 1) ) / #numerator
           (1 + gamma * j ) ; #denominator
    }
    return (f)
  }
  
  if (model=='Madsen') {
    return((1-sigma) * q ^ m + sigma * q)
  }
  
  if (model=='Morel-Neerchal') {
    return(p * ( q * ( 1 - sqrt(sigma))) ^ m + q * (q + p * sqrt(sigma))^m)
  }
}

prob.neg.pool.rand <- function(p, sigma, m, size, k, model, runs=1000, ...) {
  #Returns probability that a pool of size size has no positive units where
  #The k clusters of size m are arranged randomly within the pool

  #clust.idx is a vector of length m*k containing the indecies of each cluster
  #numbered 1 to k
  clust.idx <- as.factor(rep(1:k, each=m))

  #Pools is a matrix with 'runs' columns of lenght 'size'
  #Each column represents the the cluster indicies in a randomly arranged pool
  pools <- replicate(runs, sample(x=clust.idx,size=size), simplify=FALSE)
  counts <- lapply(pools, tabulate)

  #q is a vector of length m+1, where the m'+1th element is q_m', q_0=q[1]=1
  q=sapply(1:m, function(m) prob.neg.clust(p=p, m=m, sigma=sigma, model=model))
  q=c(1,q)

  #qpools is a vector with the probability that a pool tests negative for each run
  qpools <- sapply(counts, function(x) prod(q[x+1]))

  #return average
  mean(qpools)
}

###End General functions###
