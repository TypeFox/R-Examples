# original version, with weighted Silverman bandwidth unique option
# has been improved (mid 2014) by a k-fold CV option
# the original version is kept as "npMSL_old" in this file
################################################################ 
################################################################ 
## nonparametric algorithm for Smoothed Likelihood Maximization 
## implementing block structure
## and Silverman adaptive bandwidth
## 2014 additions: loglik stores all the sequence,
## post argument for passing an init matrix of posterior
## bwiter argument for the duration of the adaptive bw stage
## can be set to 0 for keeping an initial bw matrix when samebw=FALSE
################################################################ 
################################################################ 
npMSL_old <- function(x, mu0, blockid = 1:ncol(x),
                  bw=bw.nrd0(as.vector(as.matrix(x))), samebw = TRUE,
                  h=bw, eps=1e-8,
                  maxiter=500, bwiter = maxiter,
                  ngrid=200, post=NULL, verb = TRUE){
  bw <- h # h is alternative bandwidth argument, for backward compatibility
  x <- as.matrix(x)
  n <- nrow(x)      # number of subjects
  r <- ncol(x)      # number of measurements in each subject
  u <- match(blockid, unique(blockid))
  B <- max(u) 		# nb of blocks
  BlS <- rep(0,B)	# block sizes = C_ell in JCGS paper
  for (ell in 1:B) {
      BlS[ell] <- sum(u == ell)}      
      
  if (is.matrix(mu0)) # mu0=centers
    m <- dim(mu0)[1]  else m <- mu0  # mu0=number of clusters
    
  if(!samebw && !is.matrix(bw)) { # create initial bandwidth matrix
    cat("create initial bandwidth matrix\n")
    bw <- matrix(bw, nrow=max(u), ncol=m)  
  }
  z.hat <- matrix(0, nrow = n, ncol = m)
  tt0 <-  proc.time() # for total time
  
  ## Initial Values
  if(m == 1) z.hat <- matrix(1, nrow = n, ncol = m) 
  else if(is.null(post)) {
    kmeans <- kmeans(x, mu0)
    for(j in 1:m)
      z.hat[kmeans$cluster==j, j] <- 1
  } 
  else { 
    z.hat <- post ## Currently no error-checking is done here
  }
  
  iter <- 0
  finished <- FALSE
  lambda <- matrix(0, nrow = maxiter, ncol = m)
  loglik <- NULL; loglikseq <- rep(NA,maxiter)
  total_udfl <- 0; total_nan <- 0 # eventual NaN and underflow in C code

  tmp <- 1:n   # is this needed?
  xtra <- (max(x)-min(x))/10
  grid <- seq(min(x)-xtra, max(x)+xtra, length=ngrid)
  # f stored on a ngrid by m by B array 
  # f_{g,j,ell} = f_{j ell}(u_g)
  # f <- array(1/m/diff(grid[1:2]), c(ngrid, m, B)) 
  # this f was not normalized for being uniform over grid 
  Delta <- diff(grid[1:2])  
  f <- array(1/((ngrid-1)*Delta), c(ngrid, m, B)) 
  oldloglik <- -Inf

  orderx <- xx <- list() # preparation for adaptive bandwidth
  for(k in 1:B) {
    xx[[k]] <- as.vector(x[, u==k])
    if (!samebw) {
      orderx[[k]] = order(xx[[k]]) # only needed for IQR calculation for bw
    }
  }
  CftEstep <- ifelse(samebw, "npMSL_Estep", "npMSL_Estep_bw")
  # CftEstep <- "npMSL_Estep_bw" # temporary, for testing only the M-step  
  CftMstep <- ifelse(samebw, "npMSL_Mstep", "npMSL_Mstep_bw")
  
  while (!finished) { # Algorithm main iteration loop
    iter <- iter + 1
    bw.old <- bw	# is this needed?
    t0 <- proc.time()
    nb_udfl=0;  # nb underflows, K()*log(0) ~0 cancelled in nems_Estep.c
    nb_nan=0;  # nb nonzero K()*log(0) cancelled in nems_Estep.c

    ## Note:  Enter loop assuming E-step is done -- i.e., z.hat in place    
    ## M-Step
    lambda[iter, ] <- colMeans(z.hat)
    
    ## density estimation in M-step: WKDE-step  
      cs <- colSums(z.hat)
      z.tmp <- sweep(z.hat, 2, cs, "/")
      z.tmp[, cs==0] <- 1/NROW(z.tmp) # Just in case

    ## adaptive bandwidth update IF in adaptive bw stage
    if (!samebw && iter <= bwiter) {
      for (ell in 1:B) {
        r2 <- BlS[ell]	# block size = nb of coordinates
        wts <- apply(z.tmp, 2, function(z) rep(z/r2, r2))
        variances <- colSums(wts*outer(xx[[ell]], colSums(wts*xx[[ell]]),'-')^2)
        iqr <- apply(as.matrix(wts[orderx[[ell]],]), 2, wIQR,
                     xx[[ell]][orderx[[ell]]],
                     already.sorted=TRUE, already.normalized=TRUE)
        h <- bw[ell, ] <- 0.9 * pmin(sqrt(variances), iqr/1.34) * 
                     pmax(1,r2*n*lambda[iter, ])^(-1/5) 
                     # Note:  Doesn't allow "sample size" < 1.
  #      browser()    
        } 
    } # end of bw adaptive stage

        
    z=.C(CftMstep, as.integer(ngrid), as.integer(n),
       as.integer(m), as.integer(r), 
       as.integer(B), as.integer(BlS), as.integer(u),
       as.double(bw), as.double(x), as.double(grid), 
       new.f=as.double(f), 
       as.double(lambda[iter,]), 
       as.double(z.hat))
       
    f <- array(z$new.f, c(ngrid, m, B)) # check sum(f == 0)
# print(apply(f,2:3,sum) * Delta)
# print(max(abs(f-f2)))
#    browser()

    ## E-step (for next iteration)
    z=.C(CftEstep, as.integer(ngrid), as.integer(n),
       as.integer(m), as.integer(r), 
       as.integer(B), as.integer(u),
       as.double(bw),
       as.double(x), as.double(grid), f=as.double(f),
       as.double(lambda[iter,]), post=as.double(z.hat),
       loglik = double(1),
       nb_udfl = as.integer(nb_udfl), nb_nan = as.integer(nb_nan))

    nb_udfl = z$nb_udfl; nb_nan = z$nb_nan; 
    total_udfl <- total_udfl + nb_udfl
    total_nan <- total_nan + nb_nan
    
    z.hat <- matrix(z$post, n, m)
    if (sum(is.nan(z.hat)) > 0) cat("Error!! NaN in z.hat") # obsolete ?
    loglik <- loglikseq[iter] <- z$loglik
    loglikchange <- loglik - oldloglik
    oldloglik <- loglik
    finished <- iter >= maxiter
    if (iter>1 && max(abs(lambda[iter, ] - lambda[iter-1, ])) < eps)
      finished <- TRUE
    if (verb) {
      t1 <- proc.time()
      cat("iteration", iter, ": lambda ", round(lambda[iter, ], 4))
      cat(" obj function", round(loglik, 4))
      cat(" (change", round(loglikchange,4), ")")
      cat(" time", (t1 - t0)[3])
      if ((nb_udfl > 0) || (nb_nan >0)) cat("\n  ")
      if (nb_udfl > 0) cat("average underflows=", round(nb_udfl/(n*m*r),3)," ")
      if (nb_nan >0) cat("average NaNs=", round(nb_nan/(n*m*r),3))
# Note: these average mean nb of nan over ngrid convolution
      cat("\n")
    }
  }
  # f <- array(z$f, c(ngrid, m, r)) # obsolete in block version

 if (!samebw) {
    rownames(bw) <- paste("block", 1:max(u))
    colnames(bw) <- paste("component", 1:m)
  	}

  if (verb) {
    tt1 <- proc.time()
    cat("lambda ", round(lambda[iter, ], 4))
    cat(", total time", (tt1 - tt0)[3], "s\n")
  }
  return(structure(list(data = x, posteriors = z.hat, 
                        lambda = lambda[1:iter,], bandwidth = bw, 
                        blockid = u, lambdahat = lambda[iter,], f=f,
                        grid = grid, loglik = loglikseq[1:iter],
			meanUdfl = total_udfl/(n*m*r*iter),# average underflow
			meanNaN = total_nan/(n*m*r*iter)), # average NaN's
                    class="npEM")) # define a "npMSL" class ?
}




## updated late 2014 version
##################################################################### 
##################################################################### 
## nonparametric algorithm for Smoothed Likelihood Maximization 
## implementing block structure
## 2014 latests additions: 
## option for Silverman &  weighted k-fold CV adaptive bandwidth:
##  parameter bwmethod = "S" for Silverman's rule (the default), 
##                       "CV" for Cross-Validation
## loglik now stores all the sequence
## post argument for passing an init matrix of posterior
## bwiter argument for the duration of the adaptive bw stage
## can be set to 0 for keeping an initial bw matrix when samebw=FALSE
###################################################################### 
###################################################################### 

## preliminary function for computing the weighted version of
## k-fold Cross-Validation

# Splitting 1:n in nbsets subsets and return indices
# n = total sample size
# setsize = size of sets
# nseq = explicit 1st and last indices
splitsample <- function(n, nbsets=2) {
  k <- floor(n/nbsets)
  klast <- n - (nbsets-1)*k
  setsize <- c(rep(k, nbsets-1), klast) # n per set
  if (sum(setsize) != n) {cat("warning")}
  ni <- c(1,cumsum(setsize)+1) # indices for splitting   
  nseq <- matrix(NA, nrow=nbsets, ncol=2)
  for (j in 1:nbsets) {
    nseq[j,1] <- ni[j]     # 1st index for jth set
    nseq[j,2] <- ni[j+1]-1 # last index for jth set
  }
  a = list(setsize=setsize, nseq=nseq)
  return(a)
}

## computes CV(h) in k-fold CV case
# h  = bandwidth (first argument for optimizing)
# x = sample of data, vector
# nbsets = number of folds
# w = weights
# lower, upper = numerical integration limits, data-driven defaults
kfoldCV <- function(h, x, nbsets=2, w = rep(1, length(x)), lower=mean(x)-5*sd(x), upper=mean(x)+5*sd(x)) {
  n <- length(x)
  fold <- splitsample(n, nbsets)
  sumf <- 0
  for (k in 1:nbsets) {
    Sk <- fold$nseq[k,1]:fold$nseq[k,2] # indices of set k 
    learnset <- x[-Sk] # obs not in set k, from which f_h is "learned"
    evalset <- x[Sk]
    fk <- wkde(x = learnset, u = evalset, w = w[-Sk], bw = h)
    sumf <- sumf + sum(fk)
  }
  integrand <- function(u,...) {wkde(x, u, bw=h)^2} # computing ||f_h||^2 
  fh2 <- integrate(integrand, lower=lower, upper=upper)$value
  return(fh2 - 2*sumf/n)
}


## Weighted bw selection by k-fold CV ####
## min search done by call to optimize
# x = vector of data
# w = weights, defaults to 1, unnormalized
wbw.kCV <- function(x, nbfold=5, w = rep(1, length(x)), hmin=0.1*hmax, hmax=NULL) {
  n <- length(x)
  if (is.null(hmax)) hmax <- 1.144 * sqrt(var(x))*n^(-1/5) # default hmax as in bw.ucv
  # computing lower and upper integration limits for ||fh||^2
  wm <- weighted.mean(x, w)
  lowerf <- wm - 5*sd(x); upperf <- wm + 5*sd(x) # maybe use a weighted.sd version as well?
  fucv <- function(h) kfoldCV(h, x, nbsets=nbfold, w=w, lower=lowerf, upper=upperf)
  hopt <- optimize(fucv, lower = hmin, upper = hmax)$minimum
  return(hopt)
}


##################################################################### 
##################################################################### 
## nonparametric algorithm for Smoothed Likelihood Maximization 
## implementing block structure
## 2014 latests additions: 
## option for Silverman &  weighted k-fold CV adaptive bandwidth:
##  parameter bwmethod = "S" for Silverman's rule (the default), 
##                       "CV" for Cross-Validation
## loglik now stores all the sequence
## post argument for passing an init matrix of posterior
## bwiter argument for the duration of the adaptive bw stage
## nbfold parameter passed to wbw.kCV, for leave-[n/nbfold]-out CV
## can be set to 0 for keeping an initial bw matrix when samebw=FALSE
###################################################################### 
###################################################################### 
# ToDo: add a nbfold parameter passed to wbw.kCV
npMSL <- function(x, mu0, blockid = 1:ncol(x),
                      bw=bw.nrd0(as.vector(as.matrix(x))), samebw = TRUE,
                      bwmethod = "S", h=bw, eps=1e-8,
                      maxiter=500, bwiter = maxiter, nbfold = NULL,
                      ngrid=200, post=NULL, verb = TRUE){
  bw <- h # h is alternative bandwidth argument, for backward compatibility
  x <- as.matrix(x)
  n <- nrow(x)      # number of subjects
  r <- ncol(x)      # number of measurements in each subject
  u <- match(blockid, unique(blockid))
  B <- max(u)   	# nb of blocks
  BlS <- rep(0,B)	# block sizes = C_ell in JCGS paper
  for (ell in 1:B) {
    BlS[ell] <- sum(u == ell)}      
  
  if (is.matrix(mu0)) # mu0=centers
    m <- dim(mu0)[1]  else m <- mu0  # mu0 = number of clusters
  
  if(!samebw && !is.matrix(bw)) { # create initial bandwidth matrix h_lj
    bw <- matrix(bw, nrow=max(u), ncol=m)}
  z.hat <- matrix(0, nrow = n, ncol = m)
  tt0 <-  proc.time() # for total time  
  
  ## Initial Values
  if(m == 1) z.hat <- matrix(1, nrow = n, ncol = m) 
  else if(is.null(post)) {
    kmeans <- kmeans(x, mu0)
    for(j in 1:m)
      z.hat[kmeans$cluster==j, j] <- 1
  } 
  else { 
    z.hat <- post ## Currently no error-checking is done here
  }
  
  if (is.null(nbfold) && bwmethod == "CV") {nbfold <- 5} # default value for nbfold-CV
  
  iter <- 0
  finished <- FALSE
  lambda <- matrix(0, nrow = maxiter, ncol = m)
  loglik <- NULL
  loglikseq <- rep(NA,maxiter)
  total_udfl <- 0; total_nan <- 0 # eventual NaN and underflow in C code
  
  tmp <- 1:n
  xtra <- (max(x)-min(x))/10
  grid <- seq(min(x)-xtra, max(x)+xtra, length=ngrid)
  # f stored on a ngrid by m by B array 
  # f_{g,j,ell} = f_{j ell}(u_g)
  # f <- array(1/m/diff(grid[1:2]), c(ngrid, m, B)) 
  # this f was not normalized for being uniform over grid 
  Delta <- diff(grid[1:2])  
  f <- array(1/((ngrid-1)*Delta), c(ngrid, m, B)) 
  oldloglik <- -Inf
  
  orderx <- xx <- list() # preparation for Silverman adaptive bandwidth
  for(k in 1:B) {
    xx[[k]] <- as.vector(x[, u==k]) # data pooled from kth block
    if (!samebw && bwmethod == "S") {
      orderx[[k]] = order(xx[[k]]) # only needed for IQR calculation for bw
    }
  }
  CftEstep <- ifelse(samebw, "npMSL_Estep", "npMSL_Estep_bw")
  # CftEstep <- "npMSL_Estep_bw" # temporary, for testing only the M-step  
  CftMstep <- ifelse(samebw, "npMSL_Mstep", "npMSL_Mstep_bw")
  
  while (!finished) {
    iter <- iter + 1
    bw.old <- bw	# is this needed?
    t0 <- proc.time()
    nb_udfl=0;  # nb of underflows log(0) cancelled in nems_Estep.c
    nb_nan=0;  # nb of nonzero K()*log(0) cancelled in nems_Estep.c
    
    ## Note:  Enter loop assuming E-step is done i.e., z.hat in place    
    ## M-Step
    lambda[iter, ] <- colMeans(z.hat)
    
    ## density estimation in M-step: WKDE-step  
    cs <- colSums(z.hat)
    z.tmp <- sweep(z.hat, 2, cs, "/")
    z.tmp[, cs==0] <- 1/NROW(z.tmp) # Just in case
    
    
    
    ## adaptive bandwidth update - depends on bwmethod and adptbw stage duration in this version
    if (!samebw && iter <= bwiter) {
      # adaptive Silverman's rule
      if (bwmethod == "S") {
        for (ell in 1:B) { # for each block
          r2 <- BlS[ell]   # block size = nb of coordinates
          wts <- apply(z.tmp, 2, function(z) rep(z/r2, r2))
          variances <- colSums(wts*outer(xx[[ell]], colSums(wts*xx[[ell]]),'-')^2)
          iqr <- apply(as.matrix(wts[orderx[[ell]],]), 2, wIQR, 
                       xx[[ell]][orderx[[ell]]],
                       already.sorted=TRUE, already.normalized=TRUE)
          h <- bw[ell, ] <- 0.9 * pmin(sqrt(variances), iqr/1.34) * 
            pmax(1,r2*n*lambda[iter, ])^(-1/5) 
          # Note:  Doesn't allow "sample size" < 1.  
        }
      }
      
      # adaptive nbfold CV computation, nbfold=5 by default
      if (bwmethod == "CV") { 
        # k-fold CV      
        for (ell in 1:B) { # for each block
          r2 <- BlS[ell]   # block size = nb of coordinates
          wts <- apply(z.tmp, 2, function(z) rep(z, r2)) # replicate weights
          for (j in 1:m) bw[ell,j] <- wbw.kCV(xx[[ell]], nbfold = nbfold, w = wts[,j])
        }
      } # end of CV version  
    } # end of bw adaptive stage
    
    
    
    
    z=.C(CftMstep, as.integer(ngrid), as.integer(n),
         as.integer(m), as.integer(r), 
         as.integer(B), as.integer(BlS), as.integer(u),
         as.double(bw), as.double(x), as.double(grid), 
         new.f=as.double(f), 
         as.double(lambda[iter,]), 
         as.double(z.hat))
    
    f <- array(z$new.f, c(ngrid, m, B)) # check sum(f == 0)
    # print(apply(f,2:3,sum) * Delta)
    # print(max(abs(f-f2)))
    #    browser()
    
    ## E-step (for next iteration)
    z=.C(CftEstep, as.integer(ngrid), as.integer(n),
         as.integer(m), as.integer(r), 
         as.integer(B), as.integer(u),
         as.double(bw),
         as.double(x), as.double(grid), f=as.double(f),
         as.double(lambda[iter,]), post=as.double(z.hat),
         loglik = double(1),
         nb_udfl = as.integer(nb_udfl), nb_nan = as.integer(nb_nan))
    
    nb_udfl = z$nb_udfl; nb_nan = z$nb_nan; 
    total_udfl <- total_udfl + nb_udfl
    total_nan <- total_nan + nb_nan
    
    z.hat <- matrix(z$post, n, m)
    if (sum(is.nan(z.hat)) > 0) cat("Error!! NaN in z.hat") # obsolete ?
    loglik <- loglikseq[iter] <- z$loglik
    loglikchange <- loglik - oldloglik
    oldloglik <- loglik
    finished <- iter >= maxiter
    if (iter>1 && max(abs(lambda[iter, ] - lambda[iter-1, ])) < eps)
      finished <- TRUE
    if (verb) {
      t1 <- proc.time()
      cat("iteration", iter, ": lambda ", round(lambda[iter, ], 4))
      cat(" obj change", round(loglikchange,4))
      cat(" time", (t1 - t0)[3])
      if ((nb_udfl > 0) || (nb_nan >0)) cat("\n  ")
      if (nb_udfl > 0) cat("average underflows=", round(nb_udfl/(n*m*r),3)," ")
      if (nb_nan >0) cat("average NaNs=", round(nb_nan/(n*m*r),3))
      # Note: average nb of nan over ngrid convolution
      cat("\n")
    }
  }
  # f <- array(z$f, c(ngrid, m, r)) # obsolete in block version
  
  if (!samebw) {
    rownames(bw) <- paste("block", 1:max(u))
    colnames(bw) <- paste("component", 1:m)
  }
  
  if (verb) {
    tt1 <- proc.time()
    cat("lambda ", round(lambda[iter, ], 4))
    cat(", total time", (tt1 - tt0)[3], "s\n")
  }
  return(structure(list(data = x, posteriors = z.hat, 
                        lambda = lambda[1:iter,], bandwidth = bw, 
                        blockid = u, lambdahat = lambda[iter,], f=f,
                        grid = grid, loglik = loglikseq[1:iter],
                        meanUdfl = total_udfl/(n*m*r*iter),# average underflow
                        meanNaN = total_nan/(n*m*r*iter)), # average NaN's
                   class="npEM")) # should we define a "npMSL" class ?
}



