### multivariate npEM with multivariate conditionnally independent blocks
### functions for mvnpEM with mvwkde code all in C

######################################################################
## mvnpEM final (late 2015) version calling C codes for mvwkde's
## both for same and adaptive bandwidths
## adding a default bandwidth matrix passed as an argument
## bwdefault: a r-vector of fixed bandwidths per coordinates
##  default to Silverman's rule per coordinates, only when samebw=TRUE
## bw returns the standard deviation of the kernel density
## init= for initialization option (June 2015)
######################################################################
mvnpEM <- function (x, mu0, blockid = 1:ncol(x), samebw = TRUE,
                    bwdefault = apply(x,2,bw.nrd0), 
                    init=NULL,
                    eps=1e-8, maxiter = 500, verb = TRUE){
  x = as.matrix(x)
  n = nrow(x); r = ncol(x)
  bk = blockid #b_k=l indicate the kth coordinate belongs to lth block.
  B = max(bk) # total number of blocks 
  # ! Caution if blocks are numbered not consecutively, use unique?
  loglik = NULL
  tt0 <-  proc.time() # for total time
  # coordinate of X in lth block (block of dl-variate densities)
  # dl[l] = number of coordinates in l-th block
  # xx[[j]] = (n,dl[l]) -matrix of block l observations 
  dl=c(); for (l in 1:B){dl[l] = sum(bk==l)}
  xx = list()
  for (l in 1:B){
    xx[[l]] = as.matrix(x[,bk==l]) # coordinate of X in lth block
  }
  if (is.matrix(mu0))     m <- dim(mu0)[1]  # mu0 = initial means
  else     m <- mu0  # when mu0 = number of components
  ## NB initial means used only if kmeans used for init
  ## z.hat are posterior zij's, initialized using kmeans
  if (is.null(init)){
    z.hat <- matrix(0, nrow = n, ncol = m) 
    if (m == 1) z.hat <- matrix(1, nrow = n, ncol = m) # ?? m can never be 1
    else{
      kmeans <- kmeans(x, mu0)
      for(j in 1:m)   z.hat[kmeans$cluster==j, j] <- 1}
  }
  
  if (! is.null(init)) {z.hat <- init} # ToDo: check consistency dim (n,m)
  
  
  lambda = matrix(0, maxiter, m) # storing lambda's along iterations
  finished = FALSE 
  bw=matrix(0,m,r) # only used for adaptbw version
  
  if (samebw){ # bw computed once for all iterations and mvwkde calls
    # bw <- apply(x,2,bw.nrd0)
    # bw <- bwdefault^2*diag(r) # diagonal bw matrix appropriate for mvwkde
    bw <- bwdefault*diag(r) # diagonal bw matrix, not variances, appropriate for mvwkde
  }
  ## loop
  iter = 0 #start point
  while (!finished){
    iter = iter + 1 #step "iter + 1"
    #	t0 = proc.time()
    #################################################################
    ## M-step for the Euclidean parameter
    lambda[iter,] = colMeans(z.hat) #\lambda_j^(t+1) = 1/n* \sum_i^n (p_ij)^t
    #################################################################
    ## Weighted Kernel Density function step
    # wts[,j] = normalized weights for component j 
    cs <- colSums(z.hat) # or use n*lambda[iter,] to avoid re summing over n ? 
    wts <- sweep(z.hat, 2, cs, "/")
    wts [, cs==0] <- 1/NROW(wts)
    
    if (samebw){ # default samebw=TRUE	
      lambda.f <- matrix(NA,n,m) # lambda.f[i,j] = lambda_j f_j(x_i)
      fkernel <- matrix(1,n,m)
      # Faster new version: no loop in j, all is done in C block per block
      for (l in 1:B) {
        d = dl[l]
        tmp=as.matrix(bw[bk==l,bk==l]) # updated here, not sqrt anymore
        h=as.vector(diag(tmp));
        ans <- .C("mvwkde_samebw", n = as.integer(n),d = as.integer(d), 
                  m = as.integer(m), h = as.double(h), x=as.double(xx[[l]]), 
                  u=as.double(xx[[l]]), z=as.double(wts), f=double(n*m))
        fl = matrix(ans$f,n,m) # f_jl estimate
        fkernel <- fkernel*fl # product over blocks
      }
      lambda.f <- sweep(fkernel, 2, lambda[iter, ], "*")
    } # end of samebw case
    
    # Adaptive bandwidth case - NOT YET BEST C VERSION (for all j simultaneously ?)
    if (!samebw) {   
      for (k in 1:r) {
        # compute adaptive bw h_jk^t
        # use o <- order(x[,k]) to re-use it twice
        var <- colSums(wts*outer(x[,k], colSums(wts*x[,k]),'-')^2)
        iqr <- apply(as.matrix(wts[order(x[,k]),]),2,wIQR,
                     x[,k][order(x[,k])], 
                     already.sorted=TRUE, already.normalized=TRUE)
        bw[,k] <- 0.9*pmin(sqrt(var), iqr/1.34)*pmax(1,n*lambda[iter, ])^(-1/5)
      }
      lambda.f <- matrix(NA,n,m) # lambda.f[i,j] = lambda_j f_j(x_i)
      fkernel <- matrix(1,n,m)
      #for (j in 1:m){
      #lda.f = lambda[iter, j]
      for (l in 1:B){
        d = dl[l];
        H = as.matrix(bw[, bk==l]);
        ans <- .C("mvwkde_adaptbw", n = as.integer(n),d = as.integer(d), 
                  m = as.integer(m), H = as.double(H), x=as.double(xx[[l]]), 
                  u=as.double(xx[[l]]), z=as.double(wts), f=double(n*m))
        fl = matrix(ans$f,n,m)# f_jl estimate
        fkernel <- fkernel*fl # product over blocks
      }
      lambda.f <- sweep(fkernel, 2, lambda[iter, ], "*")
      #}
    }# end of adaptive case
    
    ################################################################
    ## E-step (for next iteration)
    z.hat = lambda.f/rowSums(lambda.f) #p_ij^t, Z_ij^t, update z.hat
    loglik <- c(loglik,sum(log(rowSums(lambda.f)))) # log-likelihood
    ## End
    finished = iter >= maxiter
    if (iter>1) {maxchange = max(abs(lambda[iter,] - lambda[iter-1,]))
    finished = finished | (maxchange < eps) }   
    if (verb) {
      # 		t1 <- proc.time()
      cat("iteration", iter, ": lambda ", round(lambda[iter, ], 4),"\n")
      #	      cat(" time", (t1 - t0)[3], "\n")   
    }
  }# End while
  ##### Output
  if (verb) {
    tt1 <- proc.time() # total time ending
    cat("# iter", iter)
    cat(", lambda ", round(lambda[iter, ], 4))
    cat(", total time", (tt1 - tt0)[3], "s\n")  }
  # final bandwidth output depending on samebw switch
  if (samebw) bw <- diag(bw)
  return(structure(list(data = x, posteriors = z.hat, lambda = lambda[1:iter,], 
                        blockid = bk, samebw=samebw, bandwidth = bw, lambdahat = lambda[iter,], 
                        loglik = loglik), 
                   class = "mvnpEM"))
} # End function.









#######################################################
# plot marginal (univariate) wkde's  from mvnpEM output
# a plot.mvnpEM method for mvnpEM class
# mu, v: true parameters, for gaussian models only
# ... passed to hist first level plotting
plot.mvnpEM <- function(x, truenorm=FALSE, mu=NULL, v=NULL,
                        lgdcex =1,
                        ...) {
  mix.object <- x
  if (!inherits(mix.object, "mvnpEM"))
  a <- x                      
  x <- a$data;  r <- ncol(x); m <- ncol(a$posteriors)
  rr <- sqrt(r); frr <- floor(rr)
  if (frr == rr) par(mfrow=c(rr,rr)) else {
    if (frr*(frr+1) >= r) 
      par(mfrow=c(frr, frr+1)) else par(mfrow=c(frr+1, frr+1))}
  #	if ((r %% 2) == 0) par(mfrow=c(floor(r/2), floor(r/2))) else {
  #		par(mfrow=c(floor(r/2)+1, floor(r/2)+1))}
  for (k in 1:r) {	# for each coord (whatever its block)
    xx <- x[,k]
    uk <- seq(min(xx),max(xx),len=100)
    tt <- paste("block", a$blockid[k],", coord",k)
    hist(xx, col=8, freq=F, xlab="", main=tt, ...)
    for (j in 1:m) {
      wj <- a$post[,j]   # weights for component j
      if (a$samebw) bw <- a$bandwidth[k] else {
        bw <- a$bandwidth[j,k]}
      f <- wkde(xx, u=uk, w=wj, bw=bw, sym=F)
      lines(uk, a$lambdahat[j]*f, col=j)
      # add true Gaussian marginal if requested/known
      if (truenorm) {
        lines(uk,lambda[j]*dnorm(uk,mean=mu[j,k], sd=sqrt(v[j,k,k])), 
              lty=2,lwd=2, col=j)}
    }
    if (a$samebw) subt <- paste("same bw:",round(a$bandwidth[k],3)) else {
      subt <- "adapt bw: "
      for (j in 1:m) subt <- paste(subt, round(a$bandwidth[j,k],3))
    }
    title(sub=subt)
  }
  lgd <- NULL
  for (j in 1:m) lgd <- c(lgd, paste("comp",j))
  legend("topright", lgd, col=1:m, lty=1, cex=lgdcex)
}





###########################################################
print.mvnpEM <- function(x,...)
{
	n <- NROW(x$data)
	r <- NCOL(x$data)
	m <- length(x$lambdahat)
	cat(paste("Observations:", n, "\n"))
	cat(paste("Coordinates per observation:", r, "\n"))
	cat(paste("Mixture components:", m, "\n"))
	if (r > 1) {
	B <- max(x$blockid)
	cat(paste("Blocks (of conditionally iid coordinates):",B, "\n\n"))
	}
	dp = match(c("data", "posteriors", "lambda", "mu"), names(x), nomatch = 0)
	print.default(structure(x[-dp], class = class(x)), ...)
	invisible(x)
}


print.summary.mvnpEM <-function (x, digits = 3, ...) 
{
  if (x$r > 1) 
    cat(paste(x$n, "observations,", x$r, "coordinates,", 
              x$m, "components, and", x$B, "blocks.\n\n"))
  else cat(paste(x$n, "univariate observations, and", x$m, 
                 "components.\n\n"))
  cat("Means (and std. deviations) for each component:\n")
  for (l in 1:x$B) {
    coords <- 1
    if (x$r > 1) {
      coords <- x$blockid == l
      cat(paste("  Block #", l, ":  Coordinate", sep = ""))
      cat(ifelse(sum(coords) > 1, "s ", " "))
      cat(which(coords))
      cat("\n")
    }
    if (sum(coords)==1) { 
      for (j in 1:x$m){
        cat(paste("   Component", j,": "))
        cat(paste(signif(x$means[j,coords], digits), 
                  " (", signif(x$variances[[j]][coords,coords],digits), ")  ", sep = ""))
      }
      cat("\n")
    }
    else {
      for (j in 1:x$m){
        cat(paste("   Component", j,": "))
        cat(paste(signif(x$means[j,coords ], digits),  sep = ""))
        cat("\n")
        print(signif(x$variances[[j]][coords,coords], digits))
        cat("\n")
      }}
  }
}

summary.mvnpEM <- function(object,...)
{
  n <- NROW(object$data);n
  r <- NCOL(object$data);r
  m <- length(object$lambdahat);m
  B <- max(object$blockid);B
  mu <- matrix(0, m, r);mu
  v <- array(0,c(m,r,r));v[1,,];v[2,,]
  normpost <- sweep(object$post, 2, sums <- colSums(object$post),"/")
  for (l in 1:B){
    coords <- object$blockid == l;coords
    sc <- sum(coords);sc
    xx <- as.matrix(object$data[, coords]);xx
    var<- list();mean<-list()
    for (j in 1:m){
      wts <- normpost[, j]
      M<-mu[j,coords] <- colSums(wts*as.matrix(xx))
      if (sc==1) v[j,coords,coords]<-
        colSums(wts*outer(as.matrix(xx), colSums(wts*as.matrix(xx)),'-')^2)
      else{
        for (t1 in 1:sc){
          for (t2 in 1:sc){
            v[j,coords,coords][t1,t2] <- v[j,coords,coords][t2,t1] <-
              colSums(wts*outer(as.matrix(xx[,t1]), colSums(wts*as.matrix(xx[,t1])),'-')
                      *outer(as.matrix(xx[,t2]), colSums(wts*as.matrix(xx[,t2])),'-'))}
        }
      }
      var[[j]] = v[j,,]
    }
  }
  rownames(mu) <-  paste("component",1:m)
  colnames(mu) <-  paste("coordinate", 1:r)
  names(var) <- paste("component",1:m)
  ans <- list(n = n, m = m, r = r, B = B, blockid = object$blockid,
              means = mu, variances = var)
  class(ans)<-"summary.mvnpEM"	
  ans
}


