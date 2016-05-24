#PostProbPois <- function(x,eta,support,mix.prop) {
#      n <- length(x)
#      T <- length(support)
#            
#      z <- .C("postprobpoisson", as.double(x), as.double(eta), 
#                  as.double(support), as.double(mix.prop),
#                  as.integer(n), as.integer(T), double(T), 
#                  post = double(n*T), loglik = double(1), PACKAGE = "rvalues")
#      ans <- list()
#      ans$postprobs <- matrix(z$post, nrow=n)
#      ans$loglik <- z$loglik
#      return(ans)
#}

PostProbPois <- function(x, eta, support, mix.prop) {
    ### Amat is an n.support x n matrix
    Amat <- exp(outer( log(support), x) - outer(support, eta) )
    B <- mix.prop*Amat
    lik <- colSums(B)  ### note this is only proportional to the likelihood
    PP <- t(B)/lik
    ans <- list()
    ans$loglik <- sum(log(lik))
    ans$postprobs <- PP
    return(ans)
}

NPestPoisson <- function(x,eta,maxiter,tol,nmix)  {
  ### This function takes an initial estimate of the mixing distribution
  ### and the loc.scale matrix and returns a final estimate of the 
  ### mixing distribution using the EM algorithm.
  ### Also, it returns an indicator of convergence.
  ### 0 - did converge.  1 - did not converge.
  
  ### densmat should have dimensions [n,T]
  ### densmat[i,j] - p(X_{j}|\theta_{i})
  
  n <- length(x)
  grid.to <- max(x/eta)
  grid.from <- min(1e-4,grid.to/2)
  if(is.null(nmix)) {
      if(n <= 200) {
           nmix <- n
      }
      else{
         ## Number of mixture components grows according to
         ## f(x) = 200 + 5*log(x)*[(x - 200)^(.3)]
  
         nn <- n - 200
         nmix <- 200 + ceiling(5*exp(.3*log(nn))*log(nn))
      }
  }
  log.lik <- rep(0,maxiter + 1)
  counter <- 0
  
  support <- c(1:nmix)*(grid.to/nmix)
  mix.prop <- rep(1/nmix,nmix)

  tmp <- PostProbPois(x,eta,support,mix.prop)
  PP <- tmp$postprobs
  log.lik[1] <- tmp$loglik
  done <- FALSE
  for(k in 1:maxiter)  {
      mix.prop <- colMeans(PP)
      support <- as.vector(crossprod(PP,x)/crossprod(PP,eta))
   
      tmp <- PostProbPois(x,eta,support,mix.prop)
      PP <- tmp$postprobs
      log.lik[k+1] <- tmp$loglik
      done <- (abs((log.lik[k+1] - log.lik[k])/log.lik[k]) < tol) 
      counter <- counter + 1
      if(done) {
          break
      }
  }
  post.mean <- PP%*%support
  log.lik <- log.lik[1:(counter+1)]
  conv <- ifelse(maxiter == counter,1,0)
  return(list(mix.prop=mix.prop,support=support,convergence = conv,log.lik=log.lik,numiter=counter, post.mean=post.mean))
}

