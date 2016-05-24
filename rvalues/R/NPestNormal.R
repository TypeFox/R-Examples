#PostProbNorm <- function(x,std_err,support,mix.prop) {
#      n <- length(x)
#      T <- length(support)
#     
#      z <- .C("postprobnormal", as.double(x), as.double(std_err), 
#                  as.double(support), as.double(mix.prop),
#                  as.integer(n), as.integer(T), double(T), 
#                  post = double(n*T), loglik = double(1), PACKAGE = "rvalues")
#      ans <- list()
#      ans$postprobs <- matrix(z$post, nrow=n)
#      ans$loglik <- z$loglik
#      return(ans)
#}


PostProbNorm <- function(x, std_err, support, mix.prop) {
    ### Amat is nsupport x n matrix
    Amat <- t(dnorm(outer(x, support, FUN="-"), sd=std_err))
  
    B <- mix.prop*Amat
    lik <- colSums(B)
    PP <- t(B)/lik
    ans <- list()
    ans$loglik <- sum(log(lik))
    ans$postprobs <- PP
    return(ans)
}

NPestNormal = function(x,std_err,maxiter,tol,nmix)  {
  ### This function takes an initial estimate of the mixing distribution
  ### and the loc.scale matrix and returns a final estimate of the 
  ### mixing distribution using the EM algorithm.
  ### Also, it returns an indicator of convergence.
  ### 0 - did converge.  1 - did not converge.
  
  ### densmat should have dimensions [n,T]
  ### densmat[i,j] - p(X_{j}|\theta_{i})
  
  n <- length(x)
  grid.to <- max(x)
  grid.from <- min(x)
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
      #num.mixcomponents <- round(n^(2/3),0) 
  }
  log.lik <- rep(0,maxiter + 1)
  counter <- 0
  
  support <- c(1:nmix)*(grid.to/nmix)
  mix.prop <- rep(1/nmix,nmix)
  ss <- 1/(std_err^2)

  tmp <- PostProbNorm(x, std_err, support, mix.prop)
  PP <- tmp$postprobs
  log.lik[1] <- tmp$loglik
  done <- FALSE
  for(k in 1:maxiter)  {
      mix.prop <- colMeans(PP)
      support <- as.vector(crossprod(PP,x*ss)/crossprod(PP,ss))
      
      tmp <- PostProbNorm(x, std_err, support, mix.prop)
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
  conv = ifelse(maxiter == counter,1,0)
  return(list(mix.prop=mix.prop,support=support,convergence = conv,log.lik=log.lik,numiter=counter,post.mean=post.mean))
}
