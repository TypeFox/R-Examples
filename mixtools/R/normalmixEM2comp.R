# fast function for two-component univariate normal mixture

normalmixEM2comp <- function(x, lambda, mu, sigsqrd, eps=1e-8, maxit=1000, verb=FALSE) {
  arbvar <- (length(sigsqrd)==2)
  mu1 <- mu[1]; mu2 <- mu[2]
  sigsqrd1 <- sigsqrd[1]; sigsqrd2 <- sigsqrd[arbvar+1]
  mx <- mean(x)
  const <- length(x) * 0.918938533204673 # i.e., times log(2*pi)/2
  dl <- 1 + eps
  iter<-0
  ll <- rep(0, maxit+1)
  a1<-(x-mu1)^2; b1<-(lambda/sqrt(sigsqrd1))*exp(-a1/2/sigsqrd1)
  a2<-(x-mu2)^2; b2<-((1-lambda)/sqrt(sigsqrd2))*exp(-a2/2/sigsqrd2)
  l <- sum(log(b1+b2))

  while (dl>eps && iter<maxit) {
    iter<-iter+1
    ll[iter] <- l
    postprobs <- b1/(b1+b2)
    lambda<-mean(postprobs)
    mu1<-mean(postprobs*x)/lambda
    mu2<-(mx-lambda*mu1)/(1-lambda)
    if (arbvar)  {   
      sigsqrd1<-mean(postprobs*a1)/lambda
      sigsqrd2<-mean((1-postprobs)*a2)/(1-lambda)
    } else {
      sigsqrd1 <- sigsqrd2 <- mean(postprobs*a1 + (1-postprobs)*a2) 
    }
    a1<-(x-mu1)^2; b1<-(lambda/sqrt(sigsqrd1))*exp(-a1/2/sigsqrd1)
    a2<-(x-mu2)^2; b2<-((1-lambda)/sqrt(sigsqrd2))*exp(-a2/2/sigsqrd2)

    oldl<-l    
    l <- sum(log(b1+b2))
    dl<-l-oldl
    if (verb) {
      cat("iteration =", iter, " log-lik diff =", dl, " log-lik =", 
          l-const, "\n")
    }
  }
  cat("number of iterations=", iter, "\n")
  iter <- iter+1
  ll[iter] <- l
  postprobs <- cbind(postprobs, 1-postprobs)
  colnames(postprobs) <- c(paste("comp", ".", 1:2, sep = ""))
  out <- list(x=x, lambda = c(lambda,1-lambda), mu = c(mu1, mu2), 
       sigma = sqrt(c(sigsqrd1, sigsqrd2)[1:(1+arbvar)]), 
       loglik = l - const, posterior = postprobs, 
       all.loglik=ll[1:iter] - const, 
       restarts=0, ft="normalmixEM")
  class(out) <- "mixEM"
  out
}

