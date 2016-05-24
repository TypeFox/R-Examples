###########################################################################
####                         predict.lcda                              ####
####  ===============================================================  ####
####  - computation of posterior probabilities                         ####
####  - determination class membership by Bayes rule                   ####
###########################################################################



predict.lcda <- function(object, newdata, ...)
{

  # object must be of class 'lcda'
  if (!inherits(object, "lcda")) 
     stop("object not of class", " 'lcda'")

  # results of the lcda procedure
  lca.w <- object$lca.w
  lca.theta <- object$lca.theta
  m <- object$m
  r <- object$r
  d <- object$d
  k <- object$k
  prior <- object$prior

  # subprocedure to convert observations to binary vectors with one bit per outcome per variable
  x.bin <- function(x=x)
  {
  res <- numeric(sum(r))
  for (i in 1:d)
  {
  for (j in 1:r[i])
  {
  if((is.na(x[i])==FALSE & x[i]==j))
  res[(sum(r[1:i]))-(r[i]-j)] <- 1
  }
  }
  return(res)
  }


  # convert lca.theta to a list of k vectors: all parameters of each latent class consecutively
  theta <- lapply(lca.theta, function(z) as.numeric(t(matrix(as.numeric(unlist(z)), nrow=m))))

  # create a list containing theta and m
  wtheta <- list()
  for (i in 1:k)
  {
  wtheta[[i]] <- list()
  wtheta[[i]][[1]] <- theta[[i]]
  wtheta[[i]][[2]] <- m[i]
  }

  # function to compute the posterior probabilites
  posterior <- function(x)
  {
  temp <- lapply(wtheta, function(z) z[[1]]^rep(x.bin(x),z[[2]]))
  temp2 <- lapply(temp, function(z) apply(matrix(z, ncol=sum(r), byrow=TRUE), 1,prod))
  tempsum <- numeric(k)
  for (i in 1:k){tempsum[i] <- (prior[i]*(sum(temp2[[i]] * lca.w[[i]])))}
  res <- tempsum/sum(tempsum)
  return(res)
  }

  # computation of the posterior probabilites and the Bayes decision
  post <- t(apply(newdata, 1, posterior))
  class <- factor(as.numeric(apply(post, 1, which.max)))

  # output of the result
  result <- list(class=class, posterior=post)
  return(result)
}
