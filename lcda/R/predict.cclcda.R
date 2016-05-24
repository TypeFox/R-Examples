###########################################################################
####                       predict.cclcda                              ####
####  ===============================================================  ####
####  - computation of posterior probabilities                         ####
####  - determination class membership by Bayes rule                   ####
###########################################################################


predict.cclcda <- function(object, newdata, ...)
{

  # object must be of class 'lcda'
  if (!inherits(object, "cclcda")) 
     stop("object not of class", " 'cclcda'")

  # results of the lcda procedure
  lca.w <- object$lca.w
  lca.theta <- object$lca.theta
  lca.wmk <- object$lca.wmk
  m <- object$m
  r <- object$r
  d <- object$d
  k <- object$k
  prior <- object$prior
  n <- ncol(newdata)

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

  # convert lca.theta to a matrix
  theta <- unlist(lca.theta, use.names=FALSE)
  theta <- matrix(theta, ncol=m, byrow=TRUE)

  # function for computation of the posterior probabilites
  posterior <- function(z)
  {
  temp <- prior * (lca.wmk%*%apply(theta^x.bin(z), 2, prod))
  return(temp/sum(temp))
  }

  # computation of the posterior probabilites and the Bayes decision
  post <- t(apply(newdata, 1, posterior))
  class <- factor(apply(post, 1, which.max))

  # output of the result
  result <- list(class=class, posterior=post)
  return(result)
}
