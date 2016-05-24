###########################################################################
####                       predict.cclcda2                             ####
####  ===============================================================  ####
####  - computation of posterior probabilities                         ####
####  - determination class membership by Bayes rule                   ####
########################################################################### 

predict.cclcda2 <- function(
                                   object,              # object of class cclcda2
                                   newdata,             # new data
                                   ...
                                  )
{

  # object must be of class 'lcda'
  if (!inherits(object, "cclcda2")) 
     stop("object not of class", " 'cclcda2'")


  x <- newdata
  # results of the cclcda procedure
  lca.w <- object$lca.w
  lca.theta <- object$lca.theta
  lca.wmk <- object$lca.wmk
  m <- object$m
  r <- object$r
  d <- object$d
  k <- object$k
  prior <- object$prior
  n <- ncol(newdata)

  theta <- unlist(lca.theta, use.names=FALSE)
  theta <- matrix(theta, ncol=m, byrow=TRUE)



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

  posterior <- function(z)
  {
  temp <- prior*apply(t(t(lca.wmk)*apply(theta^x.bin(z), 2, prod)),1,sum)
  return(temp/sum(temp))
  }

  post <- t(apply(x,1,posterior))
  class <- factor(apply(post, 1, which.max))

  # output of the result
  result <- list(class=class, posterior=post)
  return(result)

}

