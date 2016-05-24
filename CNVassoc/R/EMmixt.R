EMmixt <-
function(x, mu.surrog, sd.surrog, pi, tol=10^-5, max.iter=1000, verbose=FALSE, var.equal=TRUE)
{
  rn <- rownames(x)
  if (length(dim(x))>1) x <- x[,1]
  mean.weight <- function(x,w) sum(x*w)/sum(w)
  sd.weight <- function(x,w) sqrt(sum(w*(x-mean.weight(x,w))^2)/sum(w))
  error <- Inf
  iter <- 1
  I <- NROW(x)
  J <- NROW(mu.surrog)
  while (error > tol & iter < max.iter){
    theta.old<-c(mu.surrog, sd.surrog, pi)
    ff <- sapply(1:J,function(j) pi[j] * dnorm(x, mu.surrog[j], sd.surrog[j]))
    pi.post <- ff/rowSums(ff)
    pi <- colMeans(pi.post)
    mu.surrog <- sapply(1:J,function(j) mean.weight(x,pi.post[,j]))
    sd.surrog <- sapply(1:J,function(j) sd.weight(x,pi.post[,j]))
    if (var.equal == TRUE) sd.surrog<-rep(sqrt(sum(pi*sd.surrog^2)), J)
    theta <- c(mu.surrog, sd.surrog, pi)
    error <- sqrt(sum((theta-theta.old)^2)/length(theta))
    if (verbose) cat("-iter.",iter," error=",error,"\n")
    iter<-iter+1
  }
  if (iter == max.iter) 
    stop("Maximum iteration reached")
  out <- list()
  rownames(pi.post) <- rn
  out$z <- pi.post
  out$mu.surrog <- mu.surrog
  out$sd.surrog <- sd.surrog
  out$pi <- pi
  out$loglike <- sum(log(rowSums(ff)))
  out$numpar <- if (var.equal) J*2+1 else J*3
  out$bic<- -2*out$loglike + out$numpar*log(I)
  out$aic<- -2*out$loglike + 2*out$numpar
  out$x<-x
  out
}



       