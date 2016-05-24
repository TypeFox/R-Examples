weighted.entropy <-
function(x,w=NULL,alpha=1){
  if(is.null(w)) w<-rep(1/length(x),length(x))
  x <- as.numeric(x)
  n <- length(x)
  fx<- w/sum(w)
  if (is.null(alpha)) 
    alpha <- 1
  if (alpha == 0) { # yields mean logarithmic 
    entropy <- -sum(log(x/weighted.mean(x,w))*fx) 
  }
  else if (alpha == 1) 
    entropy <- sum(x/weighted.mean(x,w)*log(x/weighted.mean(x,w))*fx)
  else {
    entropy <- 1/(alpha * (alpha - 1))*t((x/weighted.mean(x,w))^alpha-1)%*%fx
  }
  return(entropy)
}
