weighted.atkinson <-
function (x, w=NULL, epsilon = 1,wscale=1000) {
  if(is.null(w)) w<-rep(1/length(x),length(x))
  x <- as.numeric(x)
  wr<- round(w*wscale,digits=0)
  xw<- rep(x,wr)
  if (is.null(epsilon)) 
    epsilon <- 1
  if (epsilon == 1) 
    A <- 1 - (exp(mean(log(xw)))/mean(xw))
  else {
    xw <- (xw/mean(xw))^(1 - epsilon)
    A <- 1 - mean(xw)^(1/(1 - epsilon))
  }
  return(A)
}
