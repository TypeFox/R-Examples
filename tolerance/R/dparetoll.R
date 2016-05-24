dpareto.ll <- function(x, theta = NULL, ...){
  th.table <- cbind(c(0,0.050,0.089,0.126,0.164,0.203,0.244,0.286,0.332,0.380,
                      0.431,0.486,0.545,0.609,0.678,0.753,0.835,0.925,1.025,1.135,
                      1.258,1.395,1.549,1.723,1.921,2.146,2.404,2.701,3.044,3.442,
                      3.904,4.442,5.071,5.807,6.669,7.682,8.870),seq(0,.36,by=.01))
  inf.test <- which(x==Inf)
  if(length(inf.test)>0){
    x[inf.test] <- max(x[-inf.test])
    warning("Values of x equal to 'Inf' are set to the maximum finite value.",call.=FALSE)
  }
  x.bar <- mean(x)
  if(is.null(theta)){
    if(x.bar<=max(th.table[,1])){
      ind <- max(which(th.table[,1]<x.bar))
      theta <- approx(th.table[ind:(ind+1),],xout=x.bar)$y
    } else{
      S.hat <- sapply(1:length(x),function(i) mean(x>=x[i]))
      theta <- prod(S.hat^((sum(log(1+x)))^(-1)))
    }
  }
  tmp <- which(ddpareto(x,theta=theta,log=TRUE)==-Inf)
  if(length(tmp)>0){
    x[tmp] <- max(x[-tmp])
    warning("Numerical overflow problem when calculating log-density of some x values.  The problematic values are set to the maximum finite value calculated.",call.=FALSE)
  }
  ll.f <- function(theta) -sum(ddpareto(x,theta=theta,log=TRUE))
  fit <- try(suppressWarnings(stats4::mle(ll.f,start=list(theta=theta),lower=0,upper=1,method="Brent")),silent=TRUE)
  if(class(fit)=="try-error") stop(paste("Numerical optimization of the MLE failed.  Consider trying a different starting value for theta.","\n"))
  fit
}  
