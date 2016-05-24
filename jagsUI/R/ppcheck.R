
pp.check <- function(x, actual, new){
  if(class(x)!="jagsUI"){stop('Requires jagsUI object as input')}
  devAskNewPage(ask=FALSE)
  actual <- eval(parse(text=paste('x$sims.list$',actual,sep="")))
  new <- eval(parse(text=paste('x$sims.list$',new,sep="")))
  
  bpval <- mean(actual>new)
  
  plot(x = actual, y = new, xlab="Actual Dataset", ylab="Simulated Dataset", 
       main = paste('Posterior Predictive Check','\n','Bayesian P-value =',round(bpval,2)))
  abline(1,1)
  

}