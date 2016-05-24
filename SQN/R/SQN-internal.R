SQN.demo.data <- function(){
mu0a=mu0b=rnorm(1000,5,1.8) ## control probes 
  mu1a=rnorm(5000,5,0.8);mu1b=mu1a+2  ## probes for signal
  ya=c(rnorm(1000,mu0a,.3),rnorm(5000,mu1a,.3))
  yb=c(rnorm(1000,mu0b,.3),rnorm(5000,mu1b,.3))+1 ## a systematic bias introduced
  Y=cbind(ya,yb)   #before normalization
  Ynorm=SQN(Y,ctrl.id=1:1000)  #after normalization
}
##  par(mfrow=c(1,2))
##  boxplot(Y,main="before normalization")
##  boxplot(Y[1:1000,],add=T,col=3,boxwex=.4)

##  boxplot(Ynorm,main="after normalization")
##  boxplot(Ynorm[1:1000,],add=T,col=3,boxwex=.4)
##  legend(1,10,legend=c("probes for signal","negative control probes"),text.col=c(1,3),bg="white")
