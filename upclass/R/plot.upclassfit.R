.packageName <- 'upclass'

plot.upclassfit <- function (x, ...) 
{
  if (attr(x, "class") != "upclassfit") 
    stop("Incorrect class")
  z<-x$Best$test$z
  G<-x$Best$G

 N<-nrow(z) 
  frame()
  par(usr=c(-0.01*N,N*(1.25),-0.04,1.04))
  axis(2)
  axis(1,at=pretty(1:N,10))
  matplot(z,add=TRUE,pch=1:G,col=1:G)
  legtext<-paste("Group",1:G)
  legend(N*1.02,1,legend=legtext,col=1:G,pch=1:G)
  
  title(main="Posterior Probability of Group Membership",ylab="Posterior Probability",xlab="Observation")
}
