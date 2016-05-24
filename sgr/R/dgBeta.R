#rm(list=ls())

dgBeta <-
function(x,a=min(x),b=max(x),gam=1,del=1) {
  e1 <- 1/(beta(gam,del)*(b-a)^(gam+del-1))
  d <- e1*(x-a)^(gam-1)*(b-x)^(del-1)
  return(d)  
}

## example
#x <- 1:7
#GA <- c(1,3,1.5,8); DE <- c(1,3,4,2.5)
#par(mfrow=c(2,2))
#for (j in 1:4) {
#  plot(x,dgBeta(x,gam=GA[j],del=DE[j]),type="h",
#       panel.first=points(x,gBeta(x,gam=GA[j],del=DE[j]),pch=19),
#       main=paste("gamma=",GA[j]," delta=",DE[j],sep=""),ylim=c(0,.6),
#       ylab="gBeta(x)")  
#}

