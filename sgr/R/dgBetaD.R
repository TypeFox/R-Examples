dgBetaD <-
function(x,a=min(x),b=max(x),gam=1,del=1,ct=1) {  
  X <- a:b
  t1 <- (beta(gam,del))^(-1)
  t2 <- ((b-a+2*ct)^(gam+del-1))^(-1)
  t3 <- (X-a+ct)^(gam-1)
  t4 <- (b-X+ct)^(del-1)
  l <- t1*t2*t3*t4
  L <- sum(l)
  
  t3 <- (x-a+ct)^(gam-1)
  t4 <- (b-x+ct)^(del-1)
  l <- t1*t2*t3*t4
  L <- l/L
  
  y <- ifelse((x %in% X),L,0)
  return(y)
}

## example
#x <- 1:7
#GA <- c(1,3,1.5,8); DE <- c(1,3,4,2.5)
#par(mfrow=c(2,2))
#for (j in 1:4) {
#  plot(x,gBetaD(x,gam=GA[j],del=DE[j]),type="h",
#       panel.first=points(x,gBetaD(x,gam=GA[j],del=DE[j]),pch=19),
#       main=paste("gamma=",GA[j]," delta=",DE[j],sep=""),ylim=c(0,.6),
#       ylab="gBeta(x)")  
#}