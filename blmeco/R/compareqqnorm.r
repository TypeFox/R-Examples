compareqqnorm <- function(mod){

N <- length(resid(mod))    # sample size
sigma <- summary(mod)$sigma
par(mfrow=c(3,3)) 		# divide the graphic window in 9 sub-windows
rnum<-sample(1:9, 1) 		# draw a random number between 1 and 9
for(i in 1:(rnum-1)){
  x<-rnorm(N, 0, sigma)
  qqnorm(x, main=i)
  qqline(x)
}
qqnorm(resid(mod), main=rnum)
qqline(resid(mod))
for(i in (rnum+1):9){
  x<-rnorm(N, 0, sigma)
  qqnorm(x, main=i)
  qqline(x)
}
return(rnum)
}
