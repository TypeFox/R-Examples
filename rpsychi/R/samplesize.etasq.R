samplesize.etasq <- function(k, delta, power=.80, sig.level=.05){
  delta <- sqrt(delta/(1-delta))
  u <- k - 1
  
	temp <- function(n){
  fc <- qf(p=sig.level,df1=u,df2=(u+1)*(n-1),lower.tail=FALSE)
	lamda <- (delta^2)*(n*(u+1))
	v <- (u+1)*(n-1)

	z1b <- (sqrt(2*(u+lamda)-((u+2*lamda)/(u+lamda)))-
  	sqrt((2*v-1)*((u*fc)/v)))/
  	sqrt(((u*fc)/v)+((u+2*lamda)/(u+lamda)))
	return(pnorm(z1b)-power)
	}
	ceiling(uniroot(temp,c(4,10e10))[[1]])
}
