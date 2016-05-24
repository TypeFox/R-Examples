samplesize.r <- function(delta, power=.80, sig.level=.05){
  temp <- function(n){
    	tc <- qt(p=sig.level/2,df=n-2,lower.tail=FALSE)
    	rc <- sqrt((tc^2)/((tc^2)+(n-2)))
    	fisher.z <- function(r){(1/2)*log((1+r)/(1-r))}
    	Zl <- (-fisher.z(rc)-fisher.z(delta))/(1/sqrt(n-3))
    	Zu <- (fisher.z(rc)-fisher.z(delta))/(1/sqrt(n-3))
    	return(pnorm(Zl,lower.tail=TRUE) + pnorm(Zu,lower.tail=FALSE)-power)
  }
  ceiling(uniroot(temp,c(4,10e10))[[1]])
}
