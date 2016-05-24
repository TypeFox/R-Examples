power.r <- function(n, delta, sig.level = .05){
	tc <- qt(p=sig.level/2,df=n-2,lower.tail=FALSE)
	rc <- sqrt((tc^2)/((tc^2)+(n-2)))
	fisher.z <- function(r){(1/2)*log((1+r)/(1-r))}
	Zl <- (-fisher.z(rc)-fisher.z(delta))/(1/sqrt(n-3))
	Zu <- (fisher.z(rc)-fisher.z(delta))/(1/sqrt(n-3))
	pnorm(Zl,lower.tail=TRUE) + pnorm(Zu,lower.tail=FALSE)
}