# weather generator daily

poisson.rain <- function(rate=1, ndays=30, shape=1, scale=1,plot.out=T){
 nk <- rpois(1,lambda=rate*ndays)
 tau1 <- array(); ctau1 <- array()
 for(i in 1:nk){
 tau1[i] <- ceiling(rexp(1,rate))
 ctau1[i] <- sum(tau1[1:i])
}
if (ctau1[nk]<=ndays) k <- nk else
k <- min(which(ctau1>ndays))-1
tau <- tau1[1:k]; ctau <- ctau1[1:k]

# amount of rain
x <- round(rweibull(k, shape,scale),2)
y <- cbind(tau,ctau,x)
days <- seq(1,ndays); rain <- days; rain[] <- 0
for(i in 2:ndays){
 ab <- which(ctau==days[i])
 if(length(ab)>0) rain[i] <- x[ab]
}
z <- cbind(days,rain)
 if(plot.out == T){
 panel4(size=7)
 plot(days,rain,type="s",xlab="Day",ylab="Rain (cm/day)")
 hist(rain,prob=T,main="",xlab="Rain (cm/day)")
 hist(tau,prob=T,main="",xlab="Interarrival time")
}
return(list(y=y,z=z))
}


