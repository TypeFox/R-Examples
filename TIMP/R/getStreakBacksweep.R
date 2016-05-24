"getStreakBacksweep" <- function(streakT, k, mu, x) {
	
	Y <- matrix(0,length(x),length(k))
	for(i in 1:length(k)) {
	      x1<-exp(-k[i]*(x-mu + streakT))
	      x2<-exp(-k[i]*((streakT/2) - (x-mu)))
	      x3 <- exp(-k[i]*streakT)
	      Y[,i] <- (x1+x2)/(1-x3)
	}
	Y
}
