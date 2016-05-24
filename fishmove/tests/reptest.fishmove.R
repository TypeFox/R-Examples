#
# Test how many replicates are nessary for fishmove to get kind of stable results
# 
# Author: Johannes Radinger
###############################################################################


#library(fishmove)
#large test
#reptest.fishmove(a=25,b=c(5,10,25,50,75,100,200,250,500,1000))

#small test
#reptest.fishmove(a=10,b=c(10,50,100,500))

reptest.fishmove <- function(a,b){
	
	a<-a # number of repeats of repeats (n for boxplot per group)
	b<-b #vector of number of repeats
	c <- rep(b,rep(a,length(b)))
	
	# Process calculates mean of sigma_stat and sigma_mob for different sets of repeats and writes results into ar
	ar <- array(NA, dim = c(length(c), 3))
	for (i in 1:length(c)) {
		f <- fishmove(species="Salmo trutta fario",rep=c[i])
		ar[i,1:2] <- f$pred.fishmove["fit", , , , ,]
		ar[i,3] <- c[i]
	}
	
	# plotting results to see variation in the mean (should get lesser with increasing number of repeats)
	boxstat <- boxplot(ar[,1]~ar[,3],main="Number of fishmove repeats for stat", 
			xlab="Number of repeats", ylab="mean of sigma_sat")
	boxmob <- boxplot(ar[,2]~ar[,3],main="Number of fishmove repeats for mob", 
			xlab="Number of repeats", ylab="mean of sigma_mob")
	
	par(mfrow=c(2,1))
	boxstat
	boxmob	
}

