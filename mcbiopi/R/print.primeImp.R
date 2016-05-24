print.primeImp<-function(x,...){
	cat("Prime Implicants:\n")
	for(i in 1:length(x$vec.primes))
		cat("   ",x$vec.primes[i],"\n")
}

