pHayStonLSA<-function(h,k,delta=.001){

	###Uses Hayter-Liu 1996
	our.grid<-seq(-8,8,delta)
	len<-length(our.grid)

	init.grid<-pnorm(our.grid+h)

	if(k>2){
	for(i in 3:k){
		new.grid<-cumsum(init.grid*dnorm(our.grid)*delta)+init.grid*(pnorm(our.grid+h)-pnorm(our.grid))
		init.grid<-new.grid
	}
	}
	1-sum(dnorm(our.grid)*init.grid*delta)

}