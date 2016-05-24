BayesianMass <- function(lambda,r,n){
	# n = size of discernment space.
	# 
	if(n<=1){ 
		stop('error in Bayesian Mass\n')
	}else{
		Card=c(0,1);
		for(i in 1:(n-1)){
			Card=c(Card,(Card+1));
		}
	}

	Card[1]=1; # for division

	mt=lambda*(1/Card)^r;
	mt[1]=0;
	md=mt/sum(mt);
    return(md)
}
