my1.bug<-function(){
	#likelihood: joint distribution of ys
	for(t in 1:n){
		y[t]~dnorm(0,yisigma2[t]);	
		yisigma2[t]<-1/exp(theta[t]);
	}

	#prior distributions
	mu~dnorm(0,0.1);
	phistar~dbeta(20,1.5);
	itau2~dgamma(2.5,0.025);
	beta<-exp(mu/2);
	phi<-2*phistar-1;
	tau<-sqrt(1/itau2);

	theta0 ~ dnorm(mu,itau2);
	thmean[1] <- mu + phi*(theta0-mu);
	theta[1] ~ dnorm(thmean[1],itau2);
	for (t in 2:n) {
		thmean[t] <- mu + phi*(theta[t-1]-mu);
		theta[t] ~ dnorm(thmean[t],itau2);
	}
}
