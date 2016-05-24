# Loss aversion based c* calculation routines
# for simulated marginal effects
# based on Esarey and Danneman, "A Quantitative Method for Substantive Robustness Assessment"
# November 24, 2012
# 
# Written by Justin Esarey, Rice University
#
#

cstarme<-function(sims, r){

	# translates regression results into a "loss point" using a simple kinked loss function
	# give it regression results, it returns a corresponding set of "loss points" t
	# if the point at which you begin making losses is less than t, accept evidence, else reject

	signind<-sign(mean(sims))
	if(signind==-1){sims<-(-1)*sims}
	utility<-function(x,t,a,sims){

		out<-c()
		for(i in 1:length(x)){
			k<-as.numeric((x[i]-t)>0)*2-1
			out[i]<-a^(-k)*density(sims, from=x[i],to=x[i], n=1)$y*(x[i]-t);
		}
		return(out)
	}

	maximand<-function(t,a,sims){as.numeric(integrate(utility, min(sims)-sd(sims), max(sims)+sd(sims), t=t, a=a, sims=sims, stop.on.error=FALSE)[1])}
	
	ans<-uniroot(maximand, interval=c(min(sims), max(sims)), a=r, sims=sims)$root
	if(signind==-1){ans<-(-1)*ans}

	if(signind!=sign(ans)){ans<-0}

	return(ans)
	
}


