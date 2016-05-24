## -----------------------------------------------------------------------------
## Fonction MetropolisHastings
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

MetropolisHastings = function(q=Uniform,p=dmvnorm,x0,eval_x0=-1,chain_length,modified=FALSE,modif_parameter=1,limit_fun = function(x) {-1},burnin=30,thinning=4) {
#q is the proposed PDF and y=q(x) should generate a random y|x while q(x,y) should return P(y|x)
#p is the target PDF
#x0 is the seed
#chain_length = how many samples generated, at the end the chain will be chain_length+1 long


niter = (thinning+1)*chain_length + burnin
d = length(x0)
U = matrix(nrow=d,ncol=(niter+1))
U[,1] = x0
MC = list(points=U,eval=eval_x0,acceptation=NA,Ncall=NULL,samples=matrix(nrow=d),eval_samples=NA);
tau = 0
Ncall = 0;

if(modified==FALSE){

	for (i in 1:niter) {
		candidate = q(MC$points[,i])
		acceptance = min(modif_parameter*(p(candidate)*q(candidate,MC$points[,i]))/(p(MC$points[,i])*q(MC$points[,i],candidate)),1)
		
		rand = runif(1,min=0,max=1)
		test = (rand-acceptance)<0
		if (!test){
			candidate = MC$points[,i];
			eval = MC$eval[i];
			tau = tau+1
		}
		else {
			eval = limit_fun(candidate); Ncall = Ncall+1; MC$samples = cbind(MC$samples,candidate); MC$eval_samples = c(MC$eval_samples,eval);
			indicatrice = (1-sign(eval))/2 #1 if limit_fun(x) < 0, 0 otherwise
			if (indicatrice==0 | is.nan(indicatrice)) {
				candidate = MC$points[,i];
				eval = MC$eval[i]
				tau = tau+1
			}
		}
		MC$points[,i+1] = candidate;
		MC$eval[i+1] = eval
	}
}

#modified algorithm works coordinate by coordinate
if(modified==TRUE){
	
	for (i in 1:niter) {
		candidate = NA*c(1:d)
		for(j in 1:d){
			seed = MC$points[,i][j]
			candidate[j] = q(seed)
			acceptance = min(1,
				modif_parameter*(p(candidate[j])*q(candidate[j],seed))/(p(seed)*q(seed,candidate[j])))
			rand = runif(1,min=0,max=1)
			test = (rand-acceptance)<0
			if (!test) {candidate[j] = seed;}
		}
		eval = limit_fun(candidate); Ncall = Ncall+1; MC$samples = cbind(MC$samples,candidate); MC$eval_samples = c(MC$eval_samples,eval);
		indicatrice = (1-sign(eval))/2 #1 if limit_fun(x) < 0, 0 otherwise
		if (indicatrice==0 | is.nan(indicatrice)) {
			candidate = MC$points[,i];
			eval = MC$eval[i]
# 			tau = tau+1
		}
		MC$points[,i+1] = candidate;
		MC$eval[i+1] = eval
	}
}
tau = tau/niter
MC$acceptation = 1-tau
sel_samples = burnin+1+(thinning+1)*c(0:chain_length)
MC$points = MC$points[,sel_samples]
MC$eval = MC$eval[sel_samples]
MC$Ncall = Ncall

return(MC)

}