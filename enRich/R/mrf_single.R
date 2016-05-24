mrf_single<-function(data, met, Niterations=10000, Nburnin=5000, Poisprior=c(5, 1, 0.5, 1), NBprior=c(5, 1, 1, 1, 0.5, 1, 1, 1), PoisNBprior=c(5,1,1,1, 0.5,1), var.NB=c(0.1, 0.1, 0.1, 0.1))
{	
## INPUT
## data: the counts of a single ChIP experiment, which is a vector of size n. where n is the number of regions and p is the number of experiments. 
## met: The code denote the method, be 0 for "PoisNB", 1 for "Poisson" and 2 for "NB" and it refers to the densities of the mixture distribution.
## Niterations: the number of MCMC iteration steps. 
## Nburnin: the number of burn-in steps.
## Poisprior: the gamma priors for mean parameter lambda in Poisson-Poisson mixture, the first two are priors for signal and the second two are priors for background. 
##                 Default values are (5,1, 0.5, 1). 
## NBprior: the gamma priors for mean mu and overdispersion parameters phi in NB-NB mixture, the first two are priors for mu_S for signal, the third and fourth are priors for phi_S 
##                    the fifth and sixth are priors for mu_B of background and the seventh and eighth are priors for phi_B. Default values are (5, 1, 1, 1, 0.5, 1, 1, 1).
## PoisNBprior: the gamma priors for lambda_B and mu_S, phi_S in Poisson-NB mixture, the first two are priors for mu_S, the third and the fourth are priors for phi_S, 
##                      the fifth and the sixth are priors for lambda_B. Default values are (5, 1,1,1, 0.5, 1). 
## var.NB: the variances used in Metropolis-Hastling algorithm for estimates of (mu_S, phi_S, mu_B, phi_B) for NB mixture or for estimates of (mu_S, phi_S) for poisNB mixture. 
##              Default values are (0.1, 0.1, 0.1, 0.1) or (0.1, 0.1) for NB and poisNB respectively. 

	data1=data
	qprior=c(2.0, 2.0, 2.0, 2.0)
	piprior=c(2.0,2.0)
	Na=Niterations
	Nb=Nburnin
	jumpN=floor(Na/10)
	N=ceiling((Na-Nb)/10)
	size=length(data1)
	PP=numeric(size)
	es_pi=numeric(N)
	es_q0=numeric(N)
	es_q1=numeric(N)
	es_lambda0=numeric(N)
	es_lambda1=numeric(N)
	es_mu0=numeric(N)
	es_mu1=numeric(N)
	es_phi0=numeric(N)
	es_phi1=numeric(N)
	loglikeli=numeric(N)
	acrate=numeric(4)	
	temp=.C("MRF", as.integer(data1), as.integer(size), as.integer(met), as.double(qprior), as.double(piprior), as.double(Poisprior), as.double(NBprior), as.integer(Na), as.integer(Nb), as.integer(jumpN), as.double(var.NB), PP=as.double(PP), es_pi=as.double(es_pi), es_q1=as.double(es_q1), es_q0=as.double(es_q0), es_lambda1=as.double(es_lambda1),  es_mu1=as.double(es_mu1), es_phi1=as.double(es_phi1), es_lambda0=as.double(es_lambda0), es_mu0=as.double(es_mu0), es_phi0=as.double(es_phi0), loglikeli=as.double(loglikeli), acrate=as.double(acrate), PACKAGE = "enRich")
	if (met==0)
	{
		para.sample=cbind( q1=temp$es_q1, q0=temp$es_q0, mu_S=temp$es_mu1, phi_S=temp$es_phi1, pi=temp$es_pi, lambda_B=temp$es_lambda0)
		acrate=temp$acrate[c(1:2)]
	}
	if (met==1)
	{
		para.sample=cbind(q1=temp$es_q1, q0=temp$es_q0, lambda_S=temp$es_lambda1, pi=temp$es_pi, lambda_B=temp$es_lambda0)
		acrate=NULL
	}
	if (met==2)
	{
		para.sample=cbind(q1=temp$es_q1, q0=temp$es_q0, mu_S=temp$es_mu1, phi_S=temp$es_phi1, pi=temp$es_pi,mu_B=temp$es_mu0, phi_B=temp$es_phi0)
		acrate=temp$acrate
	}
	para=round(apply(para.sample, 2, mean), 4)
	para=t(as.matrix(para))
	results=list(data=data, parameters=para, parameters.sample=para.sample, PP=temp$PP, acrate.NB=acrate)
	return(results)
}
