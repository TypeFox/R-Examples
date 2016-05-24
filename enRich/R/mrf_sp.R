mrf_sp<-function(data1, method=method, Niterations=10000, Nburnin=5000, Poisprior=NULL, NBprior=NULL, PoisNBprior=NULL, var.NB=NULL, var.q=NULL)
{	
	data1=as.matrix(data1)
	Nsp=ncol(data1)#we expect have p experiment, so data1 is in form of matrix with n row and p col.
	Na=Niterations
	Nb=Nburnin
	jumpN=floor(Na/10)
	Nsample=ceiling((Na-Nb)/10)
	N=Nsample*Nsp
	newdata1=as.vector(data1)
	size=nrow(data1)
	PP=numeric(size*Nsp)
	es_pi=numeric(N)
	es_q0=numeric(N)
	es_q1=numeric(N)
	es_lambda0=numeric(N)
	es_lambda1=numeric(N)
	es_mu0=numeric(N)
	es_mu1=numeric(N)
	es_phi0=numeric(N)
	es_phi1=numeric(N)
	loglikeli=numeric(Nsample)
	acrate=numeric(4*Nsp)
	acrate1=numeric(Nsp+1)
	qprior=rep(c(2.0, 2.0), Nsp)
	piprior=rep(c(2.0,2.0), Nsp)
	if (is.null(var.q))
	{
		var.q=c(rep(0.1, Nsp), 0.3)
	}
	if (method=="PoisNB")
	{
		met=0
		Poisprior1=NULL
		NBprior1=NULL
		for (i in 1:Nsp)
		{
			Poisprior1=c(Poisprior1, 0, 0, PoisNBprior[((i-1)*6+1):((i-1)*6+2)])
			NBprior1=c(NBprior1, PoisNBprior[((i-1)*6+3):((i-1)*6+6)], 0,0,0,0)
		}
		Poisprior=Poisprior1
		NBprior=NBprior1
		if (is.null(Poisprior)==1)
		{
			Poisprior=rep(c(0, 0, 0.5, 1), Nsp)
		}	
		if (is.null(NBprior)==1)
		{
			NBprior=rep(c(5, 1, 1, 1, 0, 0, 0, 0), Nsp)
		}
		if (is.null(var.NB)==1)
		{
			var.NB=rep(c(0.1, 0.1, 0, 0), Nsp)
		}	
		if (length(var.NB)==2*Nsp)
		{
			tvar.NB=NULL
			for (i in 1:Nsp)
			{
				tvar.NB=c(tvar.NB, var.NB[((i-1)*2+1):((i-1)*2+2)], 0, 0)
			}
			var.NB=tvar.NB 
		}
	}	
	if (method=="Poisson")
	{
		met=1
		if (is.null(Poisprior)==1)
		{
			Poisprior=rep(c(5, 1, 0.5, 1), Nsp)
		}
		NBprior=rep(0, 8*Nsp)
		var.NB=rep(0, 4*Nsp)			
	}
	if (method=="NB")
	{
		met=2
		Poisprior=rep(0, 4*Nsp)
		if (is.null(NBprior)==1)
		{
			NBprior=rep(c(5, 1, 1, 1, 0.5, 1, 1, 1), Nsp)
		}
		if (is.null(var.NB)==1)
		{
			var.NB=rep(c(0.1, 0.1, 0.1, 0.1), Nsp)
		}		
	}
	temp=.C("MRFsp", as.integer(newdata1), as.integer(Nsp), as.integer(size), as.integer(met), as.double(qprior), as.double(piprior), as.double(Poisprior), as.double(NBprior), as.integer(Na), as.integer(Nb), as.integer(jumpN), as.double(var.NB), as.double(var.q), PP=as.double(PP), es_pi=as.double(es_pi), es_q1=as.double(es_q1), es_q0=as.double(es_q0), es_lambda1=as.double(es_lambda1),  es_mu1=as.double(es_mu1), es_phi1=as.double(es_phi1), es_lambda0=as.double(es_lambda0), es_mu0=as.double(es_mu0), es_phi0=as.double(es_phi0), loglikeli=as.double(loglikeli), acrate=as.double(acrate), acrate1=as.double(acrate1))
	pi=t(matrix(temp$es_pi, Nsp, Nsample))
	mu1=t(matrix(temp$es_mu1, Nsp, Nsample))
	phi1=t(matrix(temp$es_phi1, Nsp, Nsample))
	mu0=t(matrix(temp$es_mu0, Nsp, Nsample))
	phi0=t(matrix(temp$es_phi0, Nsp, Nsample))
	lambda1=t(matrix(temp$es_lambda1, Nsp, Nsample))
	lambda0=t(matrix(temp$es_lambda0, Nsp, Nsample))
	q1=t(matrix(temp$es_q1, Nsp, Nsample))
	q0=t(matrix(temp$es_q0, Nsp, Nsample))
	PP=matrix(temp$PP, size, Nsp)
	acrate1=temp$acrate1
	para=list()
	if (met==0)
	{	
		tacrate=NULL
		for (i in 1:Nsp)
		{
			para[[i]]=cbind(q1=q1[,i], q0=q0[,i], mu1=mu1[,i], phi1=phi1[,i],pi=pi[,i], lambda0=lambda0[,i])
			tacrate=c(tacrate, temp$acrate[c(i:(i+1))])
		}
		acrate=tacrate
	}
	if (met==1)
	{
		for (i in 1:Nsp)
		{
			para[[i]]=cbind(q1=q1[,i], q0=q0[,i], lambda1=lambda1[,i],pi=pi[,i],  lambda0=lambda0[,i])
		}
		acrate=NULL
	}
	if (met==2)
	{
		for (i in 1:Nsp)
		{
			para[[i]]=cbind(q1=q1[,i], q0=q0[,i],  mu1=mu1[,i], phi1=phi1[,i], pi=pi[,i],mu0=mu0[,i], phi0=phi0[,i])
		}
		acrate=temp$acrate
	}
	results=list(parameters=para, PP=PP, method=method, acrate.NB=acrate, acrate.q=acrate1)
}

