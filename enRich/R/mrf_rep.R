mrf_rep<-function(data1, method=method, Niterations=10000, Nburnin=5000, Poisprior=NULL, NBprior=NULL, PoisNBprior=NULL, var.NB=NULL)
{	
	data1=as.matrix(data1)
	Nrep=ncol(data1)#we expect have p experiment, so data1 is in form of matrix with n row and p col.
	Na=Niterations
	Nb=Nburnin
	jumpN=floor(Na/10)
	Nsample=ceiling((Na-Nb)/10)
	N=Nsample*Nrep
	newdata1=as.vector(data1)
	size=nrow(data1)
	PP=numeric(size)
	es_pi=numeric(N)
	es_q0=numeric(Nsample)
	es_q1=numeric(Nsample)
	es_lambda0=numeric(N)
	es_lambda1=numeric(N)
	es_mu0=numeric(N)
	es_mu1=numeric(N)
	es_phi0=numeric(N)
	es_phi1=numeric(N)
	loglikeli=numeric(Nsample)
	acrate=numeric(4*Nrep)
	qprior=c(2.0, 2.0, 2.0, 2.0)
	piprior=rep(c(2.0,2.0), Nrep)
	if (method=="PoisNB")
	{
		met=0
		Poisprior1=NULL
		NBprior1=NULL
		for (i in 1:Nrep)
		{
			Poisprior1=c(Poisprior1, 0, 0, PoisNBprior[((i-1)*6+1):((i-1)*6+2)])
			NBprior1=c(NBprior1, PoisNBprior[((i-1)*6+3):((i-1)*6+6)], 0,0,0,0)
		}
		Poisprior=Poisprior1
		NBprior=NBprior1
		if (is.null(Poisprior))
		{
			Poisprior=rep(c(0, 0, 0.5, 1), Nrep)
		}	
		if (is.null(NBprior))
		{
			NBprior=rep(c(5, 1, 1, 1, 0, 0, 0, 0), Nrep)
		}
		if (is.null(var.NB))
		{
			var.NB=rep(c(0.1, 0.1, 0, 0), Nrep)
		}	
		if (length(var.NB)==2*Nrep)
		{
			tvar.NB=NULL
			for (i in 1:Nrep)
			{
				tvar.NB=c(tvar.NB, var.NB[((i-1)*2+1):((i-1)*2+2)], 0, 0)
			}
			var.NB=tvar.NB 
		}
	}	
	if (method=="Poisson")
	{
		met=1
		if (is.null(Poisprior))
		{
			Poisprior=rep(c(5, 1, 0.5, 1), Nrep)
		}
		NBprior=rep(0, 8*Nrep)
		var.NB=rep(0, 4*Nrep)
	}
	if (method=="NB")
	{
		met=2
		Poisprior=rep(0, 4*Nrep)
		if (is.null(NBprior))
		{
			NBprior=rep(c(5, 1, 1, 1, 0.5, 1, 1, 1), Nrep)
		}
		if (is.null(var.NB))
		{
			var.NB=rep(c(0.1, 0.1, 0.1, 0.1), Nrep)
		}		
	}
	temp=.C("MRFrep", as.integer(newdata1), as.integer(Nrep), as.integer(size), as.integer(met), as.double(qprior), as.double(piprior), as.double(Poisprior), as.double(NBprior), as.integer(Na), as.integer(Nb), as.integer(jumpN), as.double(var.NB), PP=as.double(PP), es_pi=as.double(es_pi), es_q1=as.double(es_q1), es_q0=as.double(es_q0), es_lambda1=as.double(es_lambda1),  es_mu1=as.double(es_mu1), es_phi1=as.double(es_phi1), es_lambda0=as.double(es_lambda0), es_mu0=as.double(es_mu0), es_phi0=as.double(es_phi0), loglikeli=as.double(loglikeli), acrate=as.double(acrate), PACKAGE = "enRich")
	pi=t(matrix(temp$es_pi, Nrep, Nsample))
	mu1=t(matrix(temp$es_mu1, Nrep, Nsample))
	phi1=t(matrix(temp$es_phi1, Nrep, Nsample))
	mu0=t(matrix(temp$es_mu0, Nrep, Nsample))
	phi0=t(matrix(temp$es_phi0, Nrep, Nsample))
	lambda1=t(matrix(temp$es_lambda1, Nrep, Nsample))
	lambda0=t(matrix(temp$es_lambda0, Nrep, Nsample))
	q1=temp$es_q1
	q0=temp$es_q0
	para=list()
	if (met==0)
	{	
		for (i in 1:Nrep)
		{
			para[[i]]=cbind(q1=q1, q0=q0, mu1=mu1[,i], phi1=phi1[,i],pi=pi[,i], lambda0=lambda0[,i])
			acrate=temp$acrate[c(1:2)]
		}
	}
	if (met==1)
	{
		for (i in 1:Nrep)
		{
			para[[i]]=cbind(q1=q1, q0=q0, lambda1=lambda1[,i],pi=pi[,i],  lambda0=lambda0[,i])
			acrate=NULL
		}
	}
	if (met==2)
	{
		for (i in 1:Nrep)
		{
			para[[i]]=cbind(q1=q1, q0=q0, mu1=mu1[,i], phi1=phi1[,i], pi=pi[,i],mu0=mu0[,i], phi0=phi0[,i])
			acrate=temp$acrate
		}
	}
	results=list(parameters=para, PP=temp$PP, method=method, acrate.NB=acrate)
}

