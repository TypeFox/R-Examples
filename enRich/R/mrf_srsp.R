mrf_srsp<-function(data1, data2, nsr1, nsr2, nsp1, nsp2, method=method, Niterations=10000, Nburnin=5000, Poisprior1=NULL, NBprior1=NULL, PoisNBprior1=NULL, Poisprior2=NULL, NBprior2=NULL, PoisNBprior2=NULL, var.NB1=NULL, var.NB2=NULL, var.q=NULL)
{	
	data1=as.matrix(data1)
	data2=as.matrix(data2)
	Na=Niterations
	Nb=Nburnin
	jumpN=floor(Na/10)
	Nsample=ceiling((Na-Nb)/10)
	Nr1=nsr1+nsp1
	N1=Nsample*Nr1
	Nr2=nsr2+nsp2
	N2=Nsample*Nr2
	Nsr=c(nsr1, nsr2)
	Nsp=c(nsp1, nsp2)
	newdata1=as.vector(data1)
	newdata2=as.vector(data2)
	size=nrow(data1)
	PP1=numeric(size*Nr1)
	es_pi1=numeric(N1)
	es_q10=numeric(N1)
	es_q11=numeric(N1)
	es_lambda10=numeric(N1)
	es_lambda11=numeric(N1)
	es_mu10=numeric(N1)
	es_mu11=numeric(N1)
	es_phi10=numeric(N1)
	es_phi11=numeric(N1)
	PP2=numeric(size*Nr2)
	es_pi2=numeric(N2)
	es_q20=numeric(N2)
	es_q21=numeric(N2)
	es_lambda20=numeric(N2)
	es_lambda21=numeric(N2)
	es_mu20=numeric(N2)
	es_mu21=numeric(N2)
	es_phi20=numeric(N2)
	es_phi21=numeric(N2)
	loglikeli=numeric(Nsample)
	acrate1=numeric(Nr1*4)
	acrate2=numeric(Nr2*4)
	indexnsr1=ifelse(nsr1>0, 1, 0)
	indexnsr2=ifelse(nsr2>0, 1, 0)
	acrate=numeric(indexnsr1+indexnsr2+nsp1+nsp2+1)
	qprior1=rep(c(2.0, 2.0), nsp1+indexnsr1)
	qprior2=rep(c(2.0, 2.0), nsp2+indexnsr2)
	piprior1=rep(c(2.0,2.0), Nr1)
	piprior2=rep(c(2.0,2.0), Nr2)
	if (is.null(var.q))
	{
		var.q=c(rep(0.1, indexnsr1+indexnsr2+nsp1+nsp2), 0.3)	
	}
	if (method=="PoisNB")
	{
		met=0
		tPoisprior1=NULL
		tNBprior1=NULL
		for (i in 1:Nr1)
		{
			tPoisprior1=c(tPoisprior1, 0, 0, PoisNBprior1[((i-1)*6+1):((i-1)*6+2)])
			tNBprior1=c(tNBprior1, PoisNBprior1[((i-1)*6+3):((i-1)*6+6)], 0,0,0,0)
		}
		Poisprior1=tPoisprior1
		NBprior1=tNBprior1
		if (is.null(Poisprior1)==1)
		{
			Poisprior1=rep(c(0, 0, 0.5, 1), Nr1)
		}	
		if (is.null(NBprior1)==1)
		{
			NBprior1=rep(c(5, 1, 1, 1, 0, 0, 0, 0), Nr1)
		}
		if (is.null(var.NB1)==1)
		{
			var.NB1=rep(c(0.1, 0.1, 0, 0), Nr1)
		}	
		if (length(var.NB1)==2*Nr1)
		{
			tvar.NB1=NULL
			for (i in 1:Nr1)
			{
				tvar.NB1=c(tvar.NB1, var.NB1[((i-1)*2+1):((i-1)*2+2)], 0, 0)
			}
			var.NB1=tvar.NB1 
		}
		tPoisprior2=NULL
		tNBprior2=NULL
		for (i in 1:Nr2)
		{
			tPoisprior2=c(tPoisprior2, 0, 0, PoisNBprior2[((i-1)*6+1):((i-1)*6+2)])
			tNBprior2=c(tNBprior2, PoisNBprior2[((i-1)*6+3):((i-1)*6+6)], 0,0,0,0)
		}
		Poisprior2=tPoisprior2
		NBprior2=tNBprior2
		if (is.null(Poisprior2)==1)
		{
			Poisprior2=rep(c(0, 0, 0.5, 1), Nr2)
		}	
		if (is.null(NBprior2)==1)
		{
			NBprior2=rep(c(5, 1, 1, 1, 0, 0, 0, 0), Nr2)
		}
		if (is.null(var.NB2)==1)
		{
			var.NB2=rep(c(0.1, 0.1, 0, 0), Nr2)
		}	
		if (length(var.NB2)==2*Nr2)
		{
			tvar.NB2=NULL
			for (i in 1:Nr2)
			{
				tvar.NB2=c(tvar.NB2, var.NB2[((i-1)*2+1):((i-1)*2+2)], 0, 0)
			}
			var.NB2=tvar.NB2 
		}
	}	
	if (method=="Poisson")
	{
		met=1
		if (is.null(Poisprior1)==1)
		{
			Poisprior1=rep(c(5, 1, 0.5, 1), Nr1)
		}
		NBprior1=rep(0, 8*Nr1)
		var.NB1=rep(0, 4*Nr1)
		if (is.null(Poisprior2)==1)
		{
			Poisprior2=rep(c(5, 1, 0.5, 1), Nr2)
		}
		NBprior2=rep(0, 8*Nr2)
		var.NB2=rep(0, 4*Nr2)
	}
	if (method=="NB")
	{
		met=2
		Poisprior1=rep(0, 4*Nr1)
		if (is.null(NBprior1)==1)
		{
			NBprior1=rep(c(5, 1, 1, 1, 0.5, 1, 1, 1), Nr1)
		}
		if (is.null(var.NB1)==1)
		{
			var.NB1=rep(c(0.1, 0.1, 0.1, 0.1), Nr1)
		}	
		Poisprior2=rep(0, 4*Nr2)
		if (is.null(NBprior2)==1)
		{
			NBprior2=rep(c(5, 1, 1, 1, 0.5, 1, 1, 1), Nr2)
		}
		if (is.null(var.NB2)==1)
		{
			var.NB2=rep(c(0.1, 0.1, 0.1, 0.1), Nr2)
		}		
	}
	temp=.C("MRFsrsp", as.integer(newdata1), as.integer(newdata2), as.integer(Nsr), as.integer(Nsp), as.integer(size), as.integer(met), as.double(qprior1), as.double(piprior1), as.double(Poisprior1), as.double(NBprior1), as.double(qprior2), as.double(piprior2), as.double(Poisprior2), as.double(NBprior2), as.integer(Na), as.integer(Nb), as.integer(jumpN), as.double(var.q), as.double(var.NB1), as.double(var.NB2), PP1=as.double(PP1), PP2=as.double(PP2), es_q11=as.double(es_q11), es_q10=as.double(es_q10), es_lambda11=as.double(es_lambda11),  es_mu11=as.double(es_mu11), es_phi11=as.double(es_phi11), es_pi1=as.double(es_pi1), es_lambda10=as.double(es_lambda10), es_mu10=as.double(es_mu10), es_phi10=as.double(es_phi10), es_q21=as.double(es_q21), es_q20=as.double(es_q20), es_lambda21=as.double(es_lambda21),  es_mu21=as.double(es_mu21), es_phi21=as.double(es_phi21), es_pi2=as.double(es_pi2), es_lambda20=as.double(es_lambda20), es_mu20=as.double(es_mu20), es_phi20=as.double(es_phi20), loglikeli=as.double(loglikeli), acrate=as.double(acrate), acrate1=as.double(acrate1), acrate2=as.double(acrate2))
	q11=t(matrix(temp$es_q11, Nr1, Nsample))
	q10=t(matrix(temp$es_q10, Nr1, Nsample))	
	pi1=t(matrix(temp$es_pi1, Nr1, Nsample))
	mu11=t(matrix(temp$es_mu11, Nr1, Nsample))
	phi11=t(matrix(temp$es_phi11, Nr1, Nsample))
	mu10=t(matrix(temp$es_mu10, Nr1, Nsample))
	phi10=t(matrix(temp$es_phi10, Nr1, Nsample))
	lambda11=t(matrix(temp$es_lambda11, Nr1, Nsample))
	lambda10=t(matrix(temp$es_lambda10, Nr1, Nsample))
	q21=t(matrix(temp$es_q21, Nr2, Nsample))
	q20=t(matrix(temp$es_q20, Nr2, Nsample))	
	pi2=t(matrix(temp$es_pi2, Nr2, Nsample))
	mu21=t(matrix(temp$es_mu21, Nr2, Nsample))
	phi21=t(matrix(temp$es_phi21, Nr2, Nsample))
	mu20=t(matrix(temp$es_mu20, Nr2, Nsample))
	phi20=t(matrix(temp$es_phi20, Nr2, Nsample))
	lambda21=t(matrix(temp$es_lambda21, Nr2, Nsample))
	lambda20=t(matrix(temp$es_lambda20, Nr2, Nsample))
	PP1=matrix(temp$PP1, size, Nr1)
	PP2=matrix(temp$PP2, size, Nr2)
	acrate=temp$acrate
	acrate1=temp$acrate1
	acrate2=temp$acrate2
	para1=list()
	para2=list()
	if (met==0)
	{	
		tacrate=NULL
		for (i in 1:Nr1)
		{
			para1[[i]]=cbind(q1=q11[,i], q0=q10[,i], mu1=mu11[,i], phi1=phi11[,i],pi=pi1[,i], lambda0=lambda10[,i])
			tacrate=c(tacrate, temp$acrate1[c(i:(i+1))])
		}
		acrate1=tacrate
		tacrate=NULL
		for (i in 1:Nr2)
		{
			para2[[i]]=cbind(q1=q21[,i], q0=q20[,i], mu1=mu21[,i], phi1=phi21[,i],pi=pi2[,i], lambda0=lambda20[,i])
			tacrate=c(tacrate, temp$acrate2[c(i:(i+1))])
		}
		acrate2=tacrate
	}
	if (met==1)
	{
		for (i in 1:Nr1)
		{
			para1[[i]]=cbind(q1=q11[,i], q0=q10[,i], lambda1=lambda11[,i],pi=pi1[,i],  lambda0=lambda10[,i])
		}
		acrate1=NULL
		for (i in 1:Nr2)
		{
			para2[[i]]=cbind(q1=q21[,i], q0=q20[,i], lambda1=lambda21[,i],pi=pi2[,i],  lambda0=lambda20[,i])
		}
		acrate2=NULL
	}
	if (met==2)
	{
		for (i in 1:Nr1)
		{
			para1[[i]]=cbind(q1=q11[,i], q0=q10[,i],  mu1=mu11[,i], phi1=phi11[,i], pi=pi1[,i],mu0=mu10[,i], phi0=phi10[,i])
		}
		acrate1=temp$acrate1	
		for (i in 1:Nr2)
		{
			para2[[i]]=cbind(q1=q21[,i], q0=q20[,i],  mu1=mu21[,i], phi1=phi21[,i], pi=pi2[,i],mu0=mu20[,i], phi0=phi20[,i])
		}
		acrate2=temp$acrate2	
	}
	acrate.NB=c(acrate1, acrate2)
	results=list(para.con1=para1, PP.con1=PP1, para.con2=para2, PP.con2=PP2, method=method, acrate.NB=acrate.NB, acrate.q=acrate)
}

