Segmentor<- function(data=numeric(), model=1, Kmax = 15, phi = numeric(), m = numeric(), keep=FALSE, compress = TRUE) UseMethod("Segmentor")

Segmentor.default <-function(data=numeric(), model=1, Kmax = 15, phi = numeric(), m = numeric(), keep=FALSE, compress = TRUE)
{
  if ((model!=1)&(model!=2)&(model!=3)&(model!=4))
    stop("Choose model=1 (Poisson), 2 (normal), 3 (Negative Binomial) or 4 (Normal-Variance)")
  if (length(data)==0)
    stop("Give me a vector of data to segment")

	if (compress)
	{
		data2<-rle(data)
		dat<-data2$values
		datasize<-data2$lengths
		n2 = length(dat)
	} else
	{
		dat<-data
		n2 = length(dat)
		datasize<-rep(1,n2)	
	}
	n<-length(data)
  breaks=matrix(0,nrow=Kmax,ncol=Kmax)
  breaks = as.vector(breaks)
  parameters=matrix(0,nrow=Kmax,ncol=Kmax)
  parameters=as.vector(parameters)
  likelihood=rep(0, Kmax)
  likelihood=as.vector(likelihood)
  compression<-n/n2
  if ((model==4) & (length(m)==0))
	m = mean(data)
  if ((model==3) & (length(phi)==0))
  {
    h<-15
		Xcum = cumsum(data)
		X2cum = cumsum(data^2)
		M = (Xcum[h:n] - c(0, Xcum[1:(n-h)])) / h
		S2 = (X2cum[h:n] - c(0, X2cum[1:(n-h)])) / (h-1) - h/(h-1)*M^2
		K = M^2 / (S2-M)
		phi = median(K[!is.na(K)])
		while ((phi<0)&(h<(n/2)))
		{
			h<-2*h
			M = (Xcum[h:n] - c(0, Xcum[1:(n-h)])) / h
			S2 = (X2cum[h:n] - c(0, X2cum[1:(n-h)])) / (h-1) - h/(h-1)*M^2
			K = M^2 / (S2-M)
			phi = median(K[!is.na(K)])   	
		}
  }
  if(!keep)
  {
		if (model==1)
			Rep<-.C("SegmentPoisson", Size = as.integer(n2),KMax = as.integer(Kmax), Data = as.integer(dat), DataComp = as.integer(datasize), Breakpoints = as.integer(breaks), Parameters = as.double(parameters), Likelihood = as.double(likelihood), PACKAGE="Segmentor3IsBack") else if (model==3)
			Rep<-.C("SegmentBinNeg", Size = as.integer(n2),KMax = as.integer(Kmax), theta = as.double(phi), Data = as.integer(dat), DataComp = as.integer(datasize), Breakpoints = as.integer(breaks), Parameters = as.double(parameters), Likelihood = as.double(likelihood),  PACKAGE="Segmentor3IsBack") else if (model==2)
			Rep<-.C("SegmentNormal", Size = as.integer(n2),KMax = as.integer(Kmax), Data = as.double(dat), DataComp = as.integer(datasize), Breakpoints = as.integer(breaks), Parameters = as.double(parameters), Likelihood = as.double(likelihood),  PACKAGE="Segmentor3IsBack") else if (model==4)
			Rep<-.C("SegmentVariance", Size = as.integer(n2),KMax = as.integer(Kmax), mu = as.double(m), Data = as.double(dat), DataComp = as.integer(datasize), Breakpoints = as.integer(breaks), Parameters = as.double(parameters), Likelihood = as.double(likelihood),  PACKAGE="Segmentor3IsBack")
	} else
	{
		  cost=matrix(0,nrow=Kmax,ncol=n2)
  		cost=as.vector(cost)
  		pos=matrix(0,nrow=Kmax,ncol=n2)
  		pos=as.vector(pos)
		if (model==1)
			Rep<-.C("SegmentPoissonKeep", Size = as.integer(n2),KMax = as.integer(Kmax), Data = as.integer(dat), DataComp = as.integer(datasize), Breakpoints = as.integer(breaks), Parameters = as.double(parameters), Likelihood = as.double(likelihood), Cost = as.double(cost), Pos = as.integer(pos), PACKAGE="Segmentor3IsBack") else if (model==3)
			Rep<-.C("SegmentBinNegKeep", Size = as.integer(n2),KMax = as.integer(Kmax), theta = as.double(phi), Data = as.integer(dat), DataComp = as.integer(datasize), Breakpoints = as.integer(breaks), Parameters = as.double(parameters), Likelihood = as.double(likelihood), Cost = as.double(cost), Pos = as.integer(pos), PACKAGE="Segmentor3IsBack") else if (model==2)
			Rep<-.C("SegmentNormalKeep", Size = as.integer(n2),KMax = as.integer(Kmax), Data = as.double(dat), DataComp = as.integer(datasize), Breakpoints = as.integer(breaks), Parameters = as.double(parameters), Likelihood = as.double(likelihood), Cost = as.double(cost), Pos = as.integer(pos), PACKAGE="Segmentor3IsBack") else if (model==4)
			Rep<-.C("SegmentVarianceKeep", Size = as.integer(n2),KMax = as.integer(Kmax), mu = as.double(m), Data = as.double(dat), DataComp = as.integer(datasize), Breakpoints = as.integer(breaks), Parameters = as.double(parameters), Likelihood = as.double(likelihood), Cost = as.double(cost), Pos = as.integer(pos), PACKAGE="Segmentor3IsBack")	
			cost = matrix(Rep$Cost,ncol=Kmax)
 		  cost = t(cost)
    	pos = matrix(Rep$Pos,ncol=Kmax)
  		pos = t(pos)
	}
  breaks=matrix(Rep$Breakpoints,ncol=Kmax)
  breaks = t(breaks)
  for (k in 1:length(breaks[,1]))
  {
  	breaks[k,1:k]<-cumsum(datasize)[breaks[k,1:k]]
  }
  parameters = matrix(Rep$Parameters,ncol=Kmax)
  parameters = t(parameters)
  rownames(breaks)<-c("1 segment", paste(2:Kmax, "segments"))
  rownames(parameters)<-c("1 segment", paste(2:Kmax, "segments"))
  colnames(breaks)<-c(paste(1:Kmax, "th break",sep=""))
  colnames(parameters)<-c(paste(1:Kmax, "th parameter",sep=""))
  likelihood=matrix(Rep$Likelihood,ncol=1)
  rownames(likelihood)<-c("1 segment", paste(2:Kmax, "segments"))
  if (!keep)
  {
		if (model==1) 
		{
		    model.dist="Poisson"
		    likelihood=likelihood+sum(lgamma(data+1))
		    Segmentor.res=new("Segmentor", data=data, model=model.dist, breaks=breaks, parameters=parameters, likelihood=likelihood, Kmax=Kmax, compression=compression)
		}
		if (model==2) 
		{
			data1<-data[-(1:3)]
			data2<-data[-c(1,2,n)]
			data3<-data[-c(1,n-1,n)]
			data4<-data[-c(n-2,n-1,n)]
			d<-c(0.1942, 0.2809, 0.3832, -0.8582)
			v2<-d[1]*data1+d[2]*data2+d[3]*data3+d[4]*data4
			v2<-v2*v2
			var<-sum(v2)/(n-3)
		  likelihood = likelihood /(2*var) + n/2*log(2*pi*var)
		  model.dist="Normal"
		  Segmentor.res=new("Segmentor", data=data, model=model.dist, breaks=breaks, parameters=parameters, likelihood=likelihood, Kmax=Kmax, compression=compression)
		}
		if (model==3) 
		{
		    model.dist="Negative binomial"
		    Segmentor.res=new("Segmentor",data=data, model=model.dist, breaks=breaks, overdispersion=phi, parameters=parameters, likelihood=likelihood, Kmax=Kmax, compression=compression)
		}
		if (model==4) 
		{
		    likelihood = likelihood + n/2*log(2*pi)
		    model.dist="Variance Segmentation"
		    Segmentor.res=new("Segmentor",data=data, model=model.dist, breaks=breaks, mean=m, parameters=parameters, likelihood=likelihood, Kmax=Kmax, compression=compression)
		}
	} else
	{
		if (model==1) 
		{
		    model.dist="Poisson"
		    likelihood=likelihood+sum(lgamma(data+1))
		    Segmentor.res=new("Segmentor", data=data, model=model.dist, breaks=breaks, parameters=parameters, likelihood=likelihood, Kmax=Kmax, Cost=cost, Pos=pos, compression=compression)
		}
		if (model==2) 
		{
			data1<-data[-(1:3)]
			data2<-data[-c(1,2,n)]
			data3<-data[-c(1,n-1,n)]
			data4<-data[-c(n-2,n-1,n)]
			d<-c(0.1942, 0.2809, 0.3832, -0.8582)
			v2<-d[1]*data1+d[2]*data2+d[3]*data3+d[4]*data4
			v2<-v2*v2
			var<-sum(v2)/(n-3)
		  likelihood = likelihood /(2*var) + n/2*log(2*pi*var)
		  model.dist="Normal"
		  Segmentor.res=new("Segmentor",data=data, model=model.dist, breaks=breaks, parameters=parameters,likelihood=likelihood,Kmax=Kmax, Cost=cost, Pos=pos, compression=compression)
		}
		if (model==3) 
		{
		    model.dist="Negative binomial"
		    Segmentor.res=new("Segmentor",data=data, model=model.dist, breaks=breaks, overdispersion=phi, parameters=parameters, likelihood=likelihood, Kmax=Kmax, Cost=cost, Pos=pos, compression=compression)
		}
		if (model==4) 
		{
		    likelihood = likelihood + n/2*log(2*pi)
		    model.dist="Variance Segmentation"
		    Segmentor.res=new("Segmentor",data=data, model=model.dist, breaks=breaks, mean=m, parameters=parameters, likelihood=likelihood, Kmax=Kmax, Cost=cost, Pos=pos, compression=compression)	
	
	}
	}
  Segmentor.res

}

############################################################################################

BestSegmentation <- function(x,K,t=numeric(),compress=TRUE)
{
	if (class(x)!="Segmentor")
		stop("x must be an object of class Segmentor returned by the Segmentor function")
	if (dim(getCost(x))[1]==0)
		stop('to use this function the data must have been processed with option keep=TRUE')
	if (K==1)
		stop("There exist only one segmentation with 0 break")
	n<-length(getData(x))
	data<-getData(x)
	s<-getModel(x)
	if (s=="Poisson") mod=1 else if (s=="Normal") mod=2 else if (s=="Negative binomial") mod=3 else if (s=="Variance Segmentation") mod=4
	if (K<=getKmax(x))
	{
		c<-getCompression(x)
		data2<-rle(data)
		dat<-data2$values
		datasize<-data2$lengths
		n2 = length(dat)
		if (compress & (c==1) &(n2<n))
		{		
			warning('Warning: Segmentor was applied with compression = FALSE, will apply new Segmentor with compression = TRUE')
			resForward <- Segmentor(data, model=mod, Kmax=K, phi=getOverdispersion(x), m=getMean(x), keep=TRUE, compress = compress)		
		} else if ((!compress) & (c>1))
		{
			warning('Warning: Segmentor was applied with compression = TRUE, will apply new Segmentor with compression = FALSE')
			resForward <- Segmentor(data, model=mod, Kmax=K, phi=getOverdispersion(x), m=getMean(x), keep=TRUE, compress = compress)	
			n2<-n
		} else
		{
			resForward <- x
			n2<-n/c
		}
	} else
	{
		warning('Warning: K is greater than Kmax, will apply new Segmentor with Kmax=K')
		resForward <- Segmentor(data, model=mod, Kmax=K, phi=getOverdispersion(x), m=getMean(x), keep=TRUE, compress = compress)
	}
	resBack <- Segmentor(rev(data), model=mod, Kmax=K, phi=getOverdispersion(x), m=getMean(x), keep=TRUE, compress = compress)
	bestCost <- t(matrix( getCost(resForward)[1:(K-1), 1:(n2-1)] + getCost(resBack)[(K-1):1, (n2-1):1], nrow=K-1))
	MF<-getPos(resForward)
	if (length(t)!=0)
	{
		MB<-getPos(resBack)
		k1<-which.min(bestCost[t,])
		k2<-K-k1; t2<-n2-t;
		TheBreakpoints<-NULL
		if (k1>1)
		{
			Prec = MF[k1,t];
			TheBreakpoints<-c(TheBreakpoints,(Prec+1));
			if (k1>2)
				for (i in (k1-1):2)
				{
					TheBreakpoints<-c(TheBreakpoints,MF[i,Prec]+1);
					Prec = MF[i,Prec];
				}
		}
		if (k2>1)
		{
			Prec = MB[k2,t2];
			TheBreakpoints<-c(TheBreakpoints,(n2-Prec-1));
			if (k2>2)
				for (i in (k2-1):2)
				{
					Prec = MB[i,Prec];
					TheBreakpoints<-c(TheBreakpoints,n2-Prec-1);	
				}
		}
		TheBreakpoints<-sort(c(TheBreakpoints,1,t,n2))
		bestCost.res<-list(bestCost=bestCost,bestSeg=TheBreakpoints)
	} else bestCost.res<-list(bestCost=bestCost)
	bestCost.res
}

SelectModel <-function(x,penalty="oracle",seuil=n/log(n),keep=FALSE,greatjump=FALSE)
{
	if ((penalty!='BIC') & (penalty!='mBIC') & (penalty!='AIC') & (penalty!='oracle'))
		stop("penalty must be BIC, mBIC, AIC or oracle")
	if (class(x)!="Segmentor")
		stop("x must be an object of class Segmentor returned by the Segmentor function")
	n<-getBreaks(x)[1,1]
	Kmax<-getKmax(x)
	sizenr<-function(k) {	vec<-sum(log(diff(c(1,getBreaks(x)[k,1:k]))))}
	saut<-function(Lv, pen,Kseq,seuil=sqrt(n)/log(n),biggest=TRUE)
	{
		J=-Lv;Kmax=length(J); k=1;kv=c();dv=c();pv=c();dmax=1
		while (k<Kmax) {
				pk=(J[(k+1):Kmax]-J[k])/(pen[k]-pen[(k+1):Kmax])
				pm=max(pk); dm=which.max(pk); dv=c(dv,dm); kv=c(kv,k); pv=c(pv,pm)
				if (dm>dmax){  
				  dmax=dm; kmax=k; pmax=pm  
				  }
				k=k+dm
		 } 
		if (biggest)
		{
			pv=c(pv,0); kv=c(kv,Kmax); dv=diff(kv); dmax=max(dv); rt=max(dv); rt=which(dv==rt)
			pmax=pv[rt[length(rt)]]
			alpha=2*pmax
			km=kv[alpha>=pv]; Kh =Kseq[km[1]] 
			return(c(Kh,alpha))
		} else
		{
			paux<-pv[which(kv<=seuil)]
			alpha<-2*min(paux)
			km=kv[alpha>=pv];	Kh =Kseq[km[1]] 
			return(c(Kh,alpha))	
		}
	}
	if (getModel(x)=="Poisson")
	{
		if(penalty=='mBIC')
			K<-which.min(crit<-getLikelihood(x)+0.5*sapply(1:Kmax,sizenr)+(1:Kmax-0.5)*log(n))
		if (penalty=='BIC')
			K<-which.min(crit<-getLikelihood(x)+1:K*log(n))
		if (penalty=='AIC')
			K<-which.min(crit<-getLikelihood(x)+1:K*2)
		if (penalty=='oracle')
		{
			Kseq=1:Kmax
			pen=Kseq*(1+4*sqrt(1.1+log(n/Kseq)))*(1+4*sqrt(1.1+log(n/Kseq)))
			if(greatjump)
			{
				K=saut(-getLikelihood(x)[Kseq],pen,Kseq)
				crit<-getLikelihood(x)[Kseq]+K[2]*pen
				K<-K[1]
			} else
			{
				K=saut(-getLikelihood(x)[Kseq],pen,Kseq,seuil,biggest=FALSE)
				crit<-getLikelihood(x)[Kseq]+K[2]*pen
				K<-K[1]
			}
		}	
	} else if (getModel(x)=="Negative binomial")
	{
		if(penalty=='mBIC')
			stop("no mBIC for Negative Binomial model")
		if (penalty=='BIC')
			K<-which.min(getLikelihood(x)+((1:Kmax)+1)*log(n))
		if (penalty=='AIC')
			K<-which.min(getLikelihood(x)+((1:Kmax)+1)*2)	
		if (penalty=='oracle')
		{
			Kseq=1:Kmax
			pen=Kseq*(1+4*sqrt(1.1+log(n/Kseq)))*(1+4*sqrt(1.1+log(n/Kseq)))
			if(greatjump)
			{
				K=saut(-getLikelihood(x)[Kseq],pen,Kseq)
				crit<-getLikelihood(x)[Kseq]+K[2]*pen
				K<-K[1]
			} else
			{
				K=saut(-getLikelihood(x)[Kseq],pen,Kseq,seuil,biggest=FALSE)
				crit<-getLikelihood(x)[Kseq]+K[2]*pen
				K<-K[1]
			}
		}	
	} else if (getModel(x)=='Normal')
	{
		if(penalty=='mBIC')
			K<-which.min(getLikelihood(x)+0.5*sapply(1:Kmax,sizenr)+(1:Kmax-0.5)*log(n))
		if (penalty=='BIC')
			K<-which.min(getLikelihood(x)+((1:Kmax)+1)*log(n))
		if (penalty=='AIC')
			K<-which.min(getLikelihood(x)+((1:Kmax)+1)*2)		
		if (penalty=='oracle')
		{
			Kseq=1:Kmax
			pen=Kseq*(2*log(n/Kseq)+5)
			if(greatjump)
			{
				K=saut(-getLikelihood(x)[Kseq],pen,Kseq)
				crit<-getLikelihood(x)[Kseq]+K[2]*pen
				K<-K[1]
			} else
			{
				K=saut(-getLikelihood(x)[Kseq],pen,Kseq,seuil,biggest=FALSE)
				crit<-getLikelihood(x)[Kseq]+K[2]*pen
				K<-K[1]
			}
		}	
	} else
	{
		if(penalty=='mBIC')
			stop("no mBIC for Variance model")
		if (penalty=='BIC')
			K<-which.min(getLikelihood(x)+((1:K)+1)*log(n))
		if (penalty=='AIC')
			K<-which.min(getLikelihood(x)+((1:K)+1)*2)		
		if(penalty=='oracle')
			stop("no oracle penalty for Variance model")
	}
	if (keep)
		res<-list(K=K,criterion=crit)
	else res<-K
	return(res)
}



print.Segmentor <-function(x,...)
{
  cat("\n Model used for the segmentation: \n")
  print(getModel(x))
  
  cat("\n Table of optimal breakpoints: \n")
  print(getBreaks(x))

  cat("\n Table of negative log-likelihood for each optimal segmentation: \n")
  print(getLikelihood(x))

  if (is.element("overdispersion",names(x)))
  {
    cat("\n Value of the overdispersion used for the segmentation: \n")
    print(getOverdispersion(x))
  }  
  if (is.element("mean",names(x)))
  {
    cat("\n Value of the mean used for the segmentation: \n")
    print(getMean(x))
  }  
  if (is.element("parameters",names(x)))
  {
    cat("\n Table of parameters: ")
    if (getModel(x)=="Poisson")
      cat("(mean of the signal in each segment) \n")
    if (getModel(x)=="Normal")
      cat("(mean of the signal in each segment) \n")
    if (getModel(x)=="Negative binomial")
      cat("(success-probability of the signal in each segment) \n")
    if (getModel(x)=="Variance Segmentation")
      cat("(Variance of the signal in each segment) \n")
    print(getParameters(x))
  }
}

