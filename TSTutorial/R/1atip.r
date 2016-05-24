########################## ATIPICS2.R ############################


outdetec<-function(object,dif=c(0,0),crit,LS=T){
	residuals<-object$residuals
	m<-length(residuals)

	piweight<--ARMAtoMA(ar=-object$model$theta, ma=-object$model$phi, lag.max=m+sum(dif))

	if (dif[1]!=0) for(i in 1:dif[1]) piweight<-c(piweight,0)-c(-1,piweight)
	if (length(dif)>1){
		for (i in 2:length(dif)){
			if (dif[i]>1) piweight<-c(piweight,rep(0,dif[i]))-c(rep(0,dif[i]-1),-1,piweight)
		}
	}
	piweight<-piweight[1:m]

	atip<-NULL

	num<-NULL
	type<-NULL
	wcoeff<-NULL
	LCrit<-NULL
	if (crit<=0) cat("The Critical value may be positive") 

	va<-mean(residuals^2)

	c<-cumsum(piweight)-1

	d<-rep(0,m)
	delta<-0.7
	d[1]<-piweight[1]-delta
	for (i in 2:m) d[i]<-delta*d[i-1]+piweight[i]

	sum1<-1+sum(piweight*piweight)
	sum2<-1+sum(c*c)
	sum3<-1+sum(d*d)

	maxL<-crit+1

	while (maxL>crit){	
		ka1<-sum1
		ks1<-sum2
		kt1<-sum3
		  maxL<-0
		for (i in 1:m){

			suma1<-sum(residuals[i:m]*c(1,-piweight[1:(m-i)]))
			suma2<-sum(residuals[i:m]*c(1,-c[1:(m-i)]))
			suma3<-sum(residuals[i:m]*c(1,-d[1:(m-i)]))

			ka1<-ka1 - piweight[m-i+1]*piweight[m-i+1]
			w_ao<-suma1/ka1
			v_ao<-va/ka1
			l_ao<-w_ao/sqrt(v_ao)

			ks1<-ks1 - c[m-i+1]*c[m-i+1]
			w_ls<-suma2/ks1
			v_ls<-va/ks1
			l_ls<-w_ls/sqrt(v_ls)

			kt1<-kt1 - d[m-i+1]*d[m-i+1]
			w_tc<-suma3/kt1
			v_tc<-va/kt1
			l_tc<-w_tc/sqrt(v_tc)

			if(abs(l_ao)>maxL){	
				maxL<-abs(l_ao)
				t<-i
				w<-w_ao
				v<-v_ao
				ts<-"AO"
			}
			if(abs(l_ls)>maxL & LS==T){	
				maxL<-abs(l_ls)
				t<-i
				w<-w_ls
				v<-v_ls
				ts<-"LS"
			}
			if(abs(l_tc)>maxL){
				maxL<-abs(l_tc)
				t<-i
				w<-w_tc
				v<-v_tc
				ts<-"TC"
			}
		}

		if(maxL > crit){
			if(ts=="AO") residuals[t:m]<-residuals[t:m]+w*c(-1,piweight[1:(m-t)])
			if(ts=="LS") residuals[t:m]<-residuals[t:m]+w*c(-1,c[1:(m-t)])
			if(ts=="TC") residuals[t:m]<-residuals[t:m]+w*c(-1,d[1:(m-t)])

			val<-mean(residuals^2)
			l<-w/sqrt(v*val/va)
			va<-val

			num<-c(num,t)
			type<-c(type,ts)
			wcoeff<-c(wcoeff,w)
			LCrit<-c(LCrit,abs(l))

			atip<-data.frame(Obs=num,type_detected=type,W_coeff=wcoeff,ABS_L_Ratio=LCrit)
		}
	}
	return(list(atip=atip,sigma2=va,resid=residuals))
}


#funció per linealitzar la serie.

lineal<-function(serie,atip){
	m<-length(serie)
	for(i in 1:nrow(atip)){
		t<-atip[i,1]
		ts<-atip[i,2]
		w<-atip[i,3]
		if(ts=="TC")	serie[t:m]<-serie[t:m]-w*c(1,0.7^(1:(m-t)))
		if(ts=="LS")	serie[t:m]<-serie[t:m]-w
		if(ts=="AO")	serie[t]<-serie[t]-w
	}
	return(serie)
}