## updated Sept 11th
## updated May 25th. Am. Nat. April 1997; Martins and Hansen
## updated last on June 15th. rewrite the ouput for new prediction and consequently for the analysis with no unknown species
## updated on August 31st. Incorporate more unknown entries in the data file for prediction (>1 data points). And solve the issue of dimension error with the position of the nodes.
##correct the degree of freedom for sigma-sq



pGLS<-function(formula,data,covarmatrix,na.action, intercept = TRUE) {
	cl <- match.call()
   	mf <- match.call(expand.dots = FALSE)
	
   	m <- match(c("formula", "data", "covarmatrix", "na.action"), names(mf), 0)

	mf <- mf[c(1, m)]
	
	mf$drop.unused.levels = TRUE
	
	mf[[1]] <- as.name("model.frame")
	
	mf <- eval(mf)
	
	
	mt <- attr(mf, "terms")
	
	y <- model.response(mf)
	

	
	x <- model.matrix(mt, mf)
	
	if(!intercept) {
		x=x[,-1]
	}
	# verify and compute
	library("MASS")
	covarmatrix_t<-covarmatrix
	y_t<-y
	x_t<-x
	
	for(i in 1:length(x[,1]) ) {
		for(j in 1:length(x[1,])) {
			if(is.na(x[i,j])) {
				stop("the predictors must not have missing!")
			}
		}
	}
	


	n_na=0;
	nalist=NA;
	for(i in length(covarmatrix[,1]):1) {
		
		if(is.na(y[i]) ) {
			n_na=n_na+1;
			nalist=rbind(nalist,i);
			if(n_na<(length(y_t)-1)) {
				covarmatrix=covarmatrix[,-i];
				covarmatrix=covarmatrix[-i,];
				y=y[-i]
				x=x[-i,]
				
			} else {
				stop("2 many NA in the dataset!")

			}
		}
		
	}
	nalist=nalist[-1];
	n=length(x[,1])
	p=length(x[1,])

	if(length(x[,1])<=(length(x[1,])+1)){ 
		stop("2 few datapoints!")
	}
	
	svdm<-svd(covarmatrix)
	# $u $d(diag) $v
	dmar=svdm$u %*% diag(sqrt(1/svdm$d)) %*% t(svdm$v)
	
	#predictor
	zmar=dmar %*% y
	
	umar=dmar %*% x
	bhat=ginv(t(umar)%*% umar) %*% (t(umar)%*% zmar)
	sigma_sq=t(zmar-umar%*%bhat)%*%(zmar-umar%*%bhat)/(n-p)
	varb=as.vector(sigma_sq) * (ginv(t(umar)%*%umar))
	yhat= x %*% bhat
	vary=x%*%varb%*%t(x)	
	svdm_t<-svd(covarmatrix_t)
	dmar_t=svdm_t$u %*% diag(sqrt(1/svdm_t$d)) %*% t(svdm_t$v)
	umar_t=dmar_t %*% x_t
	varb_t=as.vector(sigma_sq) * (ginv(t(umar_t)%*%umar_t))
	yhat_t=ginv(dmar_t) %*% (dmar_t %*% x_t %*% bhat)
	vary_t=x_t%*%varb_t%*%t(x_t)	
	ones=c(1:length(x[,1]))
	for( element in 1:length(x[,1])) {
		ones[element]=1;
	}
	w_avgy=(t(ones)%*%ginv(covarmatrix)%*%y)/(t(ones)%*%ginv(covarmatrix)%*%ones)
	a_t=c(1:length(x[,1]))
	for( element in 1:length(x[,1])) {
		a_t[element]=w_avgy;
	}
	print(yhat-y)
	print(covarmatrix)
	R_2=1-(t(as.vector(yhat-y))%*%ginv(covarmatrix)%*%as.vector(yhat-y))/(t(as.vector(a_t-y))%*%ginv(covarmatrix)%*%as.vector(a_t-y));
	if(n_na>0) {
		y_output=data.frame(cbind(1:length(x_t[,1]),1:length(x_t[,1]),1:length(x_t[,1])))
		colnames(y_output)=c("species","predicted","s(predicted)")
		for(i in 1:length(x_t[,1])) {

			mu1=yhat_t[i]

			sig1=vary_t[i,i]

			y_output[i,1]=as.character(data[i,1])

			y_output[i,2]=mu1;
			y_output[i,3]=sqrt(sig1);


		}
		t_score=(bhat)/sqrt(diag(varb))
		z=abs(t_score)
		
		t_pro=(1-pt(z,n-p))*2
		coeff=cbind(bhat,sqrt(diag(varb)),t_score,t_pro)
		
		
		colnames(coeff)=c("coefficients","std.err","T_score","p_value")
		rownames(coeff)=colnames(x)

		#conditional distribution of x|x_others; work on the conditionl dist. one datapoint at a time

		pred.new=data.frame(cbind(1,1,1:n_na));

		for(i in 1:n_na) {
			y_p=y_t;
			x_p=x_t;
			covarmatrix_p<-covarmatrix_t
			for(j in n_na:1) {
				if(i!=j) {
					covarmatrix_p=covarmatrix_p[,-nalist[j]]
					covarmatrix_p=covarmatrix_p[-nalist[j],]
					x_p=x_p[-nalist[j]]
					y_p=y_p[-nalist[j]]
				}
			}
			pred_no=0;
			for(k in length(y_p):1) {
		
				if(is.na(y_p[k]) ) {
					pred_no=k;
				}
						
						
			}
			s12=as.vector(covarmatrix_p[pred_no,])
			s12=s12[-pred_no]
			s22=covarmatrix_p[-pred_no,]
			print(pred_no)
			s22=s22[,-pred_no]
			
		
			mu=yhat_t[nalist[i]]+t(s12) %*% ginv(s22) %*% (yhat-y)
			sig=(covarmatrix_p[pred_no,pred_no]- t(s12) %*% ginv(s22) %*% s12)*sigma_sq
			pred.new1=cbind(as.character(data[nalist[i],1]),mu,sig)
			pred.new[i,]=pred.new1[1,]
		}
		
		colnames(pred.new)=c("species","predicted_value(cond.)","Std.Error(cond.)")
		z <- list(pred=y_output,coefficients=coeff,cov.coeff=varb,"sigma^2"=sigma_sq,"pred.cond."=pred.new,"R-Sq"=R_2)
		
	}else {
		y_output=data.frame(cbind(1:length(x[,1]),1:length(x[,1]),1:length(x[,1])))
		colnames(y_output)=c("species","predicted","s(predicted)")
		for(i in 1:length(x[,1])) {

			mu1=yhat[i]

			sig1=vary[i,i]

			y_output[i,1]=as.character(data[i,1])

			y_output[i,2]=mu1;
			y_output[i,3]=sqrt(sig1);


		}
		t_score=(bhat)/sqrt(diag(varb))
		z=abs(t_score)
		t_pro=(1-pt(z,n-p))*2
		coeff=cbind(bhat,sqrt(diag(varb)),t_score,t_pro)	
		colnames(coeff)=c("coefficients","std.err","T_score","p_value")
		rownames(coeff)=colnames(x)	
		z <- list(pred=y_output,coefficients=coeff,cov.coeff=varb,"sigma^2"=sigma_sq,"R-Sq"=R_2)
	}
	print(z)
}



