ltsreg<-function (x, y, Hsets, h, intercept, intadj, nitermax = NULL) { #Lightly modified from robustbase
    n <- nrow(x)
    p <- ncol(x)
    if (intercept) {
        xint <- cbind(1, x)
        p <- p + 1
    }  else {
        xint <- x
    }
    if (p < 5) {
        eps <- 1e-12
    } else {
        if (p <= 8) {
            eps <- 1e-14
        } else {
            eps <- 1e-16
        }
    }
    data <- cbind(xint, y)
    nsamp <- length(Hsets)
    coeffs <- vector("list", nsamp)
    objfun <- rep(NA, nsamp)
    Osets <- vector("list", nsamp)
    obs_in_set <- Hsets
    delta <- 1000
    objfun1 <- 1000
    niter <- 0
    while (delta > 0.001) {
        niter <- niter + 1
        xin <- xint[obs_in_set, ]
        yin <- y[obs_in_set]
        coeff <- lm.fit(xin, yin)$coefficients
        if (intadj) {
            ti <- y - xint[, 2:p] %*% coeff[2:p]
            coeff[1] <- unimcd(y = ti, h = h, len = n - h + 1)$initmean
        }
        resid <- y - xint %*% coeff
        resid <- abs(resid)
        obs_in_set <- which(resid <= quantile(resid, h/n, type = 1))
        objfun0 <- objfun1
        objfun1 <- mean(resid[obs_in_set]^2)
        delta <- log(objfun0) - log(objfun1)
        if (is.null(nitermax) != 1 && niter == nitermax) 
            delta <- 0
    }
    Osets <- obs_in_set
    objfun <- objfun1
    xin <- xint[obs_in_set, ]
    yin <- y[obs_in_set]
    coeff <- lm.fit(xin, yin)$coefficients
    coeffs <- coeff
    Coeffs <- coeffs
    bestset <- Osets
    out <- list(coef = Coeffs, best = bestset, allc = coeffs, 
        alls = Osets, objfct = objfun)
    return(out)
}
UniMCD_test<-function(y,h){
#test function for unimcd_in.cpp
    out<-list()
    quan<-h
    ncas<-length(y)
    len<-ncas-quan+1
    if(len==1){
        out$tmcd<-mean(y)
        out$smcd<-sqrt(var(y))
    }else{
        ay<-c()
        I<-order(y)
        y<-y[I]
        ay[1]<-sum(y[1:quan])
        for(samp in 2:len){
            ay[samp]<-ay[samp-1]-y[samp-1]+y[samp+quan-1]
        }
        ay2<-ay^2/quan
        sq<-c()
        sq[1]<-sum(y[1:quan]^2)-ay2[1]
        for(samp in 2:len){
            sq[samp]<-sq[samp-1]-y[samp-1]^2+y[samp+quan-1]^2-ay2[samp]+ay2[samp-1]
        }
        sqmin<-min(sq)
        Isq<-order(sq)
        ndup<-sum(sq==sqmin)
        ii<-Isq[1:ndup]
        slutn<-c()
        slutn[1:ndup]<-ay[ii]
        out$initmean<-slutn[floor((ndup+1)/2)]/quan
        out$initcov<-sqmin/(quan-1)
    }
    return(out)
}
scaleTau2_test<-function(x,c1=4.5,c2=3,consistency=TRUE,mu.too=FALSE){
#test function for scaleTau2.cpp
    n<-length(x)
	tol<-1e-8
    medx<-median(x)
    x.<-abs(x-medx)
    sigma0<-median(x.)
	if(sigma0<tol)	print(sigma0)
    mu<-if(c1>0){
        x.<-x./(sigma0*c1)
        w<-1-x.*x.
        w<-((abs(w)+w)/2)^2
        sum(x*w)/sum(w)
    }
    else medx
    x<-(x-mu)/sigma0
    rho<-x^2
    rho[rho>c2^2]<-c2^2
    if (!identical(consistency,FALSE)){
        Erho<-function(b) 2*((1-b^2)*pnorm(b)-b*dnorm(b)+b^2)-1
        Es2<-function(c2) Erho(c2*qnorm(3/4))
        nEs2<-(if (consistency=="finiteSample") 
            n-2
        else n)*Es2(c2)
    }
    else nEs2<-n 
	if(sum(rho)<1e-8)	stop('singular data.')
  c(if (mu.too) mu,sigma0*sqrt(sum(rho)/nEs2))
}
test_fxOGK<-function(x0,y0,cent_est='mean',scal_est='sd',alpha=0.5){#cent_est<-"mean"
#test function for FFOgkBasis.cpp
	x2<-cbind(x0,y0)
	n<-nrow(x0)
	p<-ncol(x2);
	U<-diag(0,p)
	tol<-1e-8
	h<-DetR::quanf(alpha,n=n,p=p)
	fun1<-get(cent_est)
	fun2<-get(scal_est)
	if(scal_est==cent_est){
		F<-apply(x2,2,fun1,mu.too=TRUE)
	} else {
		F<-rbind(apply(x2,2,fun1),apply(x2,2,fun2))
	}
	if(min(F[2,])<tol)	stop('singular data.')
	x2<-scale(x2,center=F[1,],scale=F[2,])
	for(i in 2:p){#i<-1
		sYi=x2[,i];
		for(j in 1:(i-1)){#j<-1
			sYj=x2[,j]
			sY=sYi+sYj
			dY=sYi-sYj
			U[i,j]=0.25*(fun2(sY)^2-fun2(dY)^2);
			U[j,i]=U[i,j]
		}
	}
	diag(U)<-1
	d1<-svd(U)
	x3<-x2%*%d1$v
	d4<-apply(x3,2,fun2)
	if(min(d4)<tol)	stop('singular data.')
	x3<-sweep(x3,2,d4,FUN='/')
	d5<-rowSums(x3*x3)
	ois<-which(d5<=quantile(d5,h/n));
	ois<-test_cov_Cstep(ois,x2,h)
	list(bestRaw=ois$H0,U=U)
}
test_cov_Cstep<-function(ois,x2,h){
	old_obj<-Inf
	cont<-1
	nstep<-0
	while(cont){
		nstep<-nstep+1
		a1<-mahalanobis(x2,colMeans(x2[ois,]),cov(x2[ois,]))
		ois<-which(a1<=quantile(a1,h/nrow(x2)))
		new_obj<-det(cov(x2[ois,]),log=TRUE)
		if(old_obj-new_obj>1e-3){
			old_obj<-new_obj
		} else {
			cont<-0
		}
	}
	return(list(H0=ois,nstep=nstep))
}
test_fxOGK2<-function(x0,y0,cent_est='mean',scal_est='sd'){
	x2<-cbind(x0,y0)
	p<-ncol(x2);
	n<-nrow(x0)
	U<-diag(0,p)
	tol<-1e-8
	if(is.null(h))	h<-ceiling((n+p)/2)
	fun1<-get(cent_est)
	fun2<-get(scal_est)
	F<-rbind(apply(x2,2,fun1),apply(x2,2,fun2))
	if(min(F[,2])<tol)	stop('singular data.')
	x2<-scale(x2,center=F[1,],scale=F[2,])
	lm(x2[,p]~x2[,-p]-1)$coef
}
test_Cstep<-function(x,y,h,z0){
#test function for Cstep.cpp
	n<-nrow(x)
	w0<-abs(lm(y~x,weights=as.numeric((1:n)%in%z0))$resid);
	d1<-Inf
	carry_on<-1
	citer<-0
	while(carry_on){
		z0<-which(w0<=quantile(w0,h/n))
		cf<-lm(y[z0]~x[z0,])$coef	
		w0<-y-x%*%cf[-1]
		cf[1]<-UniMCD_test(y=w0,h=h)$initmean
		d0<-d1
		citer<-citer+1
		w0<-abs(y-cbind(1,x)%*%cf)
		d1<-mean(w0[z0])
		if(log(d0)-log(d1)<1e-3)	carry_on<-0
	}
	list(bestCStep=z0,citer=citer,objfun=d1,old_objfun=d0)
}
test_function<-function(){
#Test functions for the cpp code
}
