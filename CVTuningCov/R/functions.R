
#Covariance Matrix with AR1 structure
AR1<-function(p,rho=0.5)
{
	Sigma<-matrix(0,p,p);
	for (i in 1:p)
	{
		for (j in 1:p)
		{
			Sigma[i,j]<-rho^(abs(i-j));
		}
	}
	return(Sigma);
}

#Calculate squared L2 norm
L2.norm2<-function(A)
{
	norm<-( eigen( (A)%*%t(A) )$values )[1];
	return(norm);
}

#Calculate squared F norm
F.norm2<-function(A)
{
	norm<-sum(A^2);
	return(norm);
}


#define tapering matrix
tapering<-function(p,k=1)
{
	k<-min(k,2*(p-1));
	w<-matrix(0,p,p);
	for (i in 1:p)
	{
		for(j in 1:p)
		{
			if(abs(i-j)<=k/2){ w[i,j]<-1;}
			else 
			if (k/2 < abs(i-j) & abs(i-j)<k){w[i,j]<-2 - 2*abs(i-j)/k;}
		}
	}
	return(w);
}

#define banding matrix
banding<-function(p,k=1)
{
	B<-matrix(0,p,p);
	for (i in 1:p)
	{
		for(j in 1:p)
		{
			B[i,j]<-1*(abs(i-j)<=k);
		}
	}
	return(B);
}

#hard thresholding
hard.thresholding<-function(Sigma,c=0.5)
{
	threshold.Sigma<-Sigma*(abs(Sigma)>=c)
	return(threshold.Sigma);
}

#soft thresholding
soft.thresholding<-function(Sigma,c=0.5)
{
	threshold.Sigma<-sign(Sigma)*pmax(abs(Sigma)-c,0);
	return(threshold.Sigma);
}


#CV tuning
regular.CV<-function(X,k.grid=0.5, method='Tapering',fold=5,norm='F',seed=10323)
{
	n<-nrow(X);
	p<-ncol(X);
	length.k.grid<-length(k.grid);
	#sample variance matrix
	tilde.Sigma<-var(X);

	#generate CV index
	set.seed(seed);
	CV.index<-sample(1:n,size=n, replace = FALSE)

	CV.pre.error<-matrix(0,length.k.grid,2);
	fold<-round(fold);
	fold.size<-round(n/fold);

	#evaluate all k values
	for (j in 1:length.k.grid)
	{
		k<-k.grid[j];
		if (method=='Banding'){
			W<-banding(p,k);
			Sigma.k<-W*tilde.Sigma;
		}else if (method=='HardThresholding'){
			Sigma.k<-hard.thresholding(tilde.Sigma,k);
		}else if (method=='SoftThresholding'){
			Sigma.k<-soft.thresholding(tilde.Sigma,k);
		}else {
			W<-tapering(p,k);
			Sigma.k<-W*tilde.Sigma;
		}
		

		#CV 
		CV<-rep(0,2); #1: regular CV,larger training set; 2:swith training and test
		for (cv.v in 1:fold)
		{
			if (cv.v==fold){
				cv.test<- CV.index[ ((cv.v-1)*fold.size+1):n];
				cv.train<-cv.test; #switch training and testing;
			}else{
				cv.test<-CV.index[ ((cv.v-1)*fold.size+1):(cv.v*fold.size)];
				cv.train<-cv.test; #switch training and testing;
			}
			if (norm=='L2'){
				if (method=='HardThresholding')
				{
					CV[1]<-CV[1]+L2.norm2( hard.thresholding(var(X[-(cv.test),]),k) - var(X[cv.test,]) );
					CV[2]<-CV[2]+L2.norm2( hard.thresholding(var(X[cv.train,]),k) - var(X[-cv.train,]) );
				}else	if (method=='SoftThresholding')
				{
					CV[1]<-CV[1]+L2.norm2( soft.thresholding(var(X[-(cv.test),]),k) - var(X[cv.test,]) );
					CV[2]<-CV[2]+L2.norm2( soft.thresholding(var(X[cv.train,]),k) - var(X[-cv.train,]) );
				}else{	
					CV[1]<-CV[1]+L2.norm2( W*var(X[-(cv.test),]) - var(X[cv.test,]) );
					CV[2]<-CV[2]+L2.norm2( W*var(X[cv.train,]) - var(X[-cv.train,]) );
				}
			}else{
				if (method=='HardThresholding')
				{
					CV[1]<-CV[1] + F.norm2( hard.thresholding(var(X[-(cv.test),]),k) - var(X[cv.test,]) );
					CV[2]<-CV[2] + F.norm2( hard.thresholding(var(X[cv.train,]),k) - var(X[-cv.train,]) );
				}else	if (method=='SoftThresholding')
				{
					CV[1]<-CV[1] + F.norm2( soft.thresholding(var(X[-(cv.test),]),k) - var(X[cv.test,]) );
					CV[2]<-CV[2] + F.norm2( soft.thresholding(var(X[cv.train,]),k) - var(X[-cv.train,]) );
				}else{	
					CV[1]<-CV[1]+ F.norm2( W*var(X[-(cv.test),]) - var(X[cv.test,]) );
					CV[2]<-CV[2]+ F.norm2( W*var(X[cv.train,]) - var(X[-cv.train,]) );
				}

			}	
		}
		CV.pre.error[j,]<-CV/fold;
	}

	CV.k<-rep(NA,2);
	CV.k[1]<-k.grid[which.min(CV.pre.error[,1])];
	CV.k[2]<-k.grid[which.min(CV.pre.error[,2])];

	return(list(CV.k=CV.k, k.grid=k.grid, CV.pre.error=CV.pre.error) );
	
}


#random CV tuning
random.CV<-function(X,k.grid=0.5, method='Tapering',test.size=5,norm='F',boot.num=50,seed=10323)
{
	n<-nrow(X);
	p<-ncol(X);
	length.k.grid<-length(k.grid);
	#sample variance matrix
	tilde.Sigma<-var(X);

	CV.pre.error<-matrix(0,length.k.grid,2);
	test.size<-round(test.size);

	#evaluate all k values
	for (j in 1:length.k.grid)
	{
		k<-k.grid[j];
		if (method=='Banding'){
			W<-banding(p,k);
			Sigma.k<-W*tilde.Sigma;
		}else if (method=='HardThresholding'){
			Sigma.k<-hard.thresholding(tilde.Sigma,k);
		}else if (method=='SoftThresholding'){
			Sigma.k<-soft.thresholding(tilde.Sigma,k);
		}else {
			W<-tapering(p,k);
			Sigma.k<-W*tilde.Sigma;
		}
		

		#random CV with boot.num repetitions 
		CV<-rep(0,2); #1: regular CV,larger training set; 2:swith training and test
		
		for (cv.v in 1:boot.num)
		{

			#generate random CV index
			set.seed(cv.v*seed);
			cv.test<-sample(1:n,size=test.size, replace = FALSE);
			cv.train<-cv.test; #switch training and testing;

			if (norm=='L2'){
				if (method=='HardThresholding')
				{
					CV[1]<-CV[1]+L2.norm2( hard.thresholding(var(X[-(cv.test),]),k) - var(X[cv.test,]) );
				}else if (method=='SoftThresholding')
				{
					CV[1]<-CV[1]+L2.norm2( soft.thresholding(var(X[-(cv.test),]),k) - var(X[cv.test,]) );
				}else{	
					CV[1]<-CV[1]+L2.norm2( W*var(X[-(cv.test),]) - var(X[cv.test,]) );
				}
			}else{
				if (method=='HardThresholding')
				{
					CV[1]<-CV[1] + F.norm2( hard.thresholding(var(X[-(cv.test),]),k) - var(X[cv.test,]) );
				}else if (method=='SoftThresholding')
				{
					CV[1]<-CV[1] + F.norm2( soft.thresholding(var(X[-(cv.test),]),k) - var(X[cv.test,]) );
				}else{	
					CV[1]<-CV[1]+ F.norm2( W*var(X[-(cv.test),]) - var(X[cv.test,]) );
				}

			}	
		}
		CV.pre.error[j,]<-CV/boot.num;
	}

	CV.k<-rep(0,2);
	CV.k[1]<-k.grid[which.min(CV.pre.error[,1])];

	return(list(CV.k=CV.k, k.grid=k.grid, CV.pre.error=CV.pre.error) );

}



