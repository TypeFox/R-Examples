PGM.path.neighborhood <-
function(X,Y,nlams=10,startb=0,lambda=NULL)
{
	n = nrow(X); p = ncol(X);
	if(is.null(lambda)){
		lmax = lambdaMax(X)
		lambda = exp(seq(log(lmax),log(0.001*lmax),l=nlams))
	}
	
	#lams = exp(seq(log(lmax),log(.0001),l=nlams));
	#-if(nlams==1){lams = lmax};
	thr = 1e-6; maxit = 1e6;
	Xt = cbind(t(t(rep(1,n))),X);
	L = max(eigen(t(X)%*%X,only.values=TRUE)$values);
	alphas = 0; Bmat = matrix(0,p,nlams);
	if(sum(startb)==0){Bhat = matrix(runif(p+1),p+1,1); Bhat[1] = 0}else{Bhat=startb}
	for(i in 1:nlams){
		iter = 1; ind = 1;
		while(thr<ind & iter<maxit){
			oldb = Bhat;
			tmp = Bhat - (t(Xt)%*%Y - t(Xt)%*%exp(-Xt%*%Bhat))/L;
			#Bhat = matrix(sapply(tmp - lams[i]/L,max,0),p+1,1); Bhat[1] = tmp[1];
			Bhat = matrix(sapply(tmp - lambda[i]/L,max,0),p+1,1); Bhat[1] = tmp[1];
			ind = sum((Bhat - oldb)^2);
			iter = iter + 1;
		}
		alphas[i] = -Bhat[1];
		Bmat[,i] = -Bhat[2:(p+1),drop=FALSE];
	}
  #return(list(alpha=alphas,Bmat=Bmat,lams=lams))
  return(list(alpha=alphas,Bmat=Bmat,lambda=lambda))
}
