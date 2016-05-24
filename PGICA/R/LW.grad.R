LW.grad <-
function(w,theta.l,mu.l,sigma.v,ind=500,nproc=10,n,m,X,l.b,u.b,N.s,files){
	
	#	added following lines to avoid notes in package
#	files=get("files");
#	n=get("n");
#	m=get("m");
#	X=get("X");
#	l.b=get("l.b");
#	u.b=get("u.b");
#	N.s=get("N.s");
	
	load(files[1])
	fe=matrix(0,n,m)
	feg=matrix(0,n,m)
	fegg=matrix(0,n,m)
	grad=c()
	gradnn=c()
	res=c()

	w1=t(w[[1]])
	winv=solve(w1)
	x1=t(X)
	mult=pmin(pmax(x1%*%w1,l.b),u.b)
	
	# compute the density values for each column
	for(j in 1:m){
		fapp=fhatnew(mult[,j],theta.l[[j]],mu.l[[j]],sigma.v[j])
		fe[,j]=fapp[,1]
		feg[,j]=fapp[,2]
	}
	
	fegfe=feg/fe	
# 	
# 	for(subj in 1:N.s){
# 		w1=t(w[[subj]])
# 		winv=solve(w1)
# 		load(files[subj])
# 		x1=t(X)
# 		
# 		gradn=matrix(0,m*m,n)	
# 		for(k in 1:m){
# 			for(j in 1:m){
# 				gradn[((k-1)*m+1)+j-1,]=x1[,j]*fegfe[,k]+t(winv)[j,k]
# 			}
# 		}
# 		grad=apply(gradn,1,sum)
# 		
# 		mind=floor(n/ind)
# 		hes=0
# 		for(i in 1:mind){			
# 		hes=hes+gradn[,((i-1)*ind+1):(i*ind)]%*%t(gradn[,((i-1)*ind+1):(i*ind)])}
# 		if((mind*ind)<n){
# 		hes=hes+gradn[,(mind*ind+1):n]%*%t(gradn[,(mind*ind+1):n])}
# 				
# 		res=c(res,list(solve(hes)%*%grad))
# 	}
#	if(par==0){

#  options(sge.trace=TRUE)
#  options(sge.debug=FALSE)
#	options(sge.use.cluster=FALSE);
#  	res=sge.parLapply(1:N.s,parcom,w,files,m,n,fegfe,ind,njobs=30)
#	options(sge.use.cluster=TRUE)

#	}
#	if(par==1){
#	res=sge.parLapply(1:N.s,parcom,w,files,m,n,fegfe,ind,njobs=30)
#	}
#	if(par==2){
#	cl=makeSOCKcluster(c("localhost","localhost"))
#  cl=makeSOCKcluster(rep("localhost",nproc))

  cl = makeCluster(nproc)
	res=parSapply(cl,X=1:N.s,parcom,w,files,m,n,fegfe,ind,X)
  stopCluster(cl)

#  res = sapply(1:N.s,parcom,w,files,m,n,fegfe,ind,X)
#	}		
	return(res)
}
