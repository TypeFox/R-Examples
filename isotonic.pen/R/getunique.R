getunique <-
function(smat,k){
	n=dim(smat)[1]
	nu=n;group=1:n
	for(i in 1:(n-1)){
		j=i+1
		while(all(smat[i,1:k]==smat[j,1:k])&j<n){
			group[j]=group[i]
			j=j+1
		}
	}
	gru=unique(group)
	nu=length(gru)
	w=1:nu*0+1
	if(nu<n){
		s2mat=matrix(nrow=nu,ncol=k);yu=1:nu
		for(i in 1:nu){
			s2mat[i,]=smat[gru[i],1:k]
			yu[i]=mean(smat[group==gru[i],k+1])
			w[i]=sum(smat[group==gru[i],k+2])
		}
	}else{
		s2mat=smat[,1:k]
		yu=smat[,k+1]
		w=smat[,k+2]
	}
	ans=new.env()
	ans$xmat=s2mat
	ans$y=yu
	ans$group=group
	ans$wt=w
	ans
}
