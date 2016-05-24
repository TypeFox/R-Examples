ccgrouptest <-
function(answermat,group){
	n=dim(answermat)[1];m=dim(answermat)[2];nsim=1000
## check for input errors
	if(length(group)!=n){
		print("ERROR: length of group vector not equal to number of informants in answer matrix")
	}else if(!all(sort(unique(group))[1:2]==c(1,2))|length(sort(unique(group)))!=2){
		print("ERROR: group labels not 1 and 2")
	}else if(length(unique(group))==1){
		print("ERROR: only one group")
	}else{
		answers=unique(sort(round(answermat,2)))
		nl=length(answers)
		if(all(answers==1:nl)){
			xmat=answermat
		}else{
			xmat=matrix(0,nrow=n,ncol=m)
			for(i in 1:nl){xmat[answermat==answers[i]]=i}
		}
#  compute match matrix
		matches=matrix(nrow=n,ncol=n)
		for(i in 1:n){
			for(j in 1:n){
				matches[i,j]=sum(xmat[i,]==xmat[j,])
			}
		}
		print("simulating null distribution of test statistic....")
		ans1=mlesol1(matches,nl)
		dvec1=ans1$par
		ans2a=mlesol1(matches[group==1,group==1],nl)
		dvec2a=ans2a$par
		dvec2=1:n;dvec2[group==1]=dvec2a
		if(sum(group==2)>1){
			ans2b=mlesol1(matches[group==2,group==2],nl)
			dvec2b=ans2b$par
			dvec2[group==2]=dvec2b
		}else if(sum(group==2)==1){
			dvec2b=.99995
			dvec2[group==2]=dvec2b
		}
		ans2=sum(dvec2a)+sum(dvec2b)
		key1=1:m;key2=1:m	
		pans=1:nl;obs4=1:nl
		for(i in 1:m){	
			for(l in 1:nl){pans[l]=sum(dvec2a[xmat[group==1,i]==l])}
			key1[i]=obs4[pans==max(pans)]
			for(l in 1:nl){pans[l]=sum(dvec2b[xmat[group==2,i]==l])}
			key2[i]=obs4[pans==max(pans)]
		}
		if(!all(key1==key2)){
			cdiff=ans2-sum(dvec1)
			key=1:m;pans=1:nl;obs4=1:nl
			pans=1:nl;obs4=1:nl
			for(i in 1:m){
				for(l in 1:nl){pans[l]=sum(dvec1[xmat[,i]==l])}
				key[i]=obs4[pans==max(pans)]
			}
### now find the bootstrap distribution of the lrt
			cdiffb=1:nsim*0
			xboot=matrix(nrow=n,ncol=m)
			for(isim in 1:nsim){
				for(i in 1:n){
					for(k in 1:m){
						u=runif(1)
						if(u<dvec1[i]){xboot[i,k]=key[k]}else{xboot[i,k]=1+trunc(runif(1)*nl)}
					}
				}
#  compute match matrix
				mboot=matrix(nrow=n,ncol=n)
				for(i in 1:n){
					for(j in 1:n){
						mboot[i,j]=sum(xboot[i,]==xboot[j,])
					}
				}
				ans1b=mlesol1(mboot,nl)
				ans2b=fitfun2(mboot,nl,group)
	
				cdiffb[isim]=ans2b-sum(ans1b$par)
			}
			pval=sum(cdiffb>cdiff)/nsim
		}else{pval=1}
	}
	ans=new.env()
	ans$pval=pval
	ans$key1=key1
	ans$key2=key2
	ans$comp1=dvec1
	ans$comp2=dvec2
	ans$simdist=cdiffb
	ans$diff=cdiff
	ans
}
