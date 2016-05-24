ccmle <-
function(answermat){
	n=dim(answermat)[1];m=dim(answermat)[2]
	answers=unique(sort(round(answermat,2)))
	nl=length(answers)
	if(all(answers==1:nl)){
		xmat=answermat
	}else{
		xmat=matrix(0,nrow=n,ncol=m)
		for(i in 1:nl){xmat[answermat==answers[i]]=i}
	}
	matches=matrix(nrow=n,ncol=n)
	for(i in 1:n){
		for(j in 1:n){
			matches[i,j]=sum(xmat[i,]==xmat[j,])
		}
	}
	sm=1e-5
	if(min(matches)==m){
# everyone has the same answers
		ans=new.env()
		ans$conv=0
		ans$par=1:n*0+1
		ans$val=1e18		
	}else{
		dguess=1:n
		for(i in 1:n){dguess[i]=mean(matches[,i])/m}
#  find solution
		sm=1e-5
		blower=1:n*0+sm;bupper=1:n*0+1-sm
		anso=optim(dguess,fit,gr=grfit,matches,m,nl,n,method="L-BFGS-B",lower=blower,upper=bupper)
		key=1:m
		pans=1:nl;obs4=1:nl
		for(i in 1:m){
			for(l in 1:nl){pans[l]=sum(anso$par[xmat[,i]==l])}
			key[i]=obs4[pans==max(pans)]
		}
		ans=new.env()
		ans$conv=anso$convergence
		ans$comp=round(anso$par,4)
		ans$val=anso$value
		ans$key=key
	}
	ans	
}
