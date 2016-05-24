mlesol1 <-
function(matches,nl){
	n=dim(matches)[1];m=max(matches)
	sm=1e-5
	if(min(matches)==m){
# everyone has the same answers
		ans=new.env()
		ans$conv=0
		ans$par=1:n*0+1-sm
		ans$val=1e18
		
	}else{
		dguess=1:n
		for(i in 1:n){dguess[i]=mean(matches[,i])/m}
#  find solution
		sm=1e-5
		blower=1:n*0+sm;bupper=1:n*0+1-sm
		anso=optim(dguess,fit,gr=grfit,matches,m,nl,n,method="L-BFGS-B",lower=blower,upper=bupper)
		ans=new.env()
		ans$conv=anso$convergence
		ans$par=anso$par
		ans$val=anso$value
	}
	ans
}
